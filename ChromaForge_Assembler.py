# ==============================================================================
# Script Name: ChromaForge - Automated Sanger Sequencing Assembler (Beta)
# Author: Muhammad Aldiansyah Karisma
# Affiliation: Universitas Bengkulu
# Year: 2026
# License: MIT License
# Description: Automated script for resolving contig conflicts using Phred Q-scores.
#              Outputs dynamic FASTA files based on IUPAC ambiguity detection.
# ==============================================================================

import os
import glob
import re
import time
from difflib import SequenceMatcher
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner
from Bio.Blast import NCBIWWW, NCBIXML

if os.name == 'nt':
    os.system('')

WARNA_BIRU = '\033[94m'
WARNA_KUNING = '\033[93m'
WARNA_ORANYE = '\033[38;5;208m'
WARNA_MERAH = '\033[91m'
WARNA_UNGU = '\033[95m'
RESET_WARNA = '\033[0m'

def dapatkan_kode_iupac(basa1, basa2):
    kombinasi = set([basa1.upper(), basa2.upper()])
    if kombinasi == {'A', 'G'}: return 'R'
    if kombinasi == {'C', 'T'}: return 'Y'
    if kombinasi == {'G', 'C'}: return 'S'
    if kombinasi == {'A', 'T'}: return 'W'
    if kombinasi == {'G', 'T'}: return 'K'
    if kombinasi == {'A', 'C'}: return 'M'
    return 'N'

def robust_trim_masking_dan_evaluasi(record, window_size=5, threshold_q=20):
    quals = record.letter_annotations.get("phred_quality")
    seq_str = str(record.seq)
    
    if not quals:
        return seq_str, 0.0, []
        
    length = len(quals)
    start_idx = 0
    for i in range(length - window_size + 1):
        window = quals[i:i + window_size]
        if sum(window) / window_size >= threshold_q:
            start_idx = i
            break
            
    end_idx = length
    for i in range(length - window_size, -1, -1):
        window = quals[i:i + window_size]
        if sum(window) / window_size >= threshold_q:
            end_idx = i + window_size
            break
            
    if start_idx >= end_idx:
        return seq_str, 0.0, []
        
    trimmed_seq = seq_str[start_idx:end_idx]
    trimmed_quals = quals[start_idx:end_idx]
    trim_length = len(trimmed_seq)
    
    if trim_length == 0:
        return seq_str, 0.0, []
        
    masked_chars = []
    good_count = 0
    for base, q in zip(trimmed_seq, trimmed_quals):
        if q < threshold_q:
            masked_chars.append('N')
        else:
            masked_chars.append(base)
            good_count += 1
            
    processed_seq = "".join(masked_chars)
    ratio_good = good_count / trim_length
    
    return processed_seq, ratio_good, trimmed_quals

def bangun_strict_consensus(f_seq, r_seq, alignments, f_qual, r_qual):
    alignment = alignments[0]
    t_span = alignment.aligned[0]
    q_span = alignment.aligned[1]
    
    if len(t_span) == 0 or len(q_span) == 0:
        return None
        
    f_start_overlap = t_span[0][0]
    consensus_iupac = list(f_seq[:f_start_overlap])
    consensus_pseudo = list(f_seq[:f_start_overlap])
    
    snp_count = 0
    log_mismatch = []
    
    for t_block, q_block in zip(t_span, q_span):
        for i in range(t_block[1] - t_block[0]):
            f_idx = t_block[0] + i
            q_idx = q_block[0] + i
            
            f_base = f_seq[f_idx]
            r_base = r_seq[q_idx]
            
            posisi_contig = len(consensus_iupac) + 1
            
            if f_base == r_base:
                consensus_iupac.append(f_base)
                consensus_pseudo.append(f_base)
            else:
                q_f = f_qual[f_idx] if f_qual else 0
                q_r = r_qual[q_idx] if r_qual else 0
                
                if q_f > q_r + 10:
                    consensus_iupac.append(f_base)
                    consensus_pseudo.append(f_base)
                    log_mismatch.append(f"Contig[{posisi_contig}]:F[{f_idx}]/Rc[{q_idx}]")
                elif q_r > q_f + 10:
                    consensus_iupac.append(r_base)
                    consensus_pseudo.append(r_base)
                    log_mismatch.append(f"Contig[{posisi_contig}]:F[{f_idx}]/Rc[{q_idx}]")
                else:
                    if q_f >= 30 and q_r >= 30:
                        kode_iupac = dapatkan_kode_iupac(f_base, r_base)
                        consensus_iupac.append(kode_iupac)
                        basa_terpilih = f_base if q_f >= q_r else r_base
                        consensus_pseudo.append(basa_terpilih)
                        snp_count += 1
                        log_mismatch.append(f"Contig[{posisi_contig}]:F[{f_idx}]/Rc[{q_idx}]({kode_iupac} -> PseudoPhased:{basa_terpilih})")
                    else:
                        consensus_iupac.append('N')
                        consensus_pseudo.append('N')
                        snp_count += 1
                        log_mismatch.append(f"Contig[{posisi_contig}]:F[{f_idx}]/Rc[{q_idx}](N)")
                    
    r_end_overlap = q_span[-1][1]
    consensus_iupac.extend(list(r_seq[r_end_overlap:]))
    consensus_pseudo.extend(list(r_seq[r_end_overlap:]))
    
    return "".join(consensus_iupac), "".join(consensus_pseudo), snp_count, log_mismatch

def proses_sampel_dengan_aturan_lengkap(f_rec, r_rec, id_sampel, min_length=200):
    ambang_layak = 0.85 
    ambang_minimum = 0.70
    laporan_teks = None
    
    if f_rec and not r_rec:
        f_proc, f_score, _ = robust_trim_masking_dan_evaluasi(f_rec)
        if len(f_proc) < min_length:
            print(f"    -> {WARNA_MERAH}[SKIPPED] Sekuen F tunggal terlalu pendek ({len(f_proc)} bp).{RESET_WARNA}")
            return None, None, "SKIPPED", None
        print(f"    -> {WARNA_KUNING}[Info] Hanya F tersedia. Diproses sebagai sekuen tunggal (Q: {f_score*100:.1f}%, Len: {len(f_proc)} bp).{RESET_WARNA}")
        return Seq(f_proc), Seq(f_proc), f"Single Read F (Q: {f_score*100:.1f}%)", None

    if r_rec and not f_rec:
        r_proc, r_score, _ = robust_trim_masking_dan_evaluasi(r_rec)
        if len(r_proc) < min_length:
            print(f"    -> {WARNA_MERAH}[SKIPPED] Sekuen R tunggal terlalu pendek ({len(r_proc)} bp).{RESET_WARNA}")
            return None, None, "SKIPPED", None
        print(f"    -> {WARNA_KUNING}[Info] Hanya R tersedia. Diputar menjadi R-rc sebagai sekuen tunggal (Q: {r_score*100:.1f}%, Len: {len(r_proc)} bp).{RESET_WARNA}")
        seq_r_rc = Seq(r_proc).reverse_complement()
        return seq_r_rc, seq_r_rc, f"Single Read R-rc (Q: {r_score*100:.1f}%)", None
        
    f_proc, f_score, f_qual = robust_trim_masking_dan_evaluasi(f_rec)
    r_proc, r_score, r_qual_raw = robust_trim_masking_dan_evaluasi(r_rec)
    
    f_layak = f_score >= ambang_layak
    r_layak = r_score >= ambang_layak
    
    if f_layak and r_layak:
        r_rc = str(Seq(r_proc).reverse_complement())
        r_qual_rc = r_qual_raw[::-1] if r_qual_raw else []
        
        aligner = PairwiseAligner()
        aligner.mode = 'local'
        aligner.match_score = 2
        aligner.mismatch_score = -3
        aligner.open_gap_score = -5
        aligner.extend_gap_score = -2
        
        alignments = aligner.align(f_proc, r_rc)
        if not alignments:
            if f_score >= r_score:
                if len(f_proc) < min_length:
                    print(f"    -> {WARNA_MERAH}[SKIPPED] Contig gagal, fallback F terlalu pendek ({len(f_proc)} bp).{RESET_WARNA}")
                    return None, None, "SKIPPED", None
                return Seq(f_proc), Seq(f_proc), f"Contig gagal, fallback F (Q: {f_score*100:.1f}%)", None
            else:
                if len(r_proc) < min_length:
                    print(f"    -> {WARNA_MERAH}[SKIPPED] Contig gagal, fallback R-rc terlalu pendek ({len(r_proc)} bp).{RESET_WARNA}")
                    return None, None, "SKIPPED", None
                seq_r_rc = Seq(r_rc)
                return seq_r_rc, seq_r_rc, f"Contig gagal, fallback R-rc (Q: {r_score*100:.1f}%)", None
                
        hasil_consensus = bangun_strict_consensus(f_proc, r_rc, alignments, f_qual, r_qual_rc)
        
        if hasil_consensus:
            seq_final_iupac, seq_final_pseudo, snp_total, log_mismatch = hasil_consensus
            
            if len(seq_final_iupac) < min_length:
                print(f"    -> {WARNA_MERAH}[SKIPPED] Contig terlalu pendek ({len(seq_final_iupac)} bp). Cek manual.{RESET_WARNA}")
                return None, None, "SKIPPED", None
                
            if log_mismatch:
                mismatch_str = ", ".join(log_mismatch)
                status = f"Contig Strict (Mismatch: {len(log_mismatch)} titik -> Posisi: {mismatch_str})"
                laporan_teks = f"SAMPEL: {id_sampel}\nPOSISI: {mismatch_str}\n\n{str(alignments[0])}\n\n"
            else:
                status = "Contig Strict (Identik Sempurna)"
                
            if str(seq_final_iupac) != str(seq_final_pseudo):
                print(f"    -> {WARNA_UNGU}[IUPAC] {status} (Len: {len(seq_final_iupac)} bp){RESET_WARNA}")
            else:
                print(f"    -> {WARNA_BIRU}[OK] {status} (Len: {len(seq_final_iupac)} bp){RESET_WARNA}")
                
            return Seq(seq_final_iupac), Seq(seq_final_pseudo), status, laporan_teks
            
        return Seq(f_proc), Seq(f_proc), f"Contig fallback F", None
    elif f_layak:
        if len(f_proc) < min_length:
            print(f"    -> {WARNA_MERAH}[SKIPPED] Hanya F layak namun terlalu pendek ({len(f_proc)} bp).{RESET_WARNA}")
            return None, None, "SKIPPED", None
        print(f"    -> {WARNA_KUNING}[Info] Hanya F layak (Q: {f_score*100:.1f}%, Len: {len(f_proc)} bp).{RESET_WARNA}")
        return Seq(f_proc), Seq(f_proc), f"Hanya F (F layak: {f_score*100:.1f}%)", None
    elif r_layak:
        if len(r_proc) < min_length:
            print(f"    -> {WARNA_MERAH}[SKIPPED] Hanya R layak namun terlalu pendek ({len(r_proc)} bp).{RESET_WARNA}")
            return None, None, "SKIPPED", None
        print(f"    -> {WARNA_KUNING}[Info] Hanya R layak (Q: {r_score*100:.1f}%, Len: {len(r_proc)} bp).{RESET_WARNA}")
        seq_r_rc = Seq(r_proc).reverse_complement()
        return seq_r_rc, seq_r_rc, f"Hanya R-rc (R layak: {r_score*100:.1f}%)", None
    else:
        if f_score >= r_score:
            seq_recovery = f_proc
            tipe_bacaan = "F"
            seq_final = Seq(seq_recovery)
        else:
            seq_recovery = r_proc
            tipe_bacaan = "R-rc"
            seq_final = Seq(seq_recovery).reverse_complement()
            
        if len(seq_recovery) < min_length:
            print(f"    -> {WARNA_MERAH}[SKIPPED] Kualitas sangat rendah. Sisa {tipe_bacaan} terlalu pendek ({len(seq_recovery)} bp).{RESET_WARNA}")
            return None, None, "SKIPPED", None
            
        status_metode = f"Best Single-Read ({tipe_bacaan})"
        print(f"    -> {WARNA_ORANYE}[Warning] Kualitas di bawah standar. Menggunakan {status_metode} (Len: {len(seq_recovery)} bp).{RESET_WARNA}")
        return seq_final, seq_final, status_metode, None

def ekstrak_id_cerdas(clean_name):
    core = re.sub(r'^\d+_', '', clean_name)
    parts = core.split('_')
    if len(parts) > 1 and len(parts[0]) <= 2: 
        return f"{parts[0]}_{parts[1]}"
    return parts[0]

def pasangkan_sampel_otomatis_ai(berkas_ab1):
    f_files = []
    r_files = []
    
    for path in berkas_ab1:
        nm_file = os.path.basename(path)
        nm_for_detect = re.sub(r'_(resend|repeat)', '', nm_file, flags=re.IGNORECASE).upper()
        
        if any(k in nm_for_detect for k in ['_R.', '_R_', 'BR1', 'GADPHR', 'CYTBR', '16SD_R', 'R_G', 'R1.', 'R1_']) or nm_for_detect.replace('.AB1','').endswith('R'):
            arah = 'R'
        else:
            arah = 'F'
            
        clean_name = re.sub(r'\.ab1$', '', nm_file, flags=re.IGNORECASE)
        clean_name = re.sub(r'^1st_BASE_\d+_', '', clean_name, flags=re.IGNORECASE)
        
        if arah == 'F':
            f_files.append((clean_name, path))
        else:
            r_files.append((clean_name, path))
            
    pasangan_sampel = {}
    
    for f_clean, f_path in f_files:
        best_match_r = None
        best_score = 0
        
        for r_clean, r_path in r_files:
            score = SequenceMatcher(None, f_clean.lower(), r_clean.lower()).ratio()
            if score > best_score:
                best_score = score
                best_match_r = r_path
                
        if best_match_r and best_score > 0.45: 
            id_sampel = ekstrak_id_cerdas(f_clean)
            pasangan_sampel[id_sampel] = {'F': f_path, 'R': best_match_r}
            r_files = [(c, p) for c, p in r_files if p != best_match_r]
        else:
            id_sampel = ekstrak_id_cerdas(f_clean)
            pasangan_sampel[id_sampel] = {'F': f_path, 'R': None}
            
    for r_clean, r_path in r_files:
        id_sampel = ekstrak_id_cerdas(r_clean)
        pasangan_sampel[id_sampel] = {'F': None, 'R': r_path}
        
    return pasangan_sampel

def urutkan_natural(teks):
    return [int(c) if c.isdigit() else c.lower() for c in re.split(r'(\d+)', teks)]

def koreksi_orientasi_via_blastn(seq_record):
    print(f"    -> [BLAST] Mengirim {seq_record.id} ke NCBI BLASTn... (Jeda aman 20 detik aktif, harap tunggu)")
    try:
        result_handle = NCBIWWW.qblast("blastn", "nt", seq_record.seq, hitlist_size=1)
        blast_record = NCBIXML.read(result_handle)
        
        if blast_record.alignments:
            alignment = blast_record.alignments[0]
            hsp = alignment.hsps[0]
            
            if hsp.frame[1] < 0:
                print(f"    -> [BLAST] Orientasi Minus terdeteksi (HSP Frame: {hsp.frame}). Memutar orientasi...")
                rev_seq = seq_record.seq.reverse_complement()
                return SeqRecord(rev_seq, id=seq_record.id, description=seq_record.description + " | BLAST_Orientation: Reversed")
            else:
                print(f"    -> [BLAST] Orientasi Plus/Plus sudah akurat.")
                return SeqRecord(seq_record.seq, id=seq_record.id, description=seq_record.description + " | BLAST_Orientation: OK")
        else:
            print(f"    -> [BLAST] Tidak ada hit referensi ditemukan.")
            return seq_record
    except Exception as e:
        print(f"    -> [BLAST] Gagal menghubungi NCBI ({e}).")
        return seq_record

def main():
    print("=" * 65)
    print(" CHROMAFORGE: SANGER CONTIG ASSEMBLER (Strict CAP Mode)")
    print(" Author: Muhammad Aldiansyah Karisma | Universitas Bengkulu")
    print("=" * 65)
    
    folder = input("Ketik nama folder penyimpanan file .ab1 (Enter jika di folder ini): ").strip()
    pola = os.path.join(folder, "*.ab1") if folder else "*.ab1"
    
    pakai_blast = input("Aktifkan validasi orientasi otomatis BLASTn NCBI? [y/n]: ").strip().lower() == 'y'
    
    berkas_ab1 = glob.glob(pola)
    if not berkas_ab1:
        print("[!] Tidak ada file .ab1 ditemukan.")
        return
        
    pasangan_sampel = pasangkan_sampel_otomatis_ai(berkas_ab1)
    rekaman_fasta_iupac = []
    rekaman_fasta_pseudo = []
    semua_laporan_mismatch = []
    
    sampel_iupac_terdeteksi = []
    total_sampel = len(pasangan_sampel)
    
    print(f"\nTotal kelompok sampel terdeteksi: {total_sampel}\n")

    for urutan, id_s in enumerate(sorted(pasangan_sampel.keys(), key=urutkan_natural), 1):
        path_dict = pasangan_sampel[id_s]
        file_F = path_dict.get('F')
        file_R = path_dict.get('R')
        
        base_F = os.path.basename(file_F) if file_F else "Tidak ada Forward"
        base_R = os.path.basename(file_R) if file_R else "Tidak ada Reverse"

        if file_F or file_R:
            print(f"[{urutan}/{total_sampel}] Memproses sampel: {id_s} ({base_F} - {base_R})")
            
            f_rec = SeqIO.read(file_F, "abi") if file_F else None
            r_rec = SeqIO.read(file_R, "abi") if file_R else None
            
            seq_iupac, seq_pseudo, status, laporan_mismatch = proses_sampel_dengan_aturan_lengkap(f_rec, r_rec, id_s, min_length=200)
            
            if laporan_mismatch:
                semua_laporan_mismatch.append(laporan_mismatch)
                
            if seq_iupac is None or seq_pseudo is None: continue
            
            if str(seq_iupac) != str(seq_pseudo):
                sampel_iupac_terdeteksi.append(f"[{urutan}/{total_sampel}] {id_s}")
            
            record_iupac = SeqRecord(seq_iupac, id=id_s, description=f"Method: {status} | Mode: IUPAC")
            record_pseudo = SeqRecord(seq_pseudo, id=id_s, description=f"Method: {status} | Mode: PseudoPhased")
            
            if pakai_blast:
                record_iupac = koreksi_orientasi_via_blastn(record_iupac)
                record_pseudo = SeqRecord(record_iupac.seq, id=id_s, description=record_pseudo.description + " | BLAST_Orientation: Validated")
                time.sleep(20)
                
            rekaman_fasta_iupac.append(record_iupac)
            rekaman_fasta_pseudo.append(record_pseudo)

    if rekaman_fasta_iupac:
        if sampel_iupac_terdeteksi:
            nama_file_dnasp = "MultiFASTA_ChromaForge_For_DnaSP.fasta"
            nama_file_mega = "MultiFASTA_ChromaForge_For_MEGA.fasta"
            
            SeqIO.write(rekaman_fasta_iupac, nama_file_dnasp, "fasta")
            SeqIO.write(rekaman_fasta_pseudo, nama_file_mega, "fasta")
            
            print(f"\n[SELESAI] Heterozigot ganda (IUPAC) terdeteksi pada {len(sampel_iupac_terdeteksi)} sampel:")
            for s_iupac in sampel_iupac_terdeteksi:
                print(f"  -> {s_iupac}")
            print(f"\nFile diekspor menjadi:")
            print(f"  1. {nama_file_dnasp}")
            print(f"  2. {nama_file_mega}")
        else:
            nama_file_tunggal = "MultiFASTA_ChromaForge_Final.fasta"
            
            SeqIO.write(rekaman_fasta_iupac, nama_file_tunggal, "fasta")
            
            print(f"\n[SELESAI] Tidak ada kode IUPAC. File diekspor tunggal:")
            print(f"  1. {nama_file_tunggal} (Aman digunakan untuk MEGA maupun DnaSP)")
            
    if semua_laporan_mismatch:
        out_report = "Mismatch_Alignment_Report.txt"
        with open(out_report, "w", encoding="utf-8") as f:
            f.write("".join(semua_laporan_mismatch))

if __name__ == '__main__':
    main()