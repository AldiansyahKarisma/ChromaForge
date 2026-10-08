# ==============================================================================
# Script Name: ChromaForge - Automated Sanger Sequencing Assembler (Beta)
# Author: Muhammad Aldiansyah Karisma
# Affiliation: Universitas Bengkulu
# Year: 2026
# License: MIT License
# Description: Automated script for resolving contig conflicts using Phred Q-scores
#              from Sanger sequencing chromatograms (.ab1).
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
    
    # PERBAIKAN: Menggunakan len() untuk menghindari ValueError pada array numpy Biopython versi baru
    if len(t_span) == 0 or len(q_span) == 0:
        return None
        
    f_start_overlap = t_span[0][0]
    consensus = list(f_seq[:f_start_overlap])
    
    snp_count = 0
    for t_block, q_block in zip(t_span, q_span):
        for i in range(t_block[1] - t_block[0]):
            f_idx = t_block[0] + i
            q_idx = q_block[0] + i
            
            f_base = f_seq[f_idx]
            r_base = r_seq[q_idx]
            
            if f_base == r_base:
                consensus.append(f_base)
            else:
                q_f = f_qual[f_idx] if f_qual else 0
                q_r = r_qual[q_idx] if r_qual else 0
                
                if q_f > q_r + 10:
                    consensus.append(f_base)
                elif q_r > q_f + 10:
                    consensus.append(r_base)
                else:
                    consensus.append('N')
                    snp_count += 1
                    
    r_end_overlap = q_span[-1][1]
    consensus.extend(list(r_seq[r_end_overlap:]))
    
    return "".join(consensus), snp_count

def proses_sampel_dengan_aturan_lengkap(f_rec, r_rec):
    f_proc, f_score, f_qual = robust_trim_masking_dan_evaluasi(f_rec)
    r_proc, r_score, r_qual_raw = robust_trim_masking_dan_evaluasi(r_rec)
    
    ambang_layak = 0.85 
    ambang_minimum = 0.70
    
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
                return Seq(f_proc), f"Contig gagal, fallback F (Q: {f_score*100:.1f}%)"
            else:
                return Seq(r_rc), f"Contig gagal, fallback R-rc (Q: {r_score*100:.1f}%)"
                
        hasil_consensus = bangun_strict_consensus(f_proc, r_rc, alignments, f_qual, r_qual_rc)
        
        if hasil_consensus:
            seq_final, snp_total = hasil_consensus
            status = f"Contig CAP Strict (SNP/Konflik: {snp_total} titik)"
            return Seq(seq_final), status
            
        return Seq(f_proc), f"Contig fallback F"
        
    elif f_layak and not r_layak:
        print(f"    -> [Info] R kurang layak. Pakai F saja (Q: {f_score*100:.1f}%).")
        return Seq(f_proc), f"Hanya F (F layak: {f_score*100:.1f}%)"
        
    elif not f_layak and r_layak:
        print(f"    -> [Info] F kurang layak. Pakai R-rc saja (Q: {r_score*100:.1f}%).")
        return Seq(r_proc).reverse_complement(), f"Hanya R-rc (R layak: {r_score*100:.1f}%)"
        
    else:
        tertinggi = max(f_score, r_score)
        if tertinggi < ambang_minimum:
            print(f"    -> [SKIPPED] F dan R di bawah 70%. Cek manual.")
            return None, "SKIPPED"
            
        if f_score >= r_score:
            print(f"    -> [Warning] F dan R < 85%. Pakai F (Q: {f_score*100:.1f}%).")
            return Seq(f_proc), f"F (Tertinggi: {f_score*100:.1f}%)"
        else:
            print(f"    -> [Warning] F dan R < 85%. Pakai R-rc (Q: {r_score*100:.1f}%).")
            return Seq(r_proc).reverse_complement(), f"R-rc (Tertinggi: {r_score*100:.1f}%)"

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

def pasangkan_sampel_otomatis_ai(berkas_ab1):
    f_files = []
    r_files = []
    
    for path in berkas_ab1:
        nm_file = os.path.basename(path)
        name = nm_file.replace('.ab1', '').replace('.AB1', '')
        
        name_no_resend = re.sub(r'_resend', '', name, flags=re.IGNORECASE)
        upper = name_no_resend.upper()
        
        if any(k in upper for k in ['_R', 'BR1', 'GADPHR', 'CYTBR', 'BFIB5R', 'R_G', 'R1']) or upper.endswith('R'):
            arah = 'R'
        elif any(k in upper for k in ['_F', 'BF1', 'GADPHF', 'CYTBF', 'BFIB5F', 'F_G', 'F1']) or upper.endswith('F'):
            arah = 'F'
        else:
            arah = 'R' if upper.endswith('R') else 'F'
            
        clean_name = re.sub(r'^1st_BASE_\d+_', '', name, flags=re.IGNORECASE)
        clean_name = re.sub(r'_[FR](?:_A|1)?$', '', clean_name, flags=re.IGNORECASE)
        clean_name = re.sub(r'_resend$', '', clean_name, flags=re.IGNORECASE)
        
        if arah == 'F':
            f_files.append((clean_name, path))
        else:
            r_files.append((clean_name, path))
            
    pasangan_sampel = {}
    
    for f_clean, f_path in f_files:
        best_match_r = None
        best_score = 0
        r_clean_best = ""
        
        for r_clean, r_path in r_files:
            score = SequenceMatcher(None, f_clean, r_clean).ratio()
            if score > best_score:
                best_score = score
                best_match_r = r_path
                r_clean_best = r_clean
                
        if best_match_r and best_score > 0.3:
            t1 = f_clean.split('_')
            t2 = r_clean_best.split('_')
            common = []
            for w1, w2 in zip(t1, t2):
                if w1 == w2:
                    common.append(w1)
                else:
                    break
            
            if not common:
                prefix = os.path.commonprefix([f_clean, r_clean_best])
                id_sampel = re.sub(r'[_-]+$', '', prefix)
            else:
                id_sampel = '_'.join(common)
                
            if not id_sampel:
                id_sampel = f_clean 
                
            pasangan_sampel[id_sampel] = {'F': f_path, 'R': best_match_r}
            r_files = [(c, p) for c, p in r_files if p != best_match_r]
            
    return pasangan_sampel

def main():
    print("=" * 65)
    print(" CHROMAFORGE: SANGER CONTIG ASSEMBLER (Strict CAP Mode)")
    print(" Author: Muhammad Aldiansyah Karisma | Universitas Bengkulu")
    print("=" * 65)
    
    folder = input("Ketik nama folder penyimpanan file .ab1 (Enter jika di folder ini): ").strip()
    pola = os.path.join(folder, "*.ab1") if folder else "*.ab1"
    
    pakai_blast = input("Aktifkan validasi orientasi otomatis BLASTn NCBI? (Aman tapi lambat) [y/n]: ").strip().lower() == 'y'
    
    berkas_ab1 = glob.glob(pola)
    if not berkas_ab1:
        print("[!] Tidak ada file .ab1 ditemukan.")
        return
        
    pasangan_sampel = pasangkan_sampel_otomatis_ai(berkas_ab1)
    rekaman_fasta = []
    
    print(f"\nTotal kelompok sampel terdeteksi: {len(pasangan_sampel)}\n")

    for id_s, path_dict in pasangan_sampel.items():
        file_F = path_dict.get('F')
        file_R = path_dict.get('R')
        
        base_F = os.path.basename(file_F) if file_F else "Tidak ada Forward"
        base_R = os.path.basename(file_R) if file_R else "Tidak ada Reverse"

        if file_F and file_R:
            print(f" [OK] Memproses sampel: {id_s} ({base_F} - {base_R})")
            f_rec = SeqIO.read(file_F, "abi")
            r_rec = SeqIO.read(file_R, "abi")
            
            sekuen_final, status = proses_sampel_dengan_aturan_lengkap(f_rec, r_rec)
            if sekuen_final is None: continue
            
            final_record = SeqRecord(sekuen_final, id=id_s, description=f"Method: {status}")
            
            if pakai_blast:
                final_record = koreksi_orientasi_via_blastn(final_record)
                time.sleep(20)
                
            rekaman_fasta.append(final_record)
        else:
            print(f" [!] Gagal dipasangkan / Tidak Lengkap: {id_s} ({base_F} - {base_R})")

    if rekaman_fasta:
        out_name = "MultiFASTA_ChromaForge_Strict.fasta"
        SeqIO.write(rekaman_fasta, out_name, "fasta")
        print(f"\n[SELESAI] File disimpan sebagai: {out_name}")

if __name__ == '__main__':
    main()