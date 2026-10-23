# ==============================================================================
# Script Name: ChromaForge - Automated Sanger Sequencing Assembler (Beta v2.0)
# Author: Muhammad Aldiansyah Karisma
# Affiliation: Biologi, Universitas Bengkulu
# Year: 2026
# License: MIT License
# Description: Intelligent contig assembler
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

# =====================================================================
# INISIALISASI GEMINI API
# =====================================================================
try:
    from google import genai
    from dotenv import load_dotenv
    LIBRARY_AI_TERSEDIA = True
except ImportError:
    LIBRARY_AI_TERSEDIA = False

ai_aktif = False
client = None

if LIBRARY_AI_TERSEDIA:
    load_dotenv()
    API_KEY = os.getenv("GEMINI_API_KEY")
    if API_KEY:
        client = genai.Client(api_key=API_KEY)
        ai_aktif = True

# =====================================================================
# LOGIKA UTAMA CHROMAFORGE
# =====================================================================
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
        
    good_count = sum(1 for q in trimmed_quals if q >= threshold_q)
    ratio_good = good_count / trim_length if trim_length > 0 else 0
    
    return trimmed_seq, ratio_good, trimmed_quals

def bangun_strict_consensus(f_seq, r_seq, alignments, f_qual, r_qual):
    alignment = alignments[0]
    t_span = alignment.aligned[0]
    q_span = alignment.aligned[1]
    
    if len(t_span) == 0 or len(q_span) == 0:
        return None
        
    f_idx = t_span[0][0]
    r_idx = q_span[0][0]
    
    consensus_iupac = list(f_seq[:f_idx])
    consensus_pseudo = list(f_seq[:f_idx])
    
    snp_count = 0
    log_mismatch = []
    
    for (t_start, t_end), (q_start, q_end) in zip(t_span, q_span):
        while f_idx < t_start:
            q_f = f_qual[f_idx] if (f_qual and f_idx < len(f_qual)) else 0
            if q_f >= 25: 
                consensus_iupac.append(f_seq[f_idx])
                consensus_pseudo.append(f_seq[f_idx])
            f_idx += 1
            
        while r_idx < q_start:
            q_r = r_qual[r_idx] if (r_qual and r_idx < len(r_qual)) else 0
            if q_r >= 25:
                consensus_iupac.append(r_seq[r_idx])
                consensus_pseudo.append(r_seq[r_idx])
            r_idx += 1
            
        for i in range(t_end - t_start):
            f_base = f_seq[t_start + i]
            r_base = r_seq[q_start + i]
            
            q_f = f_qual[t_start + i] if (f_qual and (t_start + i) < len(f_qual)) else 0
            q_r = r_qual[q_start + i] if (r_qual and (q_start + i) < len(r_qual)) else 0
            
            posisi_contig = len(consensus_iupac) + 1
            
            if f_base == r_base:
                consensus_iupac.append(f_base)
                consensus_pseudo.append(f_base)
            else:
                if q_f > q_r + 10:
                    consensus_iupac.append(f_base)
                    consensus_pseudo.append(f_base)
                    log_mismatch.append(f"Contig[{posisi_contig}]:F[{t_start + i}]/Rc[{q_start + i}]")
                elif q_r > q_f + 10:
                    consensus_iupac.append(r_base)
                    consensus_pseudo.append(r_base)
                    log_mismatch.append(f"Contig[{posisi_contig}]:F[{t_start + i}]/Rc[{q_start + i}]")
                else:
                    if q_f >= 30 and q_r >= 30:
                        kode_iupac = dapatkan_kode_iupac(f_base, r_base)
                        consensus_iupac.append(kode_iupac)
                        basa_terpilih = f_base if q_f >= q_r else r_base
                        consensus_pseudo.append(basa_terpilih)
                        snp_count += 1
                        log_mismatch.append(f"Contig[{posisi_contig}]:F[{t_start + i}]/Rc[{q_start + i}]({kode_iupac})")
                    else:
                        if max(q_f, q_r) >= 20:
                            basa_terpilih = f_base if q_f >= q_r else r_base
                            consensus_iupac.append(basa_terpilih)
                            consensus_pseudo.append(basa_terpilih)
                        else:
                            consensus_iupac.append('N')
                            consensus_pseudo.append('N')
                            snp_count += 1
                            log_mismatch.append(f"Contig[{posisi_contig}]:F[{t_start + i}]/Rc[{q_start + i}](N)")
                            
        f_idx = t_end
        r_idx = q_end
                    
    r_end_overlap = q_span[-1][1]
    consensus_iupac.extend(list(r_seq[r_end_overlap:]))
    consensus_pseudo.extend(list(r_seq[r_end_overlap:]))
    
    return "".join(consensus_iupac), "".join(consensus_pseudo), snp_count, log_mismatch

def proses_sampel_dengan_aturan_lengkap(f_rec, r_rec, id_sampel, min_length=200):
    ambang_layak = 0.85 
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
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -2
        
        alignments = aligner.align(f_proc, r_rc)
        if not alignments:
            if f_score >= r_score:
                if len(f_proc) < min_length: 
                    print(f"    -> {WARNA_MERAH}[SKIPPED] Contig gagal & F terlalu pendek.{RESET_WARNA}")
                    return None, None, "SKIPPED", None
                print(f"    -> {WARNA_ORANYE}[Fallback] Contig gagal tersambung, menggunakan F (Q: {f_score*100:.1f}%, Len: {len(f_proc)} bp).{RESET_WARNA}")
                return Seq(f_proc), Seq(f_proc), f"Contig gagal, fallback F (Q: {f_score*100:.1f}%)", None
            else:
                if len(r_proc) < min_length: 
                    print(f"    -> {WARNA_MERAH}[SKIPPED] Contig gagal & R terlalu pendek.{RESET_WARNA}")
                    return None, None, "SKIPPED", None
                print(f"    -> {WARNA_ORANYE}[Fallback] Contig gagal tersambung, menggunakan R-rc (Q: {r_score*100:.1f}%, Len: {len(r_proc)} bp).{RESET_WARNA}")
                return Seq(r_rc), Seq(r_rc), f"Contig gagal, fallback R-rc (Q: {r_score*100:.1f}%)", None
                
        hasil_consensus = bangun_strict_consensus(f_proc, r_rc, alignments, f_qual, r_qual_rc)
        
        if hasil_consensus:
            seq_final_iupac, seq_final_pseudo, snp_total, log_mismatch = hasil_consensus
            
            if len(seq_final_iupac) < min_length:
                print(f"    -> {WARNA_MERAH}[SKIPPED] Contig terlalu pendek ({len(seq_final_iupac)} bp).{RESET_WARNA}")
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
            
        print(f"    -> {WARNA_ORANYE}[Fallback] Consensus error, menggunakan F (Q: {f_score*100:.1f}%, Len: {len(f_proc)} bp).{RESET_WARNA}")
        return Seq(f_proc), Seq(f_proc), f"Contig fallback F", None
    elif f_layak:
        print(f"    -> {WARNA_KUNING}[Info] Kualitas R buruk ({r_score*100:.1f}%). Hanya menggunakan F (Q: {f_score*100:.1f}%, Len: {len(f_proc)} bp).{RESET_WARNA}")
        return Seq(f_proc), Seq(f_proc), f"Hanya F (F layak: {f_score*100:.1f}%)", None
    elif r_layak:
        seq_r_rc = Seq(r_proc).reverse_complement()
        print(f"    -> {WARNA_KUNING}[Info] Kualitas F buruk ({f_score*100:.1f}%). Hanya menggunakan R-rc (Q: {r_score*100:.1f}%, Len: {len(r_proc)} bp).{RESET_WARNA}")
        return seq_r_rc, seq_r_rc, f"Hanya R-rc (R layak: {r_score*100:.1f}%)", None
    else:
        if f_score >= r_score:
            seq_recovery = f_proc
            seq_final = Seq(seq_recovery)
            tipe_bacaan = "F"
        else:
            seq_recovery = r_proc
            seq_final = Seq(seq_recovery).reverse_complement()
            tipe_bacaan = "R-rc"
            
        status_metode = f"Best Single-Read ({tipe_bacaan})"
        print(f"    -> {WARNA_ORANYE}[Warning] Kualitas sangat buruk. Menggunakan {status_metode} (Len: {len(seq_recovery)} bp).{RESET_WARNA}")
        return seq_final, seq_final, status_metode, None

def ekstrak_id_cerdas(clean_name):
    core = re.sub(r'^1st_BASE_\d+_', '', clean_name, flags=re.IGNORECASE)
    core = re.sub(r'_(resend|repeat|F_A|R_A)$', '', core, flags=re.IGNORECASE)
    core = re.sub(r'_(27F|1492R|F|R)$', '', core, flags=re.IGNORECASE)
    core = re.sub(r'_[A-Ha-h]\d{2}_\d+.*$', '', core)
    core = re.sub(r'_(fastseq|pop|macrogen|apical|z)[a-zA-Z0-9_]*$', '', core, flags=re.IGNORECASE)
    core = re.sub(r'(?<!^)(F|R)$', '', core, flags=re.IGNORECASE)
    return core.strip('_')

# LOGIKA PENCARIAN & PEMASANGAN SAMPEL YANG UNIVERSAL & KETAT (MENCEGAH SAMPEL TERTIMPA/LOMPAT)
def pasangkan_sampel_otomatis_ai(berkas_ab1):
    f_files = {}
    r_files = {}
    
    for path in berkas_ab1:
        nm_file = os.path.basename(path)
        nm_for_detect = re.sub(r'_(resend|repeat)', '', nm_file, flags=re.IGNORECASE).upper()
        arah = 'R' if any(k in nm_for_detect for k in ['_R.', '_R_', 'BR1', 'GADPHR', 'CYTBR', '16SD_R', 'R_G', 'R1.', '1492R']) or nm_for_detect.replace('.AB1','').endswith('R') else 'F'
        clean_name = re.sub(r'\.ab1$', '', nm_file, flags=re.IGNORECASE)
        
        id_sampel = ekstrak_id_cerdas(clean_name)
        if arah == 'F':
            f_files[id_sampel] = path
        else:
            r_files[id_sampel] = path
            
    semua_id = sorted(list(set(f_files.keys()).union(set(r_files.keys()))), key=urutkan_natural)
    pasangan_sampel = {}
    
    for id_s in semua_id:
        pasangan_sampel[id_s] = {
            'F': f_files.get(id_s),
            'R': r_files.get(id_s)
        }
            
    return pasangan_sampel

def urutkan_natural(teks):
    return [int(c) if c.isdigit() else c.lower() for c in re.split(r'(\d+)', teks)]

def dry_run_koreksi(pasangan_sampel):
    while True:
        pasangan_sampel = {k: v for k, v in pasangan_sampel.items() if v['F'] or v['R']}
        list_id_terurut = sorted(pasangan_sampel.keys(), key=urutkan_natural)
        
        print(f"\n{WARNA_BIRU}[DRY-RUN] Evaluasi Draf Pasangan Sampel:{RESET_WARNA}")
        for i, id_s in enumerate(list_id_terurut, 1):
            f_name = os.path.basename(pasangan_sampel[id_s]['F']) if pasangan_sampel[id_s]['F'] else "None"
            r_name = os.path.basename(pasangan_sampel[id_s]['R']) if pasangan_sampel[id_s]['R'] else "None"
            print(f"  {i}. {id_s:<15} (F: {f_name:<25} | R: {r_name:<25})")
            
        print(f"\n{WARNA_KUNING}Petunjuk Koreksi Manual:{RESET_WARNA}")
        print(" - 'tukar [F/R] [no1] [no2]'  : Menukar posisi antar dua sampel")
        print(" - 'geser [F/R] [asal] [tuj]' : Memindah file ke slot tujuan")
        print(" - 'pisah [F/R] [no]'         : Mencabut file menjadi sampel mandiri")
        print(" - 'y'                        : Lanjut proses perakitan")
        
        pilihan = input("\nInput Anda: ").strip().lower()
        if pilihan == 'y':
            break
            
        parts = pilihan.split()
        if len(parts) >= 3:
            cmd = parts[0]
            arah = parts[1].upper()
            if arah not in ['F', 'R']:
                print(f"{WARNA_MERAH}[!] Arah harus F atau R.{RESET_WARNA}")
                continue
                
            if cmd == 'tukar' and len(parts) == 4:
                try:
                    idx1, idx2 = int(parts[2]) - 1, int(parts[3]) - 1
                    id1, id2 = list_id_terurut[idx1], list_id_terurut[idx2]
                    pasangan_sampel[id1][arah], pasangan_sampel[id2][arah] = pasangan_sampel[id2][arah], pasangan_sampel[id1][arah]
                    print(f"-> {WARNA_UNGU}[UPDATE] File {arah} pada sampel {idx1+1} dan {idx2+1} berhasil ditukar.{RESET_WARNA}")
                except (ValueError, IndexError):
                    print(f"{WARNA_MERAH}[!] Nomor sampel tidak valid.{RESET_WARNA}")
            elif cmd == 'geser' and len(parts) == 4:
                try:
                    idx_asal, idx_tuj = int(parts[2]) - 1, int(parts[3]) - 1
                    id_asal, id_tuj = list_id_terurut[idx_asal], list_id_terurut[idx_tuj]
                    file_dipindah = pasangan_sampel[id_asal][arah]
                    if not file_dipindah:
                        print(f"{WARNA_MERAH}[!] Tidak ada file {arah} di sampel asal.{RESET_WARNA}")
                        continue
                    pasangan_sampel[id_tuj][arah] = file_dipindah
                    pasangan_sampel[id_asal][arah] = None
                    nm_file = os.path.basename(file_dipindah)
                    print(f"-> {WARNA_UNGU}[UPDATE] File {arah} ({nm_file}) dipindahkan dari sampel {idx_asal+1} ke {idx_tuj+1}.{RESET_WARNA}")
                except (ValueError, IndexError):
                    print(f"{WARNA_MERAH}[!] Nomor sampel tidak valid.{RESET_WARNA}")
            elif cmd == 'pisah' and len(parts) == 3:
                try:
                    idx = int(parts[2]) - 1
                    id_asal = list_id_terurut[idx]
                    file_dicabut = pasangan_sampel[id_asal][arah]
                    if not file_dicabut:
                        print(f"{WARNA_MERAH}[!] Tidak ada file {arah} di sampel tersebut.{RESET_WARNA}")
                        continue
                    nm_file = os.path.basename(file_dicabut)
                    new_id = ekstrak_id_cerdas(re.sub(r'\.ab1$', '', nm_file, flags=re.IGNORECASE)) + f"_pisah_{arah.lower()}"
                    counter = 1
                    base_new_id = new_id
                    while new_id in pasangan_sampel:
                        new_id = f"{base_new_id}_{counter}"
                        counter += 1
                    pasangan_sampel[new_id] = {'F': None, 'R': None}
                    pasangan_sampel[new_id][arah] = file_dicabut
                    pasangan_sampel[id_asal][arah] = None
                    print(f"-> {WARNA_UNGU}[UPDATE] File {arah} ({nm_file}) dicabut dari sampel {idx+1} dan dibuatkan entitas baru.{RESET_WARNA}")
                except (ValueError, IndexError):
                    print(f"{WARNA_MERAH}[!] Nomor sampel tidak valid.{RESET_WARNA}")
            else:
                print(f"{WARNA_MERAH}[!] Format perintah salah.{RESET_WARNA}")
        else:
            print(f"{WARNA_MERAH}[!] Perintah tidak dikenali.{RESET_WARNA}")
            
    return pasangan_sampel

def koreksi_orientasi_via_blastn(seq_record):
    print(f"    -> [BLAST] Mengirim {seq_record.id} ke NCBI BLASTn... (Jeda aman 20 detik aktif)")
    try:
        result_handle = NCBIWWW.qblast("blastn", "nt", seq_record.seq, hitlist_size=1)
        blast_record = NCBIXML.read(result_handle)
        if blast_record.alignments:
            hsp = blast_record.alignments[0].hsps[0]
            if hsp.frame[1] < 0:
                print(f"    -> [BLAST] Orientasi Minus terdeteksi. Memutar orientasi...")
                return SeqRecord(seq_record.seq.reverse_complement(), id=seq_record.id, name="", description="")
            else:
                print(f"    -> [BLAST] Orientasi Plus/Plus sudah akurat.")
                return SeqRecord(seq_record.seq, id=seq_record.id, name="", description="")
        print(f"    -> [BLAST] Tidak ada hit referensi ditemukan.")
        return seq_record
    except Exception as e:
        print(f"    -> [BLAST] Gagal menghubungi NCBI ({e}).")
        return seq_record

def analyze_report_with_ai(report_path):
    if not ai_aktif or not client: return
    print(f"\n{WARNA_BIRU}[ChromaCore AI] Menganalisis laporan mismatch...{RESET_WARNA}")
    try:
        with open(report_path, "r", encoding="utf-8") as file:
            report_content = file.read()
        prompt = f"""Anda adalah asisten bioinformatika ChromaCore AI. Evaluasi laporan ini: {report_content}.
        Format dua poin teks murni (tanpa markdown):
        1. Ringkasan Kualitas Perakitan (Tinggi, Menengah, Rendah/Gagal).
        2. Instruksi Inspeksi In Silico manual.
        Gunakan gaya bahasa absolut."""
        
        response = client.models.generate_content(model='gemini-3.8-flash', contents=prompt)
        teks_reminder = f"\n\n{WARNA_KUNING}REMINDER:\nHalo! Sekadar mengingatkan untuk pengerjaan lab berikutnya ya. Pastikan semua parameter PCR dipantau dengan baik, dan usahakan proses ekstraksi DNA dilakukan sebersih dan sesteril mungkin. Langkah ini sangat krusial supaya kita terhindar dari kontaminasi silang dan bisa mendapatkan kualitas sekuen yang mulus tanpa banyak derau. Tetap semangat risetnya!{RESET_WARNA}"
        
        print(f"\n{WARNA_UNGU}================ EVALUASI CHROMACORE AI ================{RESET_WARNA}")
        print(response.text + teks_reminder)
        print(f"{WARNA_UNGU}========================================================{RESET_WARNA}\n")
        
        with open("AI_Analysis_Result.txt", "w", encoding="utf-8") as f:
            f.write(response.text + "\n\nREMINDER:\nHalo! Sekadar mengingatkan untuk pengerjaan lab berikutnya ya. Pastikan semua parameter PCR dipantau dengan baik, dan usahakan proses ekstraksi DNA dilakukan sebersih dan sesteril mungkin. Langkah ini sangat krusial supaya kita terhindar dari kontaminasi silang dan bisa mendapatkan kualitas sekuen yang mulus tanpa banyak derau. Tetap semangat risetnya!")
    except Exception as e:
        print(f"{WARNA_MERAH}[ERROR] Gagal melakukan analisis AI: {e}{RESET_WARNA}")

# =====================================================================
# EKSEKUSI UTAMA
# =====================================================================
def main():
    print("=" * 65)
    print(" CHROMAFORGE: SANGER CONTIG ASSEMBLER (Beta v2.0)")
    print(" Author: Muhammad Aldiansyah Karisma | Biologi, Universitas Bengkulu")
    if ai_aktif: print(f" Status ChromaCore AI: {WARNA_BIRU}AKTIF{RESET_WARNA}")
    elif LIBRARY_AI_TERSEDIA: print(f" Status ChromaCore AI: {WARNA_MERAH}TIDAK AKTIF (API Key tidak valid){RESET_WARNA}")
    else: print(f" Status ChromaCore AI: {WARNA_MERAH}TIDAK AKTIF (Pustaka belum terinstal){RESET_WARNA}")
    print("=" * 65)
    
    folder = input("Ketik nama folder penyimpanan file .ab1 (Enter jika di folder ini): ").strip()
    pola = os.path.join(folder, "*.ab1") if folder else "*.ab1"
    
    pakai_blast = input("Aktifkan validasi orientasi otomatis BLASTn NCBI? [y/n]: ").strip().lower() == 'y'
    
    berkas_ab1 = glob.glob(pola)
    if not berkas_ab1:
        print("[!] Tidak ada file .ab1 ditemukan.")
        return
        
    pasangan_sampel = pasangkan_sampel_otomatis_ai(berkas_ab1)
    
    # Eksekusi Dry-Run Interaktif
    pasangan_sampel = dry_run_koreksi(pasangan_sampel)
    
    list_id_terurut = sorted(pasangan_sampel.keys(), key=urutkan_natural)
    
    print(f"\n{WARNA_BIRU}Petunjuk Akhir: Ketik nomor sampel (contoh: 1, 2, 4, 7 atau rentang 1-8).{RESET_WARNA}")
    pilihan = input("Ketik 'all' jika ingin memproses semua sampel di atas: ").strip().lower()
    
    if pilihan != 'all' and pilihan != '':
        indeks_terpilih = set()
        for bagian in pilihan.split(','):
            bagian = bagian.strip()
            if '-' in bagian:
                try:
                    awal, akhir = map(int, bagian.split('-'))
                    indeks_terpilih.update(range(awal, akhir + 1))
                except ValueError: continue
            elif bagian.isdigit():
                indeks_terpilih.add(int(bagian))
        pasangan_sampel = {list_id_terurut[idx-1]: pasangan_sampel[list_id_terurut[idx-1]] for idx in indeks_terpilih if 1 <= idx <= len(list_id_terurut)}

    nama_custom = input("\nKetik nama file keluaran / awalan (tanpa ekstensi .fasta): ").strip()
    if not nama_custom: nama_custom = "Output_ChromaForge"

    total_sampel = len(pasangan_sampel)
    print(f"\nMemproses {total_sampel} sampel...\n")
    
    semua_laporan_mismatch = []
    koleksi_fasta_iupac = []
    koleksi_fasta_pseudo = []
    koleksi_fasta_tunggal = []
    
    for urutan, id_s in enumerate(sorted(pasangan_sampel.keys(), key=urutkan_natural), 1):
        print(f"[{urutan}/{total_sampel}] Memproses sampel: {id_s}")
        f_rec = SeqIO.read(pasangan_sampel[id_s]['F'], "abi") if pasangan_sampel[id_s]['F'] else None
        r_rec = SeqIO.read(pasangan_sampel[id_s]['R'], "abi") if pasangan_sampel[id_s]['R'] else None
        
        seq_iupac, seq_pseudo, status, laporan_mismatch = proses_sampel_dengan_aturan_lengkap(f_rec, r_rec, id_s, min_length=200)
        if laporan_mismatch: semua_laporan_mismatch.append(laporan_mismatch)
        if seq_iupac is None: continue
        
        record_iupac = SeqRecord(seq_iupac, id=id_s, name="", description="")
        record_pseudo = SeqRecord(seq_pseudo, id=id_s, name="", description="")
        
        if pakai_blast:
            record_iupac = koreksi_orientasi_via_blastn(record_iupac)
            record_pseudo = SeqRecord(record_iupac.seq, id=id_s, name="", description="")
            time.sleep(20)
            
        # Ekspor Satuan (Individu per Sampel)
        if str(seq_iupac) != str(seq_pseudo):
            SeqIO.write([record_iupac], f"{nama_custom}_{id_s}_For_DnaSP.fasta", "fasta")
            SeqIO.write([record_pseudo], f"{nama_custom}_{id_s}_For_MEGA.fasta", "fasta")
            print(f"  -> [EKSPOR SATUAN] {WARNA_UNGU}(Heterozigot){RESET_WARNA} {nama_custom}_{id_s}_For_DnaSP.fasta & {nama_custom}_{id_s}_For_MEGA.fasta (bisa digunakan untuk MEGA maupun DnaSP)")
            
            koleksi_fasta_iupac.append(record_iupac)
            koleksi_fasta_pseudo.append(record_pseudo)
        else:
            SeqIO.write([record_iupac], f"{nama_custom}_{id_s}_Final.fasta", "fasta")
            print(f"  -> [EKSPOR SATUAN] {nama_custom}_{id_s}_Final.fasta (bisa digunakan untuk MEGA maupun DnaSP)")
            
            koleksi_fasta_tunggal.append(record_iupac)

    # Ekspor Gabungan (MultiFASTA Sekaligus)
    print(f"\n{WARNA_BIRU}[EKSPOR GABUNGAN] Membuat file MultiFASTA...{RESET_WARNA}")
    if koleksi_fasta_iupac:
        SeqIO.write(koleksi_fasta_iupac, f"{nama_custom}_MultiFASTA_For_DnaSP.fasta", "fasta")
        SeqIO.write(koleksi_fasta_pseudo, f"{nama_custom}_MultiFASTA_For_MEGA.fasta", "fasta")
        print(f"  -> Gabungan IUPAC: {nama_custom}_MultiFASTA_For_DnaSP.fasta & {nama_custom}_MultiFASTA_For_MEGA.fasta (bisa digunakan untuk MEGA maupun DnaSP)")
        
    semua_rekaman = koleksi_fasta_tunggal + koleksi_fasta_iupac
    if semua_rekaman:
        SeqIO.write(semua_rekaman, f"{nama_custom}_MultiFASTA_Final.fasta", "fasta")
        print(f"  -> Gabungan Tunggal: {nama_custom}_MultiFASTA_Final.fasta (bisa digunakan untuk MEGA maupun DnaSP)")
            
    if semua_laporan_mismatch:
        out_report = "Mismatch_Alignment_Report.txt"
        with open(out_report, "w", encoding="utf-8") as f:
            f.write("".join(semua_laporan_mismatch))
        analyze_report_with_ai(out_report)

if __name__ == '__main__':
    main()