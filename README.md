# ChromaForge: Automated Sanger Sequencing Assembler

⚠️ **BETA VERSION - TEST PURPOSES ONLY**
Saat ini, ChromaForge masih berada dalam tahap pengujian beta (Beta Testing). Harap TIDAK menggunakan hasil perakitan skrip ini untuk pengolahan data primer pada skripsi, tesis, disertasi, atau publikasi jurnal hingga versi stabil (stable release) dirilis. Gunakan hanya untuk keperluan uji coba algoritma.

---

## Deskripsi
ChromaForge adalah instrumen bioinformatika berbasis Python yang mengekstraksi data biner kromatogram (`.ab1`) secara langsung. Sistem ini memangkas data kotor berdasarkan ambang batas Phred Q≥20 dan mengaplikasikan algoritma *Fuzzy Logic* untuk memasangkan sekuen *Forward* dan *Reverse* secara cerdas tanpa mewajibkan intervensi penamaan ulang manual.

---

## Fitur Utama
*   **Native .ab1 Parsing:** Membaca langsung file biner kromatogram tanpa perlu konversi format.
*   **Phred Q≥20 Auto-Trim:** Pemangkasan ujung sekuen secara otomatis berdasarkan probabilitas kesalahan matematis (Phred Q-score).
*   **AI Smart Match (Fuzzy Logic):** Memasangkan file F dan R dari direktori secara cerdas berdasarkan kemiripan nama sampel (>45% *string matching*).
*   **Integrasi ChromaCore AI:** Evaluasi pasca-perakitan otomatis menggunakan *Large Language Model* (LLM) untuk mendiagnosis kualitas sekuen.
*   **Penanganan Heterozigot (IUPAC):** Mempertahankan kode ambiguitas asli untuk DnaSP dan melakukan *pseudo-phasing* untuk kompatibilitas MEGA.
*   **Fault Tolerance:** Mengisolasi file bermasalah atau rusak tanpa memicu sistem berhenti beroperasi (*crash*).

---

## Panduan Penggunaan

### Langkah 1: Persiapan Ruang Kerja
1. Buat folder baru di komputer Anda (misal: `Proyek_DNA`).
2. Kumpulkan seluruh file kromatogram (`.ab1`) ke dalam folder tersebut.
3. Unduh dan letakkan file skrip `ChromaForge_Assembler.py` ke folder yang sama.

> **💡 Catatan AI Smart Match:** Pastikan file memiliki penanda F dan R (misal akhiran `_F` dan `_R`). Anda tidak perlu mengubah nama secara spesifik, sistem otomatis memasangkannya.

### Langkah 2: Konfigurasi ChromaCore AI (Opsional)
Untuk mengaktifkan asisten evaluasi otomatis (ChromaCore AI):
1. Buat API Key melalui [Google AI Studio](https://aistudio.google.com/).
2. Buka program **Notepad**.
3. Ketikkan kode berikut di baris pertama: `GEMINI_API_KEY=kunci_rahasia_anda_di_sini`
4. Klik **File > Save As**.
5. Pada bagian *Save as type*, ubah menjadi **All Files (*.*)**.
6. Pada *File name*, ketik persis **.env**
7. Simpan di dalam folder kerja Anda.

### Langkah 3: Instalasi Kebutuhan Sistem
1. Unduh dan instal **Python 3.x** dari [situs resmi Python](https://www.python.org/downloads/). **WAJIB** centang *Add Python to PATH* saat instalasi.
2. Buka Command Prompt (CMD) dan instal pustaka yang dibutuhkan dengan menjalankan perintah berikut:

```cmd
pip install biopython google-genai python-dotenv
```

### Langkah 4: Menjalankan Skrip
1. Buka folder kerja Anda, klik pada *Address Bar* (kolom lokasi folder di bagian atas), ketik `cmd`, lalu tekan Enter.
2. Di jendela CMD yang terbuka, jalankan perintah:

```cmd
python ChromaForge_Assembler.py
```

3. Tekan Enter saat diminta nama folder.
4. Pilih mode orientasi BLASTn: ketik `n` (Offline - Sangat Cepat) atau `y` (Online - Melakukan kueri ke server NCBI).

---

## Keluaran (Output)
Berdasarkan deteksi kode IUPAC dan status API, skrip menghasilkan hingga 4 file keluaran di folder kerja Anda:
1. `MultiFASTA_ChromaForge_For_MEGA.fasta`: Sekuen konsensus dengan *pseudo-phasing* basa ganda. Aman untuk MEGA.
2. `MultiFASTA_ChromaForge_For_DnaSP.fasta`: Sekuen yang mempertahankan kode ambiguitas IUPAC asli (jika terdapat heterozigot).
3. `Mismatch_Alignment_Report.txt`: Catatan teknis (log) posisi *mismatch* dan resolusi basa.
4. `AI_Analysis_Result.txt`: Evaluasi naratif kualitas sampel dan instruksi *in silico* dari ChromaCore AI.

---

## Umpan Balik (Feedback)
Kami sangat menghargai evaluasi kritis terhadap performa algoritma ini. Jika Anda menemukan anomali atau bug, mohon laporkan melalui tautan berikut:
[Formulir Evaluasi dan Bug Tracking ChromaForge](https://bit.ly/FormEvaluasiandBugTrackingChromaForge)

---

## Identitas Pengembang
**Muhammad Aldiansyah Karisma** | Universitas Bengkulu (2026)

## Lisensi
Proyek ini didistribusikan di bawah **MIT License**. Lihat file `LICENSE` untuk detail perlindungan hukum dan hak cipta.
