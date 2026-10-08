# ChromaForge

> [!WARNING]
> **BETA VERSION - TEST PURPOSES ONLY**
> Skrip ChromaForge saat ini murni berada dalam tahap pengujian beta (Beta Testing). Harap TIDAK menggunakan hasil perakitan dari skrip ini untuk pengolahan data primer pada skripsi, tesis, disertasi, atau publikasi jurnal hingga versi stabil resmi dirilis. Gunakan hanya untuk keperluan uji coba algoritma.

## Deskripsi
ChromaForge adalah instrumen bioinformatika berbasis Python yang mengekstraksi data biner kromatogram (.ab1) secara langsung. Sistem ini memangkas data kotor berdasarkan ambang batas Phred Q≥20 dan mengaplikasikan algoritma Fuzzy Logic untuk memasangkan sekuen Forward dan Reverse tanpa mewajibkan intervensi penamaan ulang secara manual.

## Fitur Utama
* **Native .ab1 Parsing**: Membaca langsung kromatogram tanpa perlu konversi format.
* **Phred Q≥20 Auto-Trim**: Pemangkasan ujung sekuen secara otomatis berdasarkan probabilitas kesalahan matematis.
* **AI Smart Match (Fuzzy Logic)**: Memasangkan file F dan R dari direktori sesuai pasangannya.
* **Fault Tolerance**: Mengisolasi file bermasalah atau rusak tanpa memicu sistem crash atau stop running.

## Langkah 1: Persiapan Ruang Kerja
1. Buat folder baru (misalnya "Proyek_DNA") dan pindahkan seluruh file kromatogram (.ab1) beserta skrip `ChromaForge_Assembler.py` ke dalamnya.
2. Pastikan file Forward dan Reverse Anda memiliki akhiran `_F` dan `_R` agar skrip otomatis memasangkannya.

## Langkah 2: Instalasi Kebutuhan Sistem
1. Unduh Python dari situs resminya. **SANGAT PENTING:** Saat menginstal, pastikan Anda mencentang kotak "Add Python to PATH" di bagian bawah layar, Jika tidak ada lanjutkan saja.
2. Buka Command Prompt (CMD), ketik `pip install biopython`, lalu tekan Enter untuk menginstal pustaka pengolah data.

## Langkah 3: Eksekusi Skrip
1. Buka folder "Proyek_DNA" Anda, klik kiri satu kali pada *Address Bar*, hapus teksnya, ketik `cmd`, lalu tekan Enter.
2. Pada jendela CMD yang terbuka, ketik `python ChromaForge_Assembler.py` dan tekan Enter.
3. Anda akan diminta memilih mode:
   * **Ketik n (Mode Offline)**: Sangat direkomendasikan. Memproses ratusan sampel dalam hitungan detik tanpa internet.
   * **Ketik y (Mode Online)**: Lambat karena mengirim data ke NCBI untuk mengecek arah untai. Cocok untuk 1-3 sampel saja.

## Hasil Akhir
Setelah proses selesai, skrip akan menghasilkan file `MultiFASTA_ChromaForge_BLAST.fasta` berisi seluruh DNA bersih yang telah dipotong dan digabungkan, siap dimasukkan ke aplikasi MEGA.

## Umpan Balik Pengujian (Beta Feedback)
Kami sangat menghargai evaluasi Anda terhadap performa algoritma ini. Jika Anda menemukan anomali, mohon laporkan melalui tautan berikut:
**[Formulir Evaluasi Beta ChromaForge](https://bit.ly/FormEvaluasiandBugTrackingChromaForge)**

## Identitas Pengembang
**Muhammad Aldiansyah Karisma**
Universitas Bengkulu (2026)

## Lisensi
Proyek ini didistribusikan di bawah **MIT License**. Lihat file LICENSE untuk detail perlindungan hukum dan hak cipta.
