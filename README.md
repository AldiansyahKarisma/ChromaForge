# ChromaForge
> [!WARNING]
> **BETA VERSION - TEST PURPOSES ONLY**
> Skrip ChromaForge saat ini murni berada dalam tahap pengujian beta (*Beta Testing*). Algoritma resolusi konflik Q-score dan integritas penanganan SNP masih dalam proses validasi silang dengan berbagai kondisi kromatogram. 
> **Harap TIDAK menggunakan hasil perakitan dari skrip ini untuk pengolahan data primer pada skripsi, tesis, disertasi, atau publikasi jurnal** hingga versi stabil resmi dirilis. Gunakan hanya untuk keperluan uji coba algoritma.

## Deskripsi
ChromaForge adalah instrumen bioinformatika berbasis Python untuk perakitan *contig* otomatis dari file *Sanger sequencing* (.ab1). Algoritma ini memitigasi bias pembacaan basa dengan menggunakan komparasi skor Phred (Q-score) mutlak untuk resolusi konflik dan *trimming* dinamis.

## Fitur Utama
*   **Robust Trimming & Masking**: Mengevaluasi dan memotong area derau di ujung sekuen (ambang batas rata-rata Q ≥ 20 dalam jendela baca).
*   **Strict Consensus Assembly**: Menyelesaikan konflik nukleotida pada area tumpang tindih (*overlap*) berdasarkan selisih margin kualitas > 10 poin. Area berimbang yang meragukan otomatis disamarkan menjadi 'N'.
*   **Automated BLASTn Orientation**: Validasi dan koreksi orientasi untai minus secara otomatis melalui API NCBI.
*   **Multi-FASTA Output**: Kompilasi akhir disimpan dalam format FASTA yang siap diolah lebih lanjut pada MEGA, BioEdit, atau DnaSP.

## Penggunaan (Beta)
Simpan skrip `ChromaForge_Assembler.py` satu folder dengan file sampel `.ab1` Anda, lalu jalankan via CMD atau terminal:
```bash
python ChromaForge_Assembler.py

Umpan Balik Pengujian (Beta Feedback)
Kami sangat menghargai evaluasi Anda terhadap performa algoritma ini. Jika Anda menemukan anomali resolusi konflik, kegagalan trimming, atau bug sistem saat memproses sampel kromatogram Anda, mohon laporkan melalui tautan berikut:
Formulir Evaluasi Beta ChromaForge

Identitas Pengembang
Muhammad Aldiansyah Karisma

Universitas Bengkulu (2026)

Lisensi
Proyek ini didistribusikan di bawah MIT License. Lihat file LICENSE untuk detail perlindungan hukum dan hak cipta.
