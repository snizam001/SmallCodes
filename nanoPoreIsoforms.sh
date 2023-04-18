genome="/media/txpn/nvme/Databases/hg38/GRCh38.p13.fa"

for f in OF10 OF11 OF12 OF13 OF14 OF15 OF16 OF17 OF18 OF19 OF1 OF20 OF21 OF22 OF23 OF24 OF25 OF26 OF27 OF28 OF29 OF2 OF30 OF31 OF32 OF33 OF34 OF3 OF4 OF5 OF6 OF7 OF8 OF9

do 

#alignment

minimap2 -ax splice -t 30 --secondary=no /media/txpn/nvme/Databases/hg38/GRCh38.p13.fa Sample_$f/$f.fastq.gz | samtools sort -@ 20 -O BAM -o $f.bam - 
samtools index $f.bam 
python /media/txpn/nvme/softwares/flair/bin/bam2Bed12.py -i $f.bam > $f.bed12
#correction
python /media/txpn/nvme/softwares/flair/flair.py correct -q $f.bed12 -g $genome -t 30 -o $f -f /media/txpn/nvme/Databases/hg38/gencode.v32.chr_patch_hapl_scaff.annotation.gtf
echo done: $f

done 

#collapse
cat *_all_corrected.bed > All_Corrected.bed
zcat Sample*/*.fastq.gz | gzip -c > Total.fastq.gz
#Collapsing
python /media/txpn/nvme/softwares/flair/flair.py collapse -g $genome -f /media/txpn/nvme/Databases/hg38/gencode.v32.chr_patch_hapl_scaff.annotation.gtf -t 30 -r Total.fastq.gz -q  All_Corrected.bed

#Quantification
dir="/media/txpn/My_Book_Duo/Omid/Nanopore/Nanopore_34_samples_20210629"
python /media/txpn/nvme/softwares/flair/flair.py quantify -t 30 --tpm -o All_Corrected_counts.txt -r $dir/SampleManifest.txt -i flair.collapse.isoforms.fa
