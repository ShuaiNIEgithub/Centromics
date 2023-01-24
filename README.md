
###Install
```
git clone --recurse-submodules https://github.com/zhangrengang/Centromics
cd Centromics
conda env create -f Centromics.yaml
conda activate RepCent
./install.sh
Centromics -h
```

###Workflow
# Let's assume that you want to identify the centromere region from T2T genome assembly of Arabidopsis thaliana (Col-CEN; doi:10.1126/science.abi7489).
````
PWD=`pwd`
mkdir $PWD/example_data
cd $PWD/example_data
```
##1.Download and process the input data
Here, we prepare two required data (reference genome and Pacbio or ONT long-reads) and two optional data (CENH3 ChIP-seq and Hi-C illumina short-resds).

#reference genome (FASTA format) 
```
wget -c https://github.com/schatzlab/Col-CEN/blob/main/v1.2/Col-CEN_v1.2.fasta.gz
zcat Col-CEN_v1.2.fasta.gz > ref.fa
```
#Long-reads
The amount of data is usually large for long-reads data. High quality data can be used after filtering. Of course you can also use raw data and Centromics will perform a simple quality control.
```
#Pacbio hifi long-reads
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR621/003/ERR6210723/ERR6210723.fastq.gz
mv ERR6210723.fastq.gz hifi.fa

#ONT long-reads
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR553/005/ERR5530735/ERR5530735.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR553/006/ERR5530736/ERR5530736.fastq.gz
mv ERR5530735.fastq.gz ont.1.fq.gz
mv ERR5530735.fastq.gz ont.2.fq.gz
```
#CENH3 ChIP-seq illumina short-resds. (optional)
you can download our example data directly as input data.
**bam format file**: 
**bam.bai format file**: https://figshare.com/articles/dataset/chip_bam_bai/21946337

Or you can download the raw data and analyse it for example files yourself.
```
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR443/007/SRR4430537/SRR4430537_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR443/007/SRR4430537/SRR4430537_2.fastq.gz
# Mapping the ChIP-seq data against reference genome to obatain the bam format file.
minimap2 ref.fa SRR4430537_1.fastq.gz SRR4430537_2.fastq.gz -x sr -a -t 10 | samtools view -b -F 2304 | samtools sort > chip.bam
```

#Hi-C sequence illumina short-resds. 
you can download our example .hic format file directly.


#Or you can download the raw data and analyse it for .hic format file yourself.
```
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR150/009/SRR1504819/SRR1504819_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR150/009/SRR1504819/SRR1504819_2.fastq.gz
# aligning the Hi-C data to reference genome for .merged_nodups.txt format file
juicer.sh
# converting .merged_nodups.txt to .hic format file (https://github.com/aidenlab/Juicebox/blob/master/HiCFormatV8.md).
java -jar juicebox_tools.jar pre -n aligned/merged_nodups.txt merged_nodups.hic ref.fa.chrom.sizes
```

##2.Perform Centromics 
#Different input data can be combined. Here we list all the possibilities.
```
cd $PWD
data=$PWD/example_data

# Pacbio hifi long reads
centromics -l $data/hifi.fq.gz -g $data/ref.fa -pre hifi -outdir hifi -tmpdir hifi.tmp -ncpu 10  &>hifi.log

# ONT long reads
centromics -l $data/ont*.fq.gz -g $data/ref.fa -pre ont  -outdir ont -tmpdir ont.tmp -ncpu 10  &>ont.log

# long reads + HiC data
centromics -l $data/hifi.fq.gz -g $data/ref.fa -pre hifihic -hic $data/merged_nodups.txt.gz -outdir hifihic -tmpdir hifihic.tmp -ncpu 10  &>hifihic.log
centromics -l $data/ont*.fq.gz -g $data/ref.fa -pre onthic -hic $data/merged_nodups.txt.gz -outdir onthic -tmpdir onthic.tmp -ncpu 10  &>onthic.log

# long reads + ChIP data
centromics -l $data/hifi.fq.gz -g $data/ref.fa -pre hifichip -chip $data/chip.bam  -outdir hifichip -tmpdir hifichip.tmp -ncpu 10  &>hifichip.log
centromics -l $data/ont*.fq.gz -g $data/ref.fa -pre ontchip  -chip $data/chip.bam  -outdir ontchip -tmpdir ontchip.tmp -ncpu 10  &>ontchip.log

# long reads + HiC data + ChIP data
centromics -l $data/hifi.fq.gz -g $data/ref.fa -pre hifihicchip -chip $data/chip.bam -hic $data/merged_nodups.txt.gz -outdir hifihicchip -tmpdir hifihicchip.tmp -ncpu 10 &>hifihicchip.log
centromics -l $data/ont*.fq.gz -g $data/ref.fa -pre onthicchip  -chip $data/chip.bam -hic $data/merged_nodups.txt.gz -outdir onthicchip -tmpdir onthicchip.tmp -ncpu 10 &>onthicchip.log
```
