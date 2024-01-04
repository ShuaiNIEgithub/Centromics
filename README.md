
# Introduction
Centromeres are extremely unique regions of heterochromatin on chromosomes. They typically exhibit close interchromosomal or intrachromosomal interactions in the three-dimensional genome structure, as well as being marked by specific DNA-binding proteins, including the centromere protein H3 variant (CENH3), and are enriched with high-copy tandem repeat sequences. Centromics is a pipeline used to identify centromeres in highly contiguous genome assemblies, such as the T2T genome. It utilizes PacBio/ONT long-read sequencing to identify specific clusters of high-copy tandem repeats, thereby pinpointing potential centromere regions. Centromeres can also be identified based on HI-C sequencing and CENH3 chip-seq.


# Install
```
git clone --recurse-submodules https://github.com/zhangrengang/Centromics
cd Centromics
conda env create -f Centromics.yaml
conda activate RepCent
./install.sh
Centromics -h
```

# Workflow
Let's assume that you want to identify the centromere region from T2T genome assembly of Arabidopsis thaliana (Col-CEN; doi:10.1126/science.abi7489).

```
PWD=`pwd`
mkdir $PWD/example_data
cd $PWD/example_data
```

## 1.Download and process the input data
Here, we prepare two required data (reference genome and Pacbio or ONT long-reads) and two optional data (CENH3 ChIP-seq and Hi-C illumina short-resds).

### reference genome (FASTA format) 
```
download from: https://github.com/schatzlab/Col-CEN/blob/main/v1.2/Col-CEN_v1.2.fasta.gz
zcat Col-CEN_v1.2.fasta.gz > ref.fa
```

### Long-reads data

The amount of data is usually large for long-reads data. High quality data can be used after filtering. Of course you can also use raw data and Centromics will perform a simple quality control.

```
# Pacbio hifi long-reads
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR621/003/ERR6210723/ERR6210723.fastq.gz
mv ERR6210723.fastq.gz hifi.fa

# ONT long-reads
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR553/005/ERR5530735/ERR5530735.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR553/006/ERR5530736/ERR5530736.fastq.gz
mv ERR5530735.fastq.gz ont.1.fq.gz
mv ERR5530735.fastq.gz ont.2.fq.gz
```

### CENH3 ChIP-seq illumina short-resds (optional)

**You can download our example data directly as input data.**

bam format file: https://figshare.com/articles/dataset/chip_bam/21946367

bam.bai format file: https://figshare.com/articles/dataset/chip_bam_bai/21946337


**Or you can download the raw data and analyse it for example files yourself.**
```
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR443/007/SRR4430537/SRR4430537_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR443/007/SRR4430537/SRR4430537_2.fastq.gz
# Mapping the ChIP-seq data against reference genome to obatain the bam format file.
minimap2 ref.fa SRR4430537_1.fastq.gz SRR4430537_2.fastq.gz -x sr -a -t 10 | samtools view -b -F 2304 | samtools sort > chip.bam
```

### Hi-C sequence illumina short-resds** (optional)
**You can download our example .hic format file directly.**

merged_nodups.txt.gz: https://figshare.com/articles/dataset/merged_nodups_txt_gz/21946505

**Or you can download the raw data and analyse it for .hic format file yourself.**

```
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR150/009/SRR1504819/SRR1504819_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR150/009/SRR1504819/SRR1504819_2.fastq.gz
# aligning the Hi-C data to reference genome for .merged_nodups.txt format file by juicer.sh (https://aidenlab.org/documentation.html).
# converting .merged_nodups.txt to .hic format file (https://github.com/aidenlab/Juicebox/blob/master/HiCFormatV8.md).
java -jar juicebox_tools.jar pre -n aligned/merged_nodups.txt merged_nodups.hic ref.fa.chrom.sizes
```

## 2.Perform Centromics 

Different input data can be combined. Here we list all the possibilities.

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

## How to interpret the output files
Here, we take the output results of the 'Pacbio hifi long reads + HiC data + ChIP data' data input combination as an example.

```
cd $PWD/hifihicchip
tree 
.
├── hifihicchip.candidate_peaks.bed # KEY
├── hifihicchip.chip.count          # KEY
├── hifihicchip.circos_legend.pdf
├── hifihicchip.circos_legend.txt
├── hifihicchip.circos.pdf          # KEY
├── hifihicchip.circos.png
├── hifihicchip.hic.count.100000.inter_chr
├── hifihicchip.hic.count.100000.intra_chr
├── hifihicchip.trf.count
├── hifihicchip.trf.fa #Sequences of Candidate TR Monomer
├── hifihicchip.circos #Intermediate PATH for Visualization
│   ├── circos.conf
│   ├── circos.pdf
│   ├── circos.png
│   ├── circos.svg
│   ├── colors.conf
│   ├── data
│   │   ├── chip_density.txt
│   │   ├── genome_gc.txt
│   │   ├── genome_karyotype.txt
│   │   ├── hic_inter.density.txt
│   │   └── tr_density.txt
│   ├── etc
│   │   ├── arial.ttf
│   │   ├── background.black.conf
│   │   ├── background.white.conf
│   │   ├── bands.conf
│   │   ├── BREWER
│   │   ├── brewer.all.conf
│   │   ├── brewer.conf
│   │   ├── colors.brewer.conf
│   │   ├── colors.conf
│   │   ├── colors_fonts_patterns.conf
│   │   ├── colors.hsv.conf
│   │   ├── colors.ucsc.conf
│   │   ├── colors.unix.txt
│   │   ├── fonts.bk.conf
│   │   ├── fonts.conf
│   │   ├── gddiag.conf
│   │   ├── housekeeping.conf
│   │   ├── ideogram.conf
│   │   ├── ideogram.label.conf
│   │   ├── ideogram.position.conf
│   │   ├── image.black.conf
│   │   ├── image.conf
│   │   ├── makehuesteps
│   │   ├── patterns.conf
│   │   ├── patterns.svg.conf
│   │   ├── README
│   │   ├── ticks.conf
│   │   └── tracks
│   │       ├── axis.conf
│   │       ├── connector.conf
│   │       ├── heatmap.conf
│   │       ├── highlight.bg.conf
│   │       ├── highlight.conf
│   │       ├── histogram.conf
│   │       ├── line.conf
│   │       ├── link.conf
│   │       ├── README
│   │       ├── scatter.conf
│   │       ├── text.conf
│   │       └── tile.conf
│   ├── histogram.conf
└── └── image.generic.conf
```



Firstly, we can preliminarily locate the centromere region by examining the distribution of tandem repeat clusters on the circos plot "hifihicchip.circos.pdf". The legends were saved in files "hifihicchip.circos_legend.pdf" and "hifihicchip.circos_legend.txt". The potential peaks were described in "hifihicchip.candidate_peaks.bed", extracted from the Circos plot at different omics levels.

Next, we can determine the start and end positions of potential centromere regions based on files that quantitatively describe the density distribution of specific features. These files include hifihicchip.trf.count, which represents the density of tandem repeats from hifi long reads, hifihicchip.hic.count.100000.inter_chr and hifihicchip.hic.count.100000.intra_chr, which represents the density of Hi-C contact links an interchromosomal and intrachromosomal level, and hifihicchip.chip.count, which represents the density of ChIP reads.

Lastly, taking chromosome 5 (Chr5) as an example:
The potential centromere region locatable based on PacBio HiFi long reads (hifihicchip.trf.count) is 12,330,000 to 13,690,000 bp.
The potential centromere region locatable based on ChIP-seq short reads (hifihicchip.chip.count) is 11,790,000 to 14,560,000 bp.
The potential centromere region locatable based on Hi-C short reads (hifihicchip.hic.count) is 11,830,000 to 14,290,000 bp.



