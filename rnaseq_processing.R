library(QuasR); library(Rsubread)

setwd('C:/Users/mcgma/Desktop/felix_seqs')

## adapter removal and trimming
# obtain a list of fastq file paths
fastqFiles <- list.files(path = ".", pattern = 'fastq.gz')

# defined processed fastq file names
outfiles <- paste("processed_", fastqFiles,".fastq",sep="")

# process fastq files
# trim 4 bases from the end of the reads (truncateEndBases)
# Remove TGGAATTCTCGGGTGCCAAGG patern if it occurs at the start (Lpattern)
# remove reads shorter than 23 base-pairs (minLength)
preprocessReads(fastqFiles, outfiles, 
                truncateEndBases=4,
                Lpattern="TGGAATTCTCGGGTGCCAAGG",
                minLength=23)

processed <- list.files(path = ".", pattern = '.fastq$')

ref <- system.file("extdata","reference.fa",package="Rsubread")
buildindex(basename="reference_index",reference=ref)

align.stat <- align(index = "./reference_index", annot.inbuilt = 'hg38', readfile1 = processed, 
                    output_file =  paste(processed,"subread",output_format='BAM',sep="."), phredOffset = 64)

files <- list.files(path = ".", pattern = '.BAM$')
fc <- featureCounts(files, annot.inbuilt="hg38")
