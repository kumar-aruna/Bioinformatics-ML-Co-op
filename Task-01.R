# Ensure necessary libraries are installed and loaded
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# Installing GenomicFeatures and S4Vectors from Bioconductor
BiocManager::install(c("GenomicFeatures", "S4Vectors"))
library(GenomicFeatures)
library(S4Vectors)



# Step 01: We are downloading gene annotation from the GENCODE database.
# GENCODE: provides comprehensive gene annotations. 
# For this task, we are looking for GTF file for the comprehensive gene annotation
# of the primary assembly (which includes chromosomes and scaffolding) from release 45.
# Note: Release 45 is the specific version of the gene annotation data provided by the GENCODE project.


# Step 01: We are downloading gene annotation from the GENCODE database

# Setting a longer timeout because we're about to download a big file.
options(timeout = max(1000, getOption("timeout")))

# Function to download the gtf file from given url
download_gtf <- function(url, destination) {
  # Parameters:
  # url:Type:String, The web address from where the GTF file should be downloaded. 
  # destination: The local path where the downloaded GTF file should be saved. 
  
  # The function uses download.file() from R's base package, which downloads files over the internet.
  # mode: is set to "wb" to write the file in binary mode, ensuring the file's integrity.
  download.file(url = url, destfile = destination, mode = "wb")
}
  
# Example usage with the corrected URL (removed the leading space)
gencode_url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.primary_assembly.annotation.gtf.gz"
gtf_path <- "Bioinformatics-ML-Co-op/gencode.v45.annotation.gtf.gz"

# Call the corrected function
download_gtf(gencode_url, gtf_path)




# Step 2: Create a TxDb Object

# The TxDb object is a database that stores transcript-level annotations, allowing 
# for efficient querying and analysis.

# First we need to download all required package

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("GenomicFeatures", force = TRUE)
library(GenomicFeatures)
exists("makeTxDbFromGFF", where = asNamespace("GenomicFeatures"))
# Creating a TxDb object from the GTF file
txdb <- GenomicFeatures::makeTxDbFromGFF(file = gtf_path, format = "gtf")
print(txdb)# Let's see what we've got!
# We converted the GTF file into a TxDb object, because TxDb object provides a structured and efficient 
# format for accessing gene and transcript annotations to do downsteem analysis.

# Result: > print(txdb)
# TxDb object:
  # Db type: TxDb
  # Supporting package: GenomicFeatures
  # Data source: Bioinformatics-ML-Co-op/gencode.v45.annotation.gtf.gz
  # Organism: NA
  # Taxonomy ID: NA
  # miRBase build ID: NA
  # Genome: NA
  # Nb of transcripts: 252989
  # Db created by: GenomicFeatures package from Bioconductor
  # Creation time: 2024-03-29 11:49:03 -0400 (Fri, 29 Mar 2024)
  # GenomicFeatures version at creation time: 1.54.4
# RSQLite version at creation time: 2.3.5
# DBSCHEMAVERSION: 1.2



# Step 3: Build an S4 Object Containing All Genes and Transcripts, because it's neat
genes <- genes(txdb)
transcripts <- transcriptsBy(txdb, by = "gene")

gene_transcript_association <- new("list", genes = genes, transcripts = transcripts)
print(gene_transcript_association)
# The main purpose is to organize gene and transcript data into an easily accessible format.


# Step 4: Compute Statistics and Make a Histogram
library(GenomicRanges)
# Extract genes and transcripts

genes <- genes(txdb)

# Extract transcripts and group them by gene
transcripts_by_gene <- transcriptsBy(txdb, by = "gene")


# Define the S4 class for gene-transcript association
# Organizing Genes and Transcripts Neatly
setClass("GeneTranscriptAssociation",
         slots = c(genes = "GRanges", transcripts = "GRangesList")) # It is like a folder where it organize the 
# the list of transcript and genes neatly
gene_transcript_association <- new("GeneTranscriptAssociation", 
                                   genes = genes, 
                                   transcripts = transcripts_by_gene)
print(gene_transcript_association)# Checking output

# Counting in each folder
transcript_counts <- sapply(gene_transcript_association@transcripts, length)
print(transcript_counts)

# Calculate the mean, minimum, and maximum number of transcripts per gene
mean_transcripts <- mean(transcript_counts)
min_transcripts <- min(transcript_counts)
max_transcripts <- max(transcript_counts)

cat("Mean number of transcripts per gene:", mean_transcripts, "\n")
# Mean number of transcripts per gene: 4.000395 

cat("Minimum number of transcripts per gene:", min_transcripts, "\n")
# Minimum number of transcripts per gene: 1 

cat("Maximum number of transcripts per gene:", max_transcripts, "\n")
# Maximum number of transcripts per gene: 296 


# Generate and display a histogram of the number of transcripts per gene
hist(transcript_counts, 
     breaks=seq(from=min_transcripts, to=max_transcripts, by=1), 
     main="Histogram of Transcripts Per Gene", 
     xlab="Number of Transcripts", 
     ylab="Frequency", 
     border="black", 
     col="skyblue")


# Step 05: 
# Save the S4 object for future use
saveRDS(object = gene_transcript_association, file = "~/gene_transcript_association.rds")
