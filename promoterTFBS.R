# Run AME to find enriched transcripction factor binding sites
# Match AME results (i.e. TFs) to corresponding genes

# Remove objects in workspace
rm(list = ls())

suppressMessages(require(memes))
suppressMessages(require(readr))
suppressMessages(require(dplyr))
suppressMessages(require(tidyr))
suppressMessages(require(rtracklayer))
suppressMessages(library(GenomicRanges))
suppressMessages(require(BSgenome.Hsapiens.UCSC.hg38))

# HOCOMOCO download: https://hocomoco11.autosome.org/downloads_v11
# Downloaded on March 7th 2023
options(meme_db = paste("/Users/achan/Documents/meme-5.5.1/",
		"/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme", sep = ""),
		meme_bin = "/Users/achan/Documents/meme-5.5.1/bin"
)

# ---------------------------------
# mean_meth <- "/u/achan/DW_variantOlap/resu_combp.csv.SIDAK_lt0.05.bed.hg38"

# Local file
mean_meth <- "/Users/achan/Documents/FetalHindbrain_Epigenetics/GRN_inference/Data/DMR_overlap.R/dmrPromoter_pairs.tsv"

# ---------------------------------
# Read in DMRs which overlap with promoters
dmrs <- read.delim(mean_meth, header = T)
message(sprintf("Loaded %i DMRs", nrow(dmrs)))

# Coerce promoter start and end positions into numeric to run IRanges
dmrs$promoter.start <- as.numeric(dmrs$promoter.start) 

dmrs$promoter.end   <- as.numeric(dmrs$promoter.end) 

# Convert to GenomicRanges object
promotersGR <- GRanges(dmrs$promoter.seqnames, IRanges(dmrs$promoter.start, dmrs$promoter.end))
hs.genome <- BSgenome.Hsapiens.UCSC.hg38

# ---------------------------------
# Output directory
outpath <- "~/Documents/FetalHindbrain_Epigenetics/GRN_inference/Data"
outDir <- sprintf("%s/Data/promoterTFBS.R", dirname(outpath))
if(!file.exists(outDir)) dir.create(outDir)

# ---------------------------------
# Run motif enrichment
# Interpret results: https://meme-suite.org/meme/doc/ame-output-format.html#html_results
seq <- get_sequence(promotersGR, hs.genome)
ame <- memes::runAme(seq, outdir = outDir)

# ---------------------------------
# Map enriched transcription factors to genes
inpath <- "~/Documents/FetalHindbrain_Epigenetics/GRN_inference/Data"
inDir  <- sprintf("%s/Data", dirname(inpath))

# ---------------------------------
# Select n most significant motifs overall

# Read in ame.tsv file
ame     <- read_delim(file = sprintf("%s/promoterTFBS.R/ame.tsv", inDir), delim = "\t")
ameDF   <- as.data.frame(ame)

# Pull motif names
ameCpos  <- regexpr("_", ameDF$motif_ID, fixed = TRUE) - 1
ameMotif <- substr(ameDF$motif_ID, 1, ameCpos)

# Select n most significant motifs by rank (https://meme-suite.org/meme/doc/ame-output-format.html#consensus_doc)
ameMotif <- ameMotif[1:100]
# ameMotif <- ameMotif[1:25]

# ---------------------------------
ameSeq    <- read_delim(file = sprintf("%s/promoterTFBS.R/sequences.tsv", inDir), delim = "\t")     # 150 859 instances
ameSeq_DF <- as.data.frame(ameSeq)                                        

# 227 blank rows in sequences.tsv file (created by AME as divider rows)
x <- problems(ameSeq)

# Drop false positives in column "class"
ameSeq_DF <- ameSeq_DF[ameSeq_DF$class == "tp", ]                                                   #  97 746 instances

# Drop NAs (i.e. blank rows from before)
ameSeq_DF <- ameSeq_DF %>% drop_na(seq_ID)                                                          #  97 519 instances

# ---------------------------------
# Select most confident (by PWM score) sequence per motif

# Pull motif names
motifCpos  <- regexpr("_", ameSeq_DF$motif_ID, fixed = TRUE) - 1
ameSeq_DF$motif_ID <- substr(ameSeq_DF$motif_ID, 1, motifCpos)

# Unique motif names
motifNames <- base::unique(ameSeq_DF$motif_ID, incomparables = FALSE)

# Empty data frame to store max PWM score
ameSeq_PWM <- data.frame(matrix(ncol = 3, nrow = NROW(motifNames)))
colnames(ameSeq_PWM) <- c("motif_ID", "seq_ID", "PWM_Score")
ameSeq_PWM$motif_ID <- motifNames

i = 0 

for(each in motifNames){
  # Index
  i = i + 1

  message(i)
  
  maxPWM <- max(ameSeq_DF$PWM_score[ which(ameSeq_DF$motif_ID == each) ])
  
  ameSeq_PWM$seq_ID[i]    <- ameSeq_DF$seq_ID[ which(ameSeq_DF$PWM_score == maxPWM) ]
  ameSeq_PWM$PWM_Score[i] <- format(maxPWM, scientific = FALSE) 
}

# Pull promoter sequence coordinates to match with DMR-promoters
pos   <- ameSeq_PWM$seq_ID
cpos  <- gregexpr(":|-|$", pos)					  # (anchors via https://stringr.tidyverse.org/articles/regular-expressions.html)
cpos1 <- unlist(lapply(cpos, function(x) x[1]))   # get colon
cpos2 <- unlist(lapply(cpos, function(x) x[2]))   # get hyphen (-)
cpos3 <- unlist(lapply(cpos, function(x) x[3]))   # get end of string 
start <- as.numeric(substr(pos, cpos1 + 1, cpos2 - 1))
end   <- as.numeric(substr(pos, cpos2 + 1, cpos3 - 1))

coordTFBS <- GRanges(ameSeq_PWM$motif_ID, IRanges(start, end))
coordTFBS <- as.data.frame(coordTFBS)
colnames(coordTFBS)[colnames(coordTFBS) == "seqnames"] <- "motif_ID"
colnames(coordTFBS)[colnames(coordTFBS) == "start"] <- "promoter.start"
colnames(coordTFBS)[colnames(coordTFBS) == "end"] <- "promoter.end"

# ---------------------------------
# Select n most confident (by PWM score) sequences per motif
# NOTE: Does not take into account duplicates from multiple DMRS
# within same promoter region (consider adding duplicates to edge score?)

# Pull motif names
motifCpos  <- regexpr("_", ameSeq_DF$motif_ID, fixed = TRUE) - 1
ameSeq_DF$motif_ID <- substr(ameSeq_DF$motif_ID, 1, motifCpos)

# Unique motif names to iterate through AME results
motifNames <- base::unique(ameSeq_DF$motif_ID, incomparables = FALSE)

# Empty data frame to store max PWM scores
ameSeq_PWM <- data.frame(matrix(ncol = 3, nrow = 5 * NROW(motifNames)))
colnames(ameSeq_PWM) <- c("motif_ID", "seq_ID", "PWM_Score")
ameSeq_PWM$motif_ID <- motifNames

i = 0 

for(each in motifNames){
  i = i + 1
  
  message(i)
  
  maxPWM <- vector(mode = "double", length = 5)
  
  # Sort, and take 5 highest PWM scores per TF 
  maxPWM <- as.numeric(head(sort(ameSeq_DF$PWM_score[ which(ameSeq_DF$motif_ID == each) ], decreasing = TRUE), 5))
  
  # Append PWM scores to data frame
  ameSeq_PWM$PWM_Score[ which(ameSeq_PWM$motif_ID == each) ] <- format(maxPWM, scientific = FALSE)
  
  k = 1
  
  for(j in 1:length(maxPWM)){
    # Pull associated TF binding sequences for each PWM score 
    ameSeq_PWM$seq_ID[ which(ameSeq_PWM$motif_ID == each) ][j] <- ameSeq_DF$seq_ID[ which(ameSeq_DF$PWM_score == maxPWM[k]) ]
    
    k = k + 1
  }
}

# Pull promoter sequence coordinates to match with DMR-promoters
pos   <- ameSeq_PWM$seq_ID
cpos  <- gregexpr(":|-|$", pos)					  # (anchors via https://stringr.tidyverse.org/articles/regular-expressions.html)
cpos1 <- unlist(lapply(cpos, function(x) x[1]))   # get colon
cpos2 <- unlist(lapply(cpos, function(x) x[2]))   # get hyphen (-)
cpos3 <- unlist(lapply(cpos, function(x) x[3]))   # get end of string 
start <- as.numeric(substr(pos, cpos1 + 1, cpos2 - 1))
end   <- as.numeric(substr(pos, cpos2 + 1, cpos3 - 1))

coordTFBS <- GRanges(ameSeq_PWM$motif_ID, IRanges(start, end))
coordTFBS <- as.data.frame(coordTFBS)
colnames(coordTFBS)[colnames(coordTFBS) == "seqnames"] <- "motif_ID"
colnames(coordTFBS)[colnames(coordTFBS) == "start"] <- "promoter.start"
colnames(coordTFBS)[colnames(coordTFBS) == "end"] <- "promoter.end"

# ---------------------------------
# Pull promoter sequence coordinates to match with DMR-promoters
pos   <- ameSeq_DF$seq_ID
cpos  <- gregexpr(":|-|$", pos)					  # (anchors via https://stringr.tidyverse.org/articles/regular-expressions.html)
cpos1 <- unlist(lapply(cpos, function(x) x[1]))   # get colon
cpos2 <- unlist(lapply(cpos, function(x) x[2]))   # get hyphen (-)
cpos3 <- unlist(lapply(cpos, function(x) x[3]))   # get end of string 
start <- as.numeric(substr(pos, cpos1 + 1, cpos2 - 1))
end   <- as.numeric(substr(pos, cpos2 + 1, cpos3 - 1))

# Get associated TF binding motif IDs
motifPos <- ameSeq_DF$motif_ID
mpos     <- gregexpr("_", motifPos)
mpos1    <- unlist(lapply(mpos, function(x) x[1]))
motif    <- substr(motifPos, 1, mpos1 - 1)

coordTFBS <- GRanges(motif, IRanges(start, end))
coordTFBS <- as.data.frame(coordTFBS)
colnames(coordTFBS)[colnames(coordTFBS) == "seqnames"] <- "motif_ID"
colnames(coordTFBS)[colnames(coordTFBS) == "start"] <- "promoter.start"
colnames(coordTFBS)[colnames(coordTFBS) == "end"] <- "promoter.end"

# Pair transcription factors with associated genes
geneTFs <- inner_join(coordTFBS, dmrs, by = c("promoter.start", "promoter.end"), multiple = "all")
geneTFs_names <- subset(geneTFs, select = c("motif_ID", "gene_name"))                                    

# Remove duplicates (from multiple DMR overlaps within same promoter region e.g., promoter for PLA2G2C)
geneTFs_names <- geneTFs_names %>% distinct(gene_name, motif_ID, .keep_all = TRUE)                       

# Subset for most significant motifs
geneTFs_names <- subset(geneTFs_names, motif_ID %in% ameMotif)                                           

# Significant motif IDs
sigMotifs <- base::unique(geneTFs_names$motif_ID, incomparables = FALSE)

# Output for promoterTSS_toNetwork.R
write.table(geneTFs_names,
			file = sprintf("%s/geneTFs_names.tsv", outDir),
			sep = "\t", row.names = F, col.names = T, quote = F)