# Convert promoter regulon - TSS linkages to GRNs for Cytoscape visualization

rm(list = ls())

write2Network <- function(inFile){
    geneTFs_names <- read.delim(inFile)
    
    # Write regulons
    regulon <- geneTFs_names[ , c("motif_ID", "gene_name")]
    regulon$edge <- 1
    
    write.table(regulon[ , c(1,3,2)], file = sprintf("%s.sif", inFile), 
        sep = "\t", col = F, row = F, quote = F) 
    
    # Node attributes
    nodes <- union(geneTFs_names$gene_name, geneTFs_names$motif_ID)
    nodes <- unique(nodes)
    isTF  <- rep("gene_name", length(nodes))
    isTF[which(nodes %in% geneTFs_names$motif_ID)] <- "motif_ID"
    nodeAttrs <- data.frame(node = nodes, TF = isTF)

    write.table(nodeAttrs, file = sprintf("%s.nodeAttrs.txt", inFile), 
        sep = "\t", col = T, row = F, quote = F)

    return(list(geneTFs_names = geneTFs_names, regulon = regulon, nodeAttrs = nodeAttrs))
}

# ---------------------------------
inpath <- "/Users/achan/Documents/FetalHindbrain_Epigenetics/GRN_inference/Data"
inDir  <- sprintf("%s/Data/promoterTFBS.R", dirname(inpath))
 
genePromoter_TSSlinkages <- write2Network(sprintf("%s/geneTFs_names.tsv", inDir))

# AME Hit Transcription Factor (TRUE/FALSE)
cytoscapeNodes <- read_excel("Documents/FetalHindbrain_Epigenetics/GRN_inference/Data/promoterTFBS.R/cytoscapeNodes.xlsx")
cytoscapeNodes$ameTF <- "NA"
cytoscapeNodes$ameTF[which(cytoscapeNodes$Nodes %in% geneTFs_names$motif_ID)] <- "TRUE"
cytoscapeNodes$ameTF[which(cytoscapeNodes$ameTF == "NA")] <- "FALSE"

# ---------------------------------
# Known Human TF (http://humantfs.ccbr.utoronto.ca/index.php)
# Downloaded April 14, 2023
humanTF <- read_csv("Documents/FetalHindbrain_Epigenetics/GRN_inference/Data/dmrTFBS.R/DatabaseExtract_v_1.01.csv")
humanTF_compare <- humanTF[c('HGNC symbol', 'Is TF?')]
humanTF_compare <- subset(humanTF_compare, humanTF_compare$`Is TF?` == "Yes")

cytoscapeNodes$isTF <- "NA"
cytoscapeNodes$isTF[which(cytoscapeNodes$Nodes %in% humanTF_compare$`HGNC symbol`)] <- "TRUE"
cytoscapeNodes$isTF[which(cytoscapeNodes$isTF == "NA")] <- "FALSE"

# ---------------------------------
# Known to be mutated/over expressed ("dysregulated" in Group 3/4 MB)
mbGenes <- read_excel("Documents/FetalHindbrain_Epigenetics/GRN_inference/Data/dmrTFBS.R/MB_Genes.xlsx")

cytoscapeNodes$mbGene <- "NA"
cytoscapeNodes$mbGene[which(cytoscapeNodes$Nodes %in% mbGenes$Genes)] <- "TRUE"
cytoscapeNodes$mbGene[which(cytoscapeNodes$mbGene == "NA")] <- "FALSE"

# ---------------------------------
# Known epigenetic regulator
epiRegulator <- read_csv("Documents/FetalHindbrain_Epigenetics/GRN_inference/Data/dmrTFBS.R/Proteins.csv")

cytoscapeNodes$epiReg <- "NA"
cytoscapeNodes$epiReg[which(cytoscapeNodes$Nodes %in% epiRegulator$HGNC_symbol)] <- "TRUE"
cytoscapeNodes$epiReg[which(cytoscapeNodes$epiReg == "NA")] <- "FALSE"

# Import table as Node attributes in Cytoscape
write.table(cytoscapeNodes, 
            file = "/Users/achan/Documents/FetalHindbrain_Epigenetics/GRN_inference/Data/promoterTFBS.R/cytoscapeNodes.tsv", 
            sep = "\t", col = T, row = F, quote = F)