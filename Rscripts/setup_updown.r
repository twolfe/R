# ec$DAVID   # terms <- read.delim("GoOut_all_test3.txt")
# EC$genelist # genes <- read.delim("FDR005_sp4_test3_forGO.txt", sep="")

# TERMS = ec$DAVID
# Category         ID                             Term
## 1       BP GO:0007507                heart development
## 2       BP GO:0001944          vasculature development
## 3       BP GO:0001568         blood vessel development
## 4       BP GO:0048729             tissue morphogenesis
## 5       BP GO:0048514       blood vessel morphogenesis
## 6       BP GO:0051336 regulation of hydrolase activity
## 
## Genes == number of DE genes, not the total number of genes from the ref genome ! 
## 1       DLC1, NRP2, NRP1, EDN1, PDLIM3, GJA1, TTN, GJA5, ZIC3, TGFB2, CERKL, GATA6, COL4A3BP, GAB1, SEMA3C, MKL2, SLC22A5, MB, PTPRJ, RXRA, VANGL2, MYH6, TNNT2, HHEX, MURC, MIB1, FOXC2, FOXC1, ADAM19, MYL2, TCAP, EGLN1, SOX9, ITGB1, CHD7, HEXIM1, PKD2, NFATC4, PCSK5, ACTC1, TGFBR2, NF1, HSPG2, SMAD3, TBX1, TNNI3, CSRP3, FOXP1, KCNJ8, PLN, TSC2, ATP6V0A1, TGFBR3, HDAC9
## 2 GNA13, ACVRL1, NRP1, PGF, IL18, LEPR, EDN1, GJA1, FOXO1, GJA5, TGFB2, WARS, CERKL, APOE, CXCR4, ANG, SEMA3C, NOS2, MKL2, FGF2, RAPGEF1, PTPRJ, RECK, EFNB2, VASH1, PNPLA6, THY1, MIB1, NUS1, FOXC2, FOXC1, CAV1, CDH2, MEIS1, WT1, CDH5, PTK2, FBXW8, CHD7, PLCD1, PLXND1, FIGF, PPAP2B, MAP2K1, TBX4, TGFBR2, NF1, TBX1, TNNI3, LAMA4, MEOX2, ECSCR, HBEGF, AMOT, TGFBR3, HDAC7
## 3        GNA13, ACVRL1, NRP1, PGF, IL18, LEPR, EDN1, GJA1, FOXO1, GJA5, TGFB2, WARS, CERKL, APOE, CXCR4, ANG, SEMA3C, NOS2, MKL2, FGF2, RAPGEF1, PTPRJ, RECK, VASH1, PNPLA6, THY1, MIB1, NUS1, FOXC2, FOXC1, CAV1, CDH2, MEIS1, WT1, CDH5, PTK2, FBXW8, CHD7, PLCD1, PLXND1, FIGF, PPAP2B, MAP2K1, TBX4, TGFBR2, NF1, TBX1, TNNI3, LAMA4, MEOX2, ECSCR, HBEGF, AMOT, TGFBR3, HDAC7
## 4                                   DLC1, ENAH, NRP1, PGF, ZIC2, TGFB2, CD44, ILK, SEMA3C, RET, AR, RXRA, VANGL2, LEF1, TNNT2, HHEX, MIB1, NCOA3, FOXC2, FOXC1, TGFB1I1, WNT5A, COBL, BBS4, FGFR3, TNC, BMPR2, CTNND1, EGLN1, NR3C1, SOX9, TCF7L1, IGF1R, FOXQ1, MACF1, HOXA5, BCL2, PLXND1, CAR2, ACTC1, TBX4, SMAD3, FZD3, SHANK3, FZD6, HOXB4, FREM2, TSC2, ZIC5, TGFBR3, APAF1
## 5                                                                                            GNA13, CAV1, ACVRL1, NRP1, PGF, IL18, LEPR, EDN1, GJA1, CDH2, MEIS1, WT1, TGFB2, WARS, PTK2, CERKL, APOE, CXCR4, ANG, SEMA3C, PLCD1, NOS2, MKL2, PLXND1, FIGF, FGF2, PTPRJ, TGFBR2, TBX4, NF1, TBX1, TNNI3, PNPLA6, VASH1, THY1, NUS1, MEOX2, ECSCR, AMOT, HBEGF, FOXC2, FOXC1, HDAC7
## 6                                                                               CAV1, XIAP, AGFG1, ADORA2A, TNNC1, TBC1D9, LEPR, ABHD5, EDN1, ASAP2, ASAP3, SMAP1, TBC1D12, ANG, TBC1D14, MTCH1, TBC1D13, TBC1D4, TBC1D30, DHCR24, HIP1, VAV3, NOS1, NF1, MYH6, RICTOR, TBC1D22A, THY1, PLCE1, RNF7, NDEL1, CHML, IFT57, ACAP2, TSC2, ERN1, APAF1, ARAP3, ARAP2, ARAP1, HTR2A, F2R
##      adj_pval
## 1 0.000002170
## 2 0.000010400
## 3 0.000007620
## 4 0.000119000
## 5 0.000720000
## 6 0.001171166

# GENES = EC$genelist
##        ID    logFC   AveExpr        t  P.Value adj.P.Val        B
## 1 Slco1a4 6.645388 1.2168670 88.65515 1.32e-18  2.73e-14 29.02715
## 2 Slc19a3 6.281525 1.1600468 69.95094 2.41e-17  2.49e-13 27.62917
## 3     Ddc 4.483338 0.8365231 65.57836 5.31e-17  3.65e-13 27.18476
## 4 Slco1c1 6.469384 1.3558865 59.87613 1.62e-16  8.34e-13 26.51242
## 5  Sema3c 5.515630 2.3252117 58.53141 2.14e-16  8.81e-13 26.33626
## 6 Slc38a3 4.761755 0.9218670 54.11559 5.58e-16  1.76e-12 25.70308




### create an intermediary files that contain all the GO terms and all the genes that were used in the DE analysis (gene universe)
#setwd("~/Documents/rnaseq-bromeliads/GOplot/Aco_total_GO/")
# fit the A.comosus tables (GO terms and candidates description)
# total genes list with GO terms
#tot.2 <- read.delim("Aco_TotalGenes_20170406.csv", sep="\t")
# total genes list with description and candidates 
#Aco <- read.delim("Aco_list.txt", sep="")
#aco.merge <- merge(tot.2, Aco, by = "genes", all = T)
# create one line / GO term : 
### split function !!!

################ TOPGO ####################################################
library("topGO", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3")
library(ALL)
data(ALL)
library(xtable)

setwd("/data/phdData/orchis/counts/de.featureCounts/enrichments")

dat.filtered <- read.table("/home/botanik/Documents/phd/research/environmental_study/deEnviroCorrelations/2018/datFilteredCounts.txt", header = T)
geneID2GO <- readMappings(file = "/data/phdData/orchis/Orchis.transcriptome.go.annotation.topGO.map")
geneID2GO <- readMappings(file = "/data/phdData/orchis/counts/de.featureCounts/retained.after.filtering.go.map.txt")
GO2geneID <- inverseList(geneID2GO)
geneNames <- names(geneID2GO)
head(geneNames)
#resEdgeRtopTags is the table found after running cleanDE.R
#myInterestingGenes <- rownames(resEdgeRtopTags$table[resEdgeRtopTags$table$FDR <= 0.05,])
#write.table(myInterestingGenes, "/data/phdData/orchis/counts/de.featureCounts/both.transcripts.de.txt", sep="\t", quote = F)#to modify some more for compatibility
myInterestingGenes <- read.table("/data/phdData/orchis/counts/de.featureCounts/both.de.results.txt", header = T)
#myInterestingGenes <- read.table("/data/phdData/orchis/hylite/HyLiTE_analysis/genesWithDiagnosticSNP.txt")
myInterestingGenes <- as.vector(myInterestingGenes$V1)
geneID2GO.de <- readMappings(file = "/data/phdData/orchis/counts/de.featureCounts/both.transcripts.de.topGO.map")
#geneID2GO.de <- readMappings(file = "/data/phdData/orchis/hylite/HyLiTE_analysis/topGO.diagnostic.genes.map")
GO2geneID.de <- inverseList(geneID2GO.de)

geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames
str(geneList)

GOdata.MF <- new("topGOdata", ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata.BP <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata.CC <- new("topGOdata", ontology = "CC", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)

#TESTS
#resultFis.BP <- runTest(GOdata.BP, algorithm = "classic", statistic = "fisher")
#resultKS.BP <- runTest(GOdata.BP, algorithm = "elim", statistic = "fisher")
#resultWeight.BP <- runTest(GOdata.BP, algorithm = "weight", statistic = "fisher")
resultWeight01.BP <- runTest(GOdata.BP, statistic = "fisher")
resultWeight01.MF <- runTest(GOdata.MF, statistic = "fisher")
resultWeight01.CC <- runTest(GOdata.CC, statistic = "fisher")

allRes.BP <- GenTable(GOdata.BP, weight01_pval=resultWeight01.BP, orderBy = "weight01", ranksOf = "weight01",topNodes = 100)
allRes.BP <- cbind(allRes.BP,"BP")
colnames(allRes.BP) <- c("GO.ID","Term","Annotated","Significant","Expected","weight01_pval", "branch")

allRes.MF <- GenTable(GOdata.MF, weight01_pval=resultWeight01.MF, orderBy = "weight01", ranksOf = "weight01",topNodes = 100)
allRes.MF <- cbind(allRes.MF,"MF")
colnames(allRes.MF) <- c("GO.ID","Term","Annotated","Significant","Expected","weight01_pval", "branch")

allRes.CC <- GenTable(GOdata.CC, weight01_pval=resultWeight01.CC, orderBy = "weight01", ranksOf = "weight01",topNodes = 100)
allRes.CC <- cbind(allRes.CC,"CC")
colnames(allRes.CC) <- c("GO.ID","Term","Annotated","Significant","Expected","weight01_pval", "branch")

allRes <- rbind(allRes.BP,allRes.MF)
allRes <- rbind(allRes, allRes.CC)

head(allRes)
x.big <- xtable(allRes, caption = "A \\code{longtable} spanning several pages")
print(x.big, hline.after=c(-1, 0), tabular.environment = "longtable", scalebox = 0.7)
term.genes <- genes(GOdata, allRes$GO.ID)

### Extract significant results 
allRes.sign <- allRes[allRes$weight01_pval <= format(0.05, scientific = TRUE),]

##### link here
allGO.BP = genesInTerm(GOdata.BP)
allGO.MF = genesInTerm(GOdata.MF)
allGO.CC = genesInTerm(GOdata.CC)
allGO = c(allGO.BP, allGO.MF, allGO.CC)

SAM_ANOTATION = lapply(allGO,function(x) x[x %in% myInterestingGenes] )
enriched_go_with_my_genes <- lapply(SAM_ANOTATION[allRes.sign[,1]], paste0, collapse = ", ")

enriched_go_with_my_genes.list <- c()
for (i in 1:length(enriched_go_with_my_genes)){
  enriched_go_with_my_genes.list <- c(enriched_go_with_my_genes.list, enriched_go_with_my_genes[[i]])
}
############### GOPLOT #################################################
library(splitstackshape)
library("GOplot", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3")

########################################################################
### very important: change this in the function before running it : ####
########################################################################
circle_dat <- function(terms, genes){
  
  colnames(terms) <- tolower(colnames(terms))
  terms$genes <- toupper(terms$genes)
  genes$ID <- toupper(genes$ID)
  tgenes <- strsplit(as.vector(terms$genes), ', ')
  if (length(tgenes[[1]]) == 1) tgenes <- strsplit(as.vector(terms$genes), ',')
  count <- sapply(1:length(tgenes), function(x) length(tgenes[[x]]))
  logFC <- sapply(unlist(tgenes), function(x) genes$logFC[match(x, genes$ID)])
  if(class(logFC) == 'factor'){
    logFC <- gsub(",", ".", gsub("\\.", "", logFC))
    logFC <- as.numeric(logFC)
  }
  s <- 1; zsc <- c()
  for (c in 1:length(count)){
    value <- 0
    e <- s + count[c] - 1
    value <- logFC[s:e]
    #value <- sapply(logFC[s:e], function(x) ifelse(x > 0, 1, -1))
    zsc <- c(zsc, sum(value, na.rm = F) / sqrt(count[c])) #### HERE : na.rm = TRUE, takes all the genes into account !
    s <- e + 1
  }
  if (is.null(terms$id)){
    df <- data.frame(category = rep(as.character(terms$category), count), term = rep(as.character(terms$term), count),
                     count = rep(count, count), genes = as.character(unlist(tgenes)), logFC = logFC, adj_pval = rep(terms$adj_pval, count),
                     zscore = rep(zsc, count), stringsAsFactors = FALSE)
  }else{
    df <- data.frame(category = rep(as.character(terms$category), count), ID = rep(as.character(terms$id), count), term = rep(as.character(terms$term), count),
                     count = rep(count, count), genes = as.character(unlist(tgenes)), logFC = logFC, adj_pval = rep(terms$adj_pval, count),
                     zscore = rep(zsc, count), stringsAsFactors = FALSE)
  }
  return(df)
}

both.go <- read.table("/data/phdData/orchis/counts/de.featureCounts/enrichments/both.blast2go.enrichments.txt"
                      , header = T, sep = "\t")

#go.dataframe <- data.frame("Category" = both.go$GO.Category, "ID" = both.go$GO.ID, "Term" = both.go$GO.Name, "Genes" = both.go$TestSet.Sequences, "adj_pval" = both.go$FDR)
go.dataframe <- data.frame("Category" = allRes.sign$branch, "ID" = allRes.sign$GO.ID, "Term" = allRes.sign$Term, "Genes" = as.vector(enriched_go_with_my_genes.list), "adj_pval" = as.numeric(allRes.sign$weight01_pval))
#EC.genelist <- data.frame("ID" = rownames(resEdgeRtopTags$table), "logFC" = resEdgeRtopTags$table$logFC, "AveExpr" = rowMeans(dat.filtered[rownames(resEdgeRtopTags$table),]), "P.Value" = resEdgeRtopTags$table$PValue, "adj.P.Val" = resEdgeRtopTags$table$FDR)
EC.genelist <- data.frame("ID" = rownames(myInterestingGenes), "logFC" = myInterestingGenes$logFC, "AveExpr" = rowMeans(dat.filtered[rownames(myInterestingGenes),]), "P.Value" = myInterestingGenes$PValue, "adj.P.Val" = myInterestingGenes$FDR)
resEdgeRtopTags$table[myInterestingGenes,]

circ <- circle_dat(go.dataframe, EC.genelist)
GOplot::GOBar(circ, zsc.col = c("blue", "white", "red"), display = 'multiple')

GOplot::GOCircle(subset(circ, category == 'BP' | category == 'MF'), lfc.col = c("blue", "red"))

# Interesting categorie for publication
IDs <- c("GO:0015979", "GO:0015994", "GO:0015996", "GO:0009963", "GO:0009873")
GOplot::GOCircle(circ, lfc.col = c("blue", "red"), nsub = IDs)

circ.BP.CC.MF <- rbind(subset(circ, category == "BP"), subset(circ, category == "CC"), subset(circ, category == "MF"))

reduced_circ <- reduce_overlap(circ.BP.CC.MF, overlap = 0.75)
GOBubble(reduced_circ, labels = 4.5, display = 'multiple', table.legend = T)
reduced_circ[c(1,2,3,4,5,6,7,24,36,38,44,45,54,55,56),]

write.table(circ.day, "circ.day.txt", sep = "\t")
write.table(circ.night, "circ.night.txt", sep = "\t")