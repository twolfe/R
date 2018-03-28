library("adegenet", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3")
library("pegas", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3")

obj <- read.genetix("/data/phdData/HyLiTE/2017.03.08/vcf/all.HyLiTE.filtered.small.call.gtx")
obj <- read.genepop("/data/phdData/HyLiTE/2017.03.08/vcf/all.HyLiTE.filtered.small.gen")
obj <- read.genetix("/data/phdData/HyLiTE/2017.03.08/2017.27.07/call.filteredDP100MQ18.gtx")

#vcf stats
info <- VCFloci("/data/phdData/HyLiTE/2017.03.08/2017.27.07/call.filteredDP100MQ18.vcf")
SNP <- is.snp(info)
table(SNP)

object.size(obj)

obj2 <- new("genlight", list("tB1830_1"=tab(obj)[1,], "tB1830_2"=tab(obj)[2,], "tB1833_1"=tab(obj)[3,], "tB1833_2"=tab(obj)[4,],
            "mA1568"=tab(obj)[5,], "mA1573"=tab(obj)[6,], "mA1661"=tab(obj)[7,], "mA1775"=tab(obj)[8,],
            "mP1722_1"=tab(obj)[9,], "mP1722_2"=tab(obj)[10,], "mP1744"=tab(obj)[11,],
            "mS1757"=tab(obj)[12,], "mS1765_1"=tab(obj)[13,], "mS1765_2"=tab(obj)[14,],
            "tA1553"=tab(obj)[15,], "tA1641"=tab(obj)[16,], "tA1670"=tab(obj)[17,],
            "tB1798_1"=tab(obj)[18,], "tB1798_2"=tab(obj)[19,], "tB1805_1"=tab(obj)[20,], "tB1805_2"=tab(obj)[21,], "tB1812_1"=tab(obj)[22,], "tB1812_2"=tab(obj)[23,],
            "tS1901"=tab(obj)[24,], "tS1902"=tab(obj)[25,], "tS1920_1"=tab(obj)[26,], "tS1920_2"=tab(obj)[27,],
            "fB1804"=tab(obj)[28,], "fB1855"=tab(obj)[29,], "fP1707_1"=tab(obj)[30,], "fP1707_2"=tab(obj)[31,], "fP9375"=tab(obj)[32,],
            "iA1586"=tab(obj)[33,], "iB1176"=tab(obj)[34,], "iB1870"=tab(obj)[35,], "iS1904"=tab(obj)[36,], "iS1908"=tab(obj)[37,]))
ploidy(obj2) <- c(rep(4,27), rep(2,10))
locNames(obj2) <- colnames(tab(obj))
pop(obj2) <- c(rep("t",4), rep("m",10), rep("t",13), rep("f",5), rep("i",5))

#OBJ cleaning
toRemove <- is.na(glMean(obj2, alleleAsUnit = FALSE)) # TRUE where NA
which(toRemove)
obj2.naRemoved <- obj2[,!toRemove]

#Allele frequencies
myFreq <- glMean(obj2.naRemoved)
hist(myFreq, proba=TRUE, col="gold", xlab="Allele frequencies",
     main="Distribution of (second) allele frequencies")

#PCA
pca <- glPca(obj2.naRemoved)
scatter(pca, posi="bottomright")
myCol <- colorplot(pca$scores,pca$scores, transp=TRUE, cex=4)
abline(h=0,v=0, col="grey")

#NJ tree
tre <- nj(dist(as.matrix(obj2.naRemoved)))
plot(tre, typ="fan", show.tip=T)
tiplabels(pch=20, col=myCol, cex=4)

#DAPC
dapc <- dapc(obj2.naRemoved, n.pca=10, n.da=1)
