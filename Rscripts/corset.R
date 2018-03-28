# ----+ LIBRARIES
library("gplots", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3")
library("RUVSeq")
library("edgeR")
library("HTSFilter")
library("RColorBrewer")
library("adegenet")
library("statmod")
library("cgwtools")                                        #---- -

#red <- transp("red",.4)
red = rgb(255, 0, 0, max = 255, alpha = 125, names = "red")
black <- rgb(0, 0, 0, max = 255, alpha = 125, names = "black")
## GLOBAL VARIABLES & FUNCTION FOR DATA ORGANISATION
  #data = read.table("/home/botanik/Documents/phd/research/differential.expression/2016.11/withExtra2RNAlibraries/without.mP1748.mA1567/modelsContrasts/trimmed.counts.txt/trdata.txt", header = T, sep = " ")
  data = read.table("/home/botanik/Documents/phd/research/differential.expression/2016.11/withExtra2RNAlibraries/without.mP1748.mA1567/19.01.2017.counts.txt", header = TRUE, row.names = 1)
  data = read.table("/home/botanik/Documents/phd/research/siRNA/differential.expression/counts_sRNA20-22.txt", header = TRUE, row.names = 1)
  data = read.table("/data/phdData/HyLiTE/2017.03.08/2017.08.09.HyLiTE/HyLiTE.100/HyLiTE.100.expression.txt", header = TRUE, row.names = 1)
  data = read.table("/data/phdData/HyLiTE/2017.03.08/all.HyLiTE.analysis/counts.txt", header = TRUE, row.names = 1)
  #data = read.table("~/Documents/phd/research/aglaia/counts.matrix", header = TRUE, row.names = 1)
## Combine replicate

combineRep = function(data){
    data = transform(data,  tB1830=tB1830_1+tB1830_2, tB1833=tB1833_1+tB1833_2, ## ebudensis
              fP1707=fP1707_1+fP1707_2, ##fuchsii
              mA1567=mA1567_1+mA1567_2,  mP1722= mP1722_1+mP1722_2, mP1748=mP1748_1+mP1748_2, mS1765=mS1765_1+mS1765_2, ##majalis
              tB1798=tB1798_1+tB1798_2, tB1805=tB1805_1+tB1805_2, tB1812=tB1812_1+tB1812_2, tS1920=tS1920_1+tS1920_2) ##traunteineri
    data = data[, -grep("_", colnames(data))]
    return(data)
}

species <- c()
#dataC = dataC[,-grep("_", colnames(dataC))]

# ----+ functions: makeVectors
makeSpeciesVector = function(data){
    indm = with(data, grepl("m", colnames(data))); indm[indm==1] = 2
    indt = with(data, grepl("t", colnames(data))); indt[indt==1] = 3
    inde = with(data, grepl("e", colnames(data))); inde[inde==1] = 3
    indf = with(data, grepl("f", colnames(data))); indf[indf==1] = 5
    indi = with(data, grepl("i", colnames(data))); indi[indi==1] = 6 #modifie here if you want to make ebu => traun
    ind = indm+indt+inde+indf+indi; ind = sub(2, "majalis", ind); ind = sub(3, "traunsteineri", ind); ind = sub(5, "fuchsii", ind); ind = sub(6, "incarnata", ind)
    species <<- ind
    return(ind)
}
makeGeographiesVector = function(data){
    indA = with(data, grepl("A", colnames(data))); indA[indA==1] = 2
    indP = with(data, grepl("P", colnames(data))); indP[indP==1] = 3
    indB = with(data, grepl("B", colnames(data))); indB[indB==1] = 4
    indS = with(data, grepl("S", colnames(data))); indS[indS==1] = 5
    ind = indA+indP+indB+indS; ind = sub(2, "alps", ind); ind = sub(3, "pyrenees", ind); ind = sub(4, "britain", ind); ind = sub(5, "scandinavia", ind)
    #geographies <<- ind
    return(ind)
}
# Homoeolog function helper
makeOriginsVector = function(data){
    indf = with(data, grepl("f", colnames(data))); indf[indf==1] = 1
    indi = with(data, grepl("i", colnames(data))); indi[indi==1] = 2
    ind = indf+indi; ind = sub(1, "fuchsiiSubGenome", ind); ind = sub(2, "incarnataSubGenome", ind)
    return(ind)
}
makeBatchMatrix = function(data){
    batchIndexVector = with(data, grepl("_", colnames(data)))
    batchVec = which(batchIndexVector == TRUE)
    return(matrix(batchVec, byrow = TRUE, ncol = 2))
}
#---- -

## MAPS
geoMap = c("alps"=17, "britain" = 19, "pyrenees" = 15, "scandinavia" = 18)
colorMap = c("ebudensis"="darkgreen","fuchsii"="purple","incarnata"="orange","majalis"="red","traunsteineri"="blue")
##Could add a Batch Map
##ploidies = c(rep("4",4), rep("2", 10), rep("4", 27))
##ploidyMap = c("2"="blue", "4"="red")



## HELPER FUNCTIONS

# ----+ functions: selectFac
selectSpeciesFac = function(ind){
  return(as.factor(species[ind]))
}
selectGeographiesFac = function(ind){
  return(as.factor(geographies[ind]))
}
selectOriginsFac = function(ind){
  return(as.factor(origins[ind]))
}
##selectPloidiesFac = function(ind){
##  return(as.factor(ploidies[ind]))
##}
#---- -
# ----+ function: Volcano Plot
VolcanoPlot = function(resEdgeRTable, lfc, pval){
  tab = data.frame(logFC = resEdgeRTable[, 1], negLogPval = -log10(resEdgeRTable[, 2]))
  par(mar = c(5, 4, 4, 4))
  plot(tab, pch = 16, cex = 0.6, xlab = expression(log[2]~fold~change),
       ylab = expression(-log[10]~pvalue), main = "Volcano plot")
  ## Log2 fold change and p-value cutoffs
  lfc = lfc
  pval = pval
  ## Selecting interest genes
  signGenes = (abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval))
  ## Identifying the selected genes
  points(tab[signGenes, ], pch = 16, cex = 0.8, col = "red")
  abline(h = -log10(pval), col = "green3", lty = 2)
  abline(v = c(-lfc, lfc), col = "blue", lty = 2)
  mtext(paste("pval =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1)
  mtext(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc),
        cex = 0.8, line = 0.5)
}
#---- -



## FUNCTIONS

# ----+ function:  selectSpecies
selectSpecies = function(data, spec1, spec2, spec3, spec4, spec5){
  index1 = with(data, grepl(spec1, colnames(data)))
  index2 = with(data, grepl(spec2, colnames(data)))
  index3 = with(data, grepl(spec3, colnames(data)))
  index4 = with(data, grepl(spec4, colnames(data)))
  index5 = with(data, grepl(spec5, colnames(data)))
  ind = (index1 | index2 | index3 | index4 | index5) # To remove duplicated columns when there is for instance a comparaison between "m" and "A"
  species <<- makeSpeciesVector(data)
  geographies <<- makeGeographiesVector(data)
  origins <<- makeOriginsVector(data)
  specFac <<- selectSpeciesFac(ind) # I made sure the order of the species fac & geographies fac are correct
  geoFac <<- selectGeographiesFac(ind)
  oriFac <<- selectOriginsFac(ind)
  geoNum <<- geoMap[as.character(geoFac)]
  specieColors <<- colorMap[as.character(specFac)]
  batchMatrix <<- makeBatchMatrix(data = data[ind])
  #ploidyFac <<- selectPloidiesFac(ind)
  return(data[ind])
}
#---- -
## spec1 & spec2 are the first character of the chosen species
## eg. fuchsii = "f" or "fP" or "A", one spec should be " " if only two comparaisons
## (Should be put in the variable dat)


## RUV (REMOVAL OF UNWANTED VARIANCE)

# ----+ function: makeRUVset
makeRUVset = function (dat){
  #plotRLE(as.matrix(dat),outline=FALSE,ylim=c(-4,4),col=specieColors, main = "Unwanted Variance not removed")
  uq = betweenLaneNormalization(as.matrix(dat), which = "full")
  set = newSeqExpressionSet(uq)
  plotRLE(set,outline=FALSE,ylim=c(-4,4),col=specieColors, main = "Unwanted Variance not removed")

  plotPCA(set, col = specieColors, k=2, pch = geoNum, cex = 1.2, main = "First two PCs'", labels = T)#pch = geoNum

  return(set)
}
#---- -
## Quick analysis without any normalisation
## (Should be put in variable set)


## TODO STUFF, NORMALISATION WITH REPLICATE SAMPLE
# ----+ function: makeRUVnormalisedSet
## Normalise with genes, the "empirical" ones found with findNonDEgenes
makeRUVnormalisedSet = function(set, empircalNonDEgenes, numFactors){
  normalisedSet = RUVg(as.matrix(counts(set)), empircalNonDEgenes, k=numFactors)
  plotRLE(normalisedSet$normalizedCounts,outline=F, col=specieColors, cex = 0.5, main = "Unwanted Variance removed")
  plotPCA(normalisedSet$normalizedCounts, k=2, pch=geoNum, col = specieColors, main = "First two PCs'")
  return(normalisedSet)
}
#---- -
## (Should be put in variable normSet)

#Using replicate samples to remove variance
makeRUVrepNormalisedSet = function(set, empircalNonDEgenes, numFactors, replicates){
  normalisedSet = RUVs(as.matrix(counts(set)), cIdx=empircalNonDEgenes, k=numFactors, scIdx=replicates)
  plotRLE(normalisedSet$normalizedCounts,outline=F, col=specieColors, cex = 0.5, main = "Unwanted Variance removed")
  plotPCA(normalisedSet$normalizedCounts, k=2, pch=geoNum, col = specieColors, main = "First two PCs'")
  return(normalisedSet)
}
                                        #---- -



# ----+ function: findNonDEgenes:
findNonDEgenes = function(model, set, numberOfGenes){
  design = model.matrix(model)
  y = DGEList(counts = counts(set), group = specFac)
  #y = calcNormFactors(y, method = "upperquartile")
  y = estimateGLMCommonDisp(y, design)
  y = estimateGLMTagwiseDisp(y, design)
  fit = glmFit(y, design)
  lrt = glmLRT(fit)
  top = topTags(lrt, n=nrow(set))$table
  #empirical <<- rownames(set)[which(!(rownames(set)%in%rownames(top)[1:numberOfGenes]))]
  return(rownames(set)[which(!(rownames(set)%in%rownames(top)[1:numberOfGenes]))]) # Should be put in variable empirical genes
  #return(rownames(top)[top$FDR >= numberOfGenes])
}
#---- -
## Using genes / set has to be used for normalisation.
## You can then use set3 normCounts for DE study The non DE genes are found with EdgeR
## eg. model = ~specFac+geoFac, set is the set returnd by makeRUVset
## (Should be put in variable empirical)
                                      

# name of the transcripts with p.vale smaller than... GLM method: rownames(top[top[, "PValue"] <0.01,])


# ----+ function: findEdgeRgenes
#findEdgeRgenes(~specFac+normSet$W-1, group = specFac, set = set, 0.05, c(-1,1,0,0)) model-1 is important
dat = selectSpecies(data[rowSums(cpm(data[,c(5,6,7,8,9)])>0.7) >= 2 & rowSums(cpm(data[,c(10,11,12,13,14)])>0.7) >= 2,], "f", "i", "m","t"," ")
findEdgeRgenes = function(model, group, set, pvalue, contrasts = c(0,0,0,0)){

    y = DGEList(counts = counts(set), group = specFac) # my stuff: counts(set), group = group) # Don't forget, model is taking into account the RUV normalised matrix of counts

    design = model.matrix(model)
    #colnames(design) <- levels(specFac)
    y = calcNormFactors(y, method = "upperquartile")
    y = estimateGLMCommonDisp(y, design)
    y = estimateGLMTagwiseDisp(y, design)
    plotMDS.DGEList(y, dim.plot = c(1,2), col = specieColors, cex = 1.5)
    fit <<- glmFit(y, design)
    lrt = glmLRT(fit, contrast = c(contrasts, rep(0,2)))
    fitQL <- glmQLFit(y, design, robust = T)
    lrt = glmQLFTest(fitQL, contrast = c(contrasts, rep(0,2)))
    resEdgeRtopTags <<- topTags(lrt, n=nrow(set))
    resEdgeRglm <<- lrt
    de = decideTestsDGE(lrt, adjust.method="fdr", p.value = pvalue)
    de.genes = rownames(lrt)[as.logical(de)]
    set.seed(1)
    randRes = resEdgeRtopTags[rev(rownames(resEdgeRtopTags)),]
    plotSmear(lrt, de.tags = de.genes, cex=0.5, main="MA-plot (p-value <= 0.05)")
    plot(randRes$table$logCPM, randRes$table$logFC, col=ifelse(randRes$table$FDR <= 0.05, "red", "black"), 
         cex=ifelse(randRes$table$FDR <= 0.05, 0.4, 0.3), pch = ifelse(randRes$table$FDR <= 0.05, 19, 1), xlab = "logCPM", ylab = "logFC")
    plot(resEdgeRtopTags$table$logCPM, resEdgeRtopTags$table$logFC, col=ifelse(resEdgeRtopTags$table$FDR <= 0.05, red, black), cex=ifelse(resEdgeRtopTags$table$FDR <= 0.05, 0.4, 0.3), pch = ifelse(resEdgeRtopTags$table$FDR <= 0.05, 19, 1))
    plotBCV(y)
    abline(h=c(-2,2), col="blue")
    abline(h=0)
    #VolcanoPlot(resEdgeRtopTags$table, 2, pvalue)
    #return(resEdgeR) # return a TopTags-class object
    #return(fit) # return the DGEGLM-class object
}
#striproundtwo <- rownames(fitQL$coefficients[fitQL$coefficients[,1] > -17 & fitQL$coefficients[,2] > -17,])

#LOOPS
cpmval = seq(0.5, 4, by = 0.25)
for(j in c(3,4))
{
for(i in cpmval)
#for(k in c(1,3,5,7))
{
  dat = selectSpecies(data[rowSums(cpm(data)>3) >= 4,], " ", " ", "m","t"," ")
  dat = selectSpecies(data, " ", " ", "m","t"," ")
  dat = selectSpecies(data[rowSums(cpm(data)>i) >= j,], "f", "i", "m","t"," ")
  #dat = selectSpecies(data[rowSums(cpm(data[,c(5,6,7,8,9)])>i) >= j | rowSums(cpm(data[,c(10,11,12,13,14)])>i) >= j,], "f", "i", "m","t"," ")
  set = makeRUVset(dat = dat)
  normSet = makeRUVrepNormalisedSet(set, rownames(counts(set)), 2, batchMatrix)
  y = DGEList(counts = counts(set), group = specFac) # my stuff: counts(set), group = group) # Don't forget, model is taking into account the RUV normalised matrix of counts
  design = model.matrix(~specFac+normSet$W-1)
  y = calcNormFactors(y, method = "upperquartile")
  y = estimateGLMCommonDisp(y, design)
  y = estimateGLMTagwiseDisp(y, design)
  fitQL <- glmQLFit(y, design, robust = T)
  lrt = glmQLFTest(fitQL, contrast = c(-1,1,rep(0,2)))
  resEdgeRtopTags <- topTags(lrt, n=nrow(set))
  retained = c(retained, nrow(dat))
  de = c(de, nrow(resEdgeRtopTags$table[resEdgeRtopTags$table$FDR <= 0.05,]))
  stripround = c(stripround, length(rownames(fitQL$coefficients[fitQL$coefficients[,1] < -15 | fitQL$coefficients[,2] < -15 | fitQL$coefficients[,3] < -15 | fitQL$coefficients[,4] < -15,])))
}
plot(cpmval, retained)
plot(cpmval, stripround)
plot(cpmval, de)
retained = c()
stripround = c()
de = c()
}

EdgeRresiduals = function(fit, residualType){ #residualType is "pearson","deviance" or "midp-quantile" depending on wish
  #fit is a DGEGLM-class object retured by my function findEdgeRgenes
  y = fit$counts 
  mu = fit$fitted.values
  phi = fit$dispersion 
  v = mu*(1+phi*mu)
  if(residualType == "pearson"){return((y-mu) / sqrt(v))}
  if(residualType == "deviance"){d = nbinomUnitDeviance(y,mu,phi); return(sign(y-mu) * sqrt(d))}
  if(residualType == "midp-quantile"){return(zscoreNBinom(y,mu,size=1/phi))}
}

findDEseqGenes = function(model,set, pvalue, normalised){
  if(normalised==TRUE){
    y <<- normCounts(set) # normCount should be changed to count
  }
  else{
    y <<- counts(set) # normCount should be changed to count
  }
  
  y = estimateSizeFactors(y)
  y = estimateDispersions(y)
  y = nbin
}
#---- -
## EdgeR DE study MODIFY IF NEEDED !!!!!
## Redo the DE study but without a GLM to estimate the NB dispersion parameters,
## todo directly on RUV normalised counts
## (Should be put in variable edgeRres)



# Visualising with adegenet (should be put in variable Obj)
adegenetVis = function(set, EdgeRgenes, significanceMeasure, significanceValue){ # significanceMeasure == "PValue" or "FDR", significanceValue = 0.01, 0.05
  setT=t(counts(set)[rownames(EdgeRgenes[EdgeRgenes[, significanceMeasure] <= significanceValue,]),])
  #ploidy = as.numeric(as.vector(ploidyFac))
  obj = df2genind(setT, type = "PA",sep=" ")
  other(obj)$geo = geoNum
  pop(obj) = as.vector(twoMajGroups)
  numPop = length(levels(pop(obj)))
  
  # Play with colors
  #temp <- as.integer(pop(obj))
  myCol <- transp(specieColors,.7)
  popColors = specieColors[levels(specFac)[1:length(levels(specFac))]] # To make sure the order of the colors follow the order of species apperance
  
  # Do PCA
  #pca1 <- dudi.pca(setT,scannf=FALSE)#, nf=3)
  #plot(pca1$li, col=myCol, pch=other(obj)$geo, cex=0.7)
  #grp <- find.clusters(obj, max.n.clust=9)
  
  # Do DAPC
  dapc1 <- dapc(obj, pop(obj), n.pca=3, n.da=(numPop - 1))
  if((numPop - 1) >= 3){
    scatter(dapc1, xax=1, yax=2, scree.da=T, bg="white", cell=2, cstar=0, col=popColors, solid=.6, cex=1,clab=0, leg=TRUE, txt.leg=levels(pop(obj)))
    scatter(dapc1, xax=1, yax=3, scree.da=T, bg="white", cell=2, cstar=0, col=popColors, solid=.6, cex=1,clab=0, leg=TRUE, txt.leg=levels(pop(obj)))
    scatter(dapc1, xax=2, yax=3, scree.da=T, bg="white", cell=2, cstar=0, col=popColors, solid=.6, cex=1,clab=0, leg=TRUE, txt.leg=levels(pop(obj)))
  }
  else{
    scatter(dapc1, xax=1, yax=2, scree.pca=F, posi.pca="topleft", scree.da=T, bg="white", cell=2, cstar=0, col=popColors, solid=.6, cex=1,clab=0, leg=TRUE, txt.leg=levels(pop(obj)))
  }
  contrib <<- loadingplot(dapc1$var.contr, axis=1, thres=.007, lab.jitter = 1) # Check wich contig contributes to the dapc the most
  return(obj)
}

                                        # MAIN
resave <- function(..., list = character(), file) {
   previous  <- load(file)
   var.names <- c(list, as.character(substitute(list(...)))[-1L])
   for (var in var.names) assign(var, get(var, envir = parent.frame()))
   save(list = unique(c(previous, var.names)), file = file)
}

#DO ANALYSIS FUNCTION
doAnalysis = function (Data) {
    pdf(file = paste("/home/botanik/Documents/phd/research/differential.expression/2016.11/withExtra2RNAlibraries/without.mP1748.mA1567/modelsContrasts/trimmed.counts.txtimages/", levels(specFac)[1], levels(specFac)[2], "continent" ,"pdf", sep="."))
    #par(mfrow=c(4,2))
    set = makeRUVset(Data)
    #empirical = findNonDEgenes(~specFac, set, 5000)
    normSet = makeRUVrepNormalisedSet(set, rownames(counts(set)), 2, batchMatrix) #2 cannot be changed with less than 2 replicates. eg fuchsii
    #edgeRres = findEdgeRgenes(~specFac, normSet, 0.05)
    dev.off()
    pdf(file = paste("/home/botanik/Documents/phd/research/differential.expression/2016.11/withExtra2RNAlibraries/without.mP1748.mA1567/modelsContrasts/trimmed.counts.txtimages/", levels(specFac)[1], levels(specFac)[2], "continent", "stats", "pdf", sep="."))
    findEdgeRgenes(~specFac+normSet$W, specFac, set,0.05)
    hist(resEdgeRtopTags$table$FDR, breaks = 100)
    hist(resEdgeRtopTags$table$logFC, breaks = 100)
    dev.off()
    #Obj = adegenetVis(normSet, edgeRres, "PValue", 1)
    FDRres = topTags(resEdgeRtopTags, n=nrow(dat))$table
    write.table(FDRres, paste("/home/botanik/Documents/phd/research/differential.expression/2016.11/withExtra2RNAlibraries/without.mP1748.mA1567/modelsContrasts/trimmed.counts.txt/DElists/", levels(specFac)[1], levels(specFac)[2], "continent", "txt", sep="."), sep="\t")
    write.table(FDRres[FDRres$FDR <= 0.05,], paste("/home/botanik/Documents/phd/research/differential.expression/2016.11/withExtra2RNAlibraries/without.mP1748.mA1567/modelsContrasts/trimmed.counts.txt/DElists/", levels(specFac)[1], levels(specFac)[2], "continent", "txt", sep="."), sep="\t", quote = F)
    save(resEdgeRtopTags, file = paste("/home/botanik/Documents/phd/research/differential.expression/2016.11/withExtra2RNAlibraries/without.mP1748.mA1567/modelsContrasts/trimmed.counts.txt/DEanalysisObjects/", "no.low.expression.majalis", "traunsteineri", "RData", sep="."))
    resave(resEdgeRglm, file = paste("/home/botanik/Documents/phd/research/differential.expression/2016.11/withExtra2RNAlibraries/without.mP1748.mA1567/modelsContrasts/trimmed.counts.txt/DEanalysisObjects/", "no.low.expression.majalis", "traunsteineri", "RData", sep="."))
    resave(set, file = paste("/home/botanik/Documents/phd/research/differential.expression/2016.11/withExtra2RNAlibraries/without.mP1748.mA1567/modelsContrasts/trimmed.counts.txt/DEanalysisObjects/", "no.low.expression.majalis", "traunsteineri", "RData", sep="."))
    resave(normSet, file = paste("/home/botanik/Documents/phd/research/differential.expression/2016.11/withExtra2RNAlibraries/without.mP1748.mA1567/modelsContrasts/trimmed.counts.txt/DEanalysisObjects/", "no.low.expression.majalis", "traunsteineri", "RData", sep="."))
}

for(s1 in 1:length(listSpecies)){
  for(s2 in 1:length(listSpecies)){
    if(s1 < s2 & s2 <= length(listSpecies)){
      #print(paste(listSpecies[s1],listSpecies[s2]))
      dat = selectSpecies(data, listSpecies[s1],listSpecies[s2], " ", " ", " ")
      #geoNum = geoMap[as.character(geoFac)]
      #specieColors = colorMap[as.character(specFac)] 
      #ploidyColors = ploidyMap[as.character(ploidyFac)]
      doAnalysis(Data = dat)
    }
  }
}

#Contrasts loop
contrast = list(c(-1,1,0,0,0,0), c(-1,0,1,0,0,0), c(-1,0,0,1,0,0), c(0,-1,1,0,0,0),c(0,-1,0,1,0,0),c(0,0,-1,1,0,0))
for(i in seq(1,6,by=1))
{
  lrt = glmQLFTest(fitQL, contrast = contrast[[i]])
  resEdgeRtopTags <- topTags(lrt, n=nrow(set))
  FDRres = topTags(resEdgeRtopTags, n=nrow(dat))$table
  write.table(FDRres[FDRres$FDR <= 0.05,],  paste("/home/botanik/Documents/phd/research/differential.expression/2016.11/withExtra2RNAlibraries/without.mP1748.mA1567/modelsContrasts/trimmed.counts.cpm/DElists/", contrast[i], "txt", sep="."), sep="\t", quote = F)
  save(resEdgeRtopTags, file = paste("/home/botanik/Documents/phd/research/differential.expression/2016.11/withExtra2RNAlibraries/without.mP1748.mA1567/modelsContrasts/trimmed.counts.cpm/DEanalysisObjects/", contrast[i], "RData", sep="."))
}
  
batchMaj = c("2013","2013","20140423","20140423","20140725","20140612","20140130",
          "2013","2013","20140423","20140612","20140130")

batchTrau = c("20140707","20140504","20140612","20140130","20140725","2013","20140423","20140612",
              "20140130","2013","20140612","20140612","20140130","2013","20140612","20140612","20140130")

batchFuc = c("20140423", "2013", "20140423", "20140423", "2013")
batchInc = c("2013", "20140423", "2013", "20140423", "2013")

#HyLiTE replicate matrix
batchMatrix = matrix(data=c(9,11,10,12,17,19,18,20,27,19,28,30,31,33,32,34,35,37,36,38,39,41,40,42,43,45,44,46,51,53,52,54), ncol=2, byrow=TRUE)
replaceHeader <- function(data){
    names(data) <- c(
        "mA1568","mA1573","mA1661","mA1775","mP1722_1","mP1722_2","mP1744","mS1757","mS1765_1","mS1765_2",
        "eB1830_1","eB1830_2","eB1833_1","eB1833_2",
        "tA1553","tA1641","tA1670","tB1798_1","tB1798_2","tB1805_1","tB1805_2","tB1812_1","tB1812_2","tS1901","tS1902","tS1920_1","tS1920_2",
        "fB1804","fB1855","fP1707_1","fP1707_2","fP9375",
        "iA1586","iB1176","iB1870","iS1904","iS1908")
        #"mS1765","eB1830_1","eB1830_2","eB1833_1",
        #             "eB1833_2","tA1553","tA1641","tA1670","tB1798_1",
        #             "tB1798_2","mA1573","tB1805_1","tB1805_2","tB1812_1",
        #             "tB1812_2","tS1901","tS1902","tS1920_1","tS1920_2",
        #             "mA1661","mA1775","mP1722_1","mP1722_2","mP1744",
        #             "mS1757","mS1765","fB1804","fB1855","fP1707_1",
        #             "fP1707_2","fP937","iA1586","iB1176","iB1870",
        #             "iS1904","iS1908")
}
