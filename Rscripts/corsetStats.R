totalRecall <- read.table("~/Documents/phd/scripts/totalRecall1.txt", sep = ";", header = T)
head(dat$cluster_name)
totalRecallclust <- totalRecall[as.vector(totalRecall$sequence_name) != as.vector(totalRecall$cluster_name),]
uniqueCLusterTable <- as.data.frame(table(totalRecallclust$cluster_name))

# Histogram of seuqence length such that the corresponding cluster contains only one sequence.
hist(totalRecallclust[uniqueCLusterTable[uniqueCLusterTable$Freq == 5,]$Var1,]$length, breaks = 100)

hist(totalRecallclust[uniqueCLusterTable$Var1,]$length, breaks = 100, col = "red")
hist(totalRecallclust[uniqueCLusterTable[uniqueCLusterTable$Freq == 0,]$Var1,]$length, breaks = 100, col = "blue", add = T)
hist(totalRecallclust[uniqueCLusterTable[uniqueCLusterTable$Freq == 1,]$Var1,]$length, breaks = 100, col = "blue", add = T)
hist(totalRecallclust[uniqueCLusterTable[uniqueCLusterTable$Freq == 2,]$Var1,]$length, breaks = 100, col = "yellow", add = T)

hist(totalRecallclust[uniqueCLusterTable[uniqueCLusterTable$Freq == 2,]$Var1,]$length, breaks = 100, col = "yellow")
hist(totalRecallclust[uniqueCLusterTable[uniqueCLusterTable$Freq == 3,]$Var1,]$length, breaks = 100, col = "orange", add = T)
hist(totalRecallclust[uniqueCLusterTable[uniqueCLusterTable$Freq == 4,]$Var1,]$length, breaks = 100, col = "magenta", add = T)
hist(totalRecallclust[uniqueCLusterTable[uniqueCLusterTable$Freq == 5,]$Var1,]$length, breaks = 100, col = "green", add = T)

hist(totalRecallclust[uniqueCLusterTable[uniqueCLusterTable$Freq >= 0,]$Var1,]$length, breaks = 100, col = "red")
hist(totalRecallclust[uniqueCLusterTable[uniqueCLusterTable$Freq >= 1,]$Var1,]$length, breaks = 100, col = "blue", add = T)
hist(totalRecallclust[uniqueCLusterTable[uniqueCLusterTable$Freq >= 2,]$Var1,]$length, breaks = 100, col = "yellow", add = T)

# Check weather the annoptations are the same for all sequences of cluster?
