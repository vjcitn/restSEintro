
ii = installed.packages()
pks = ii[,"Package"]
if (!("BiocManager" %in% pks)) install.packages("BiocManager")
if (!("fission" %in% pks)) BiocManager::install("fission", ask=FALSE, update=FALSE)

library(fission)
data(fission)
fission

fission[1:3, 1:5]

assay(fission[1:3, 1:5])

rowData(fission)

rownames(fission) = rowData(fission)$symbol
assay(fission["mcm5", 1:5])

head(colData(fission))

table(fission$strain)

boxplot(split(assay(fission["mcm5", fission$strain=="wt"]), fission$minute[fission$strain=="wt"]), main="MCM5 trajectory in Wild Type")

boxplot(split(assay(fission["mcm5", fission$strain=="mut"]), fission$minute[fission$strain=="mut"]), main="MCM5 trajectory in Mutant")

if (!("restfulSE" %in% pks)) BiocManager::install("restfulSE", ask=FALSE, update=FALSE)
library(restfulSE)

mubr = se1.3M()

mubr

assay(mubr)

apply(as.matrix(assay(mubr[,1:5])), 2, sum)

H3f3a = which(rowData(mubr)$symbol == "H3f3a")
Ptma = which(rowData(mubr)$symbol == "Ptma")
ana = function(x) as.numeric(assay(x))
plot(ana(mubr[H3f3a, 1:5000]), ana(mubr[Ptma, 1:5000]))

if (!("BiocOncoTK" %in% pks)) BiocManager::install("BiocOncoTK", ask=FALSE, update=FALSE)
if (!("bigrquery" %in% pks)) BiocManager::install("bigrquery", ask=FALSE, update=FALSE)
library(BiocOncoTK)

bq = pancan_BQ(billing="vince1-168719")
library(bigrquery)
dbListTables(bq)

BLCArna = buildPancanSE(bq, acronym="BLCA", assay="RNASeqv2")
BLCArna

BLCArnan = replaceRownames(BLCArna)
BLCArnan

names(colData(BLCArnan))

summary(BLCArnan$age_at_initial_pathologic_diagnosis)

table(BLCArnan$icd_10)

boxplot(split(log(ana(BLCArnan["FGFR3",])+1), cut(BLCArnan$age_at_initial_pathologic_diagnosis, seq(30,95,10))))
