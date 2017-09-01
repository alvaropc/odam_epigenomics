#Epigenomic Practise
#Álvaro Ponce Cabrera

library(minfi)
setwd("~/Agomez_practica/")

#Loading Data

idat.folder <- "~/Agomez_practica/idats" 
targets <- read.450k.sheet(base=idat.folder)
rgset <- read.450k.exp(targets = targets)

 
#The class of RGSet is a RGChannelSet object. 
#This is the initial object of a minfi analysis that contains the raw intensities in the green and red channels. 
#Note that this object contains the intensities of the internal control probes as well. 
phenoData <- pData(rgset)
names(phenoData)
phenoData$Status
#The RGChannelSet stores also a manifest object that contains the probe design information of the array:
manifest <- getManifest(rgset)

#Probes information
getProbeInfo(manifest)

#A MethylSet objects contains only the methylated and unmethylated signals
MSet <- preprocessRaw(rgset) 


#getMeth and getUnmeth can be used to get the methylated and unmethylated intensities matrices:
head(getUnmeth(MSet))
head(getMeth(MSet))


#A RatioSet object is a class designed to store Beta values and/or M values instead of the methylated and unmethylated signals. 
#An optional copy number matrix, CN, the sum of the methylated and unmethylated signals, can be also stored. 
#Mapping a MethylSet to a RatioSet may be irreversible, i.e. one cannot be guranteed to retrieve the methylated and unmethylated signals from a RatioSet.
#A RatioSet can be created with the function ratioConvert:
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
#The functions getBeta, getM and getCN return respectively the Beta value matrix, M value matrix and the Copy Number matrix.
beta <- getBeta(RSet)

#GenomicRatioSet

#The function mapToGenome applied to a RatioSet object will add genomic coordinates to each probe together with some additional annotation information.
#The output object is a GenomicRatioSet (class holding M or/and Beta values together with associated genomic coordinates). 
#It is possible to merge the manifest object with the genomic locations by setting the option mergeManifest to TRUE.
GRset <- mapToGenome(RSet)
beta <- getBeta(GRset)
sampleNames <- sampleNames(GRset)
probeNames <- featureNames(GRset)
pheno <- pData(GRset)

#Full annotation
annotation <- getAnnotation(GRset)
names(annotation)

#Quality control
#minfi provides a simple quality control plot that uses the log median intensity in both the methylated (M) and unmethylated (U) channels.
#When plotting these two medians against each other, it has been observed that good samples cluster together, 
#while failed samples tend to separate and have lower median intensities. 
#In order to obtain the methylated and unmethylated signals, we need to convert the RGChannelSet to an object containing the methylated and unmethylated signals using the function preprocessRaw. It takes as input a RGChannelSet and converts the red and green intensities to methylated and unmethylated signals according to the special 450K probe design, and returns the converted signals in a new object of class MethylSet. It does not perform any normalization.
#The accessors getMeth and getUnmeth can be used to get the methylated and unmethylated intensities matrices:

head(getMeth(MSet)
head(getUnmeth(MSet))
#The functions getQC and plotQC are designed to extract and plot the quality control information from the MethylSet:
qc <- getQC(MSet)
plotQC(qc)
     
phenoData <- pData(rgset)
     
     
#To further explore the quality of the samples, it is useful to look at the Beta value densities of the samples, with the option to color the densities by group:
densityPlot(MSet, sampGroups = phenoData$Status)
     
#The 450k array contains several internal control probes that can be used to assess the quality control of different sample 
#preparation steps (bisulfite conversion, hybridization, etc.). 
#The values of these control probes are stored in the initial RGChannelSet 
#and can be plotted by using the function controlStripPlot and by specifying the control probe type:
     
controlStripPlot(rgset, controls="BISULFITE CONVERSION II")
     
     
#########
     #1
#########
     
     
#Sex prediction
#By looking at the median total intensity of the X chromosome-mapped probes, denoted med(X)med(X), 
#and the median total intensity of the Y-chromosome-mapped probes, denoted med(Y)med(Y),
#one can observe two different clusters of points corresponding to which gender the samples belong to.
#getSex needs to be a GenomicMethylSet or a GenomicRatioSet.
    
predictedSex <- getSex(GRset, cutoff = -2)$predictedSex
head(predictedSex)

#To choose the cutoff to separate the two gender clusters, 
#one can plot med(Y)med(Y) against med(Y)med(Y) with the function plotSex:
plotSex(getSex(GRset, cutoff = -2))

     
     
#########
     #2
#########

#Detection pvals

detP = detectionP(rgset)


#Checking if there are some sample with a pval below than 0.05
keep = colMeans(detP) < 0.05
keep #All the cases are TRUE, so we do not need to remove any saple  

###########
     #3
###########

gRatioSet.quantile <- preprocessQuantile(rgset)##SQN

#Anotation
ann450k<-getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(ann450k)
keep <- !(featureNames(gRatioSet.quantile) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
table(keep)

gRatioSet.quantile <- gRatioSet.quantile[keep,]

# Remove probes with SNPs at CpG or SBE site:
gRatioSet.quantile<- dropLociWithSnps(gRatioSet.quantile)

betas<-getBeta(gRatioSet.quantile)

colnames(betas)<-targets$Status
head(betas)


#DMP finding

#dmpFinder Function needs only 2 types or names in this case, so we change targets$Status
targets$Status<- c("Adult","Adult","Adult","Fetal","Fetal","Fetal")
gc()
dmp <- dmpFinder(betas, pheno = targets$Status  , type = "categorical")
gc()     
#We need to filter dmp by betas>0.2 and qval <0.01
#First we create a matrix with the rowmeans (beta values mean) of adults and fetals samples
betas_adult <- rowMeans(betas[,grep("^Adu", colnames(betas))], na.rm = T)
betas_fetal <- rowMeans(betas[,grep("^Fet", colnames(betas))], na.rm = T)
m_betas <- cbind(betas_adult, betas_fetal)
colnames(m_betas) <- c("Adults","Fetals")
head(m_betas)

#Check if rownames of fmp and rownames of m_betas are in the same order
identical(rownames(dmp),rownames(m_betas))
#There aren't, so we need to make them to be in the same order
head(m_betas)
head(dmp)
Match<-match(rownames(dmp),rownames(m_betas))
m_betas<-m_betas[Match,]
identical(rownames(dmp),rownames(m_betas))

#Now they are, we can start the filter
     
res<-dmp[abs(m_betas[,"Adults"]-m_betas[,"Fetals"])>.2,]

res<-subset(res,res$qval<.01)
     
dim(res)

#########
     #4
########
ann450k<-as.data.frame(ann450k)
     
dades <- merge(betas,ann450k[,c("UCSC_RefGene_Name","Relation_to_Island","UCSC_RefGene_Group")], by="row.names")
rownames(dades)<-dades$Row.names
dades<-dades[,-1]
dades<-dades[match(rownames(res),rownames(dades)),]
identical(rownames(res),rownames(dades))

#To look for promotors into dades data, we check if we have res rownames and dades rownames
#in the same order


#Now we look for promoters
proms<-dades[grep("TSS1500|TSS200|5'UTR|1stExon",dades$UCSC_RefGene_Group),]
dim(proms)
gc()
#########
      #5
#########

#Look for "islands" inside proms
proms.island <- proms[grep("Island", proms$Relation_to_Island),]
dim(proms.island)
head(proms.island)
     
################
     #6
###############

#We need to have betas data inside data adult and fetal resume in only one column
dades_adult <- rowMeans(dades[,grep("^Adu", colnames(betas))], na.rm = T)
dades_fetal <- rowMeans(dades[,grep("^Fet", colnames(betas))], na.rm = T)
dades_mean<-cbind(dades_adult,dades_fetal)
identical(rownames(dades),rownames(dades_mean))
dades1<-cbind(dades_mean,dades[,-c(1:6)])
head(dades1)

HyperFetal <- proms.island[proms.island$dades_fetal > 0.66,]
dim(HyperFetal)
HypoAdult <- HyperFetal[HyperFetal$dades_adult < 0.33,]
dim(HypoAdult)
     
     
#############
     #7
############
gc()
DiffExpr <- read.csv("diff.expressed.AdultvsFetal.csv")
names(DiffExpr)
DiffExpr <- DiffExpr[DiffExpr$logFC > 2,]
head(DiffExpr)
     
#Obtaining Gene Symbols into DiffExpr and comparing againts Gene Symbol inside proms.island
head(proms.island$UCSC_RefGene_Name)
DiffExpr<-DiffExpr[!DiffExpr$Gene.symbol=="",]
DiffExpr_GeneSymbol<-unlist(strsplit(as.character(DiffExpr$Gene.symbol), "///"))
     
head(DiffExpr_GeneSymbol)
proms.island_GeneSymbol<- unlist(strsplit(as.character(proms.island$UCSC_RefGene_Name), ";"))
                                 
head(proms.island_GeneSymbol)
length(proms.island_GeneSymbol)
length(DiffExpr_GeneSymbol)

DiffExprsMatch <- match(DiffExpr_GeneSymbol, proms.island_GeneSymbol)
head(DiffExprsMatch)

DiffExprsMatch<- DiffExprsMatch[!is.na(DiffExprsMatch)]
head(DiffExprsMatch)
length(DiffExprsMatch)

DiffExprr_genes <- proms.island[DiffExprsMatch,]
dim(DiffExprr_genes)
gc()     
############
     #8
###########
     
#GoStatsAnalysis
     library(org.Hs.eg.db)
     library(GOstats)
     library(GO.db)
     library(annotate)
     
# Genes in promoter and island differentially methylated CpGs:
head(proms.island)
head(proms.island_GeneSymbol)

geneIDs <- unlist(mget(proms.island_GeneSymbol, ifnotfound=NA, revmap(org.Hs.egSYMBOL)))
   
universe <- Lkeys(org.Hs.egGO)
     
BP <- new("GOHyperGParams", geneIds=geneIDs, universeGeneIds = universe, annotation="org.Hs.eg.db", ontology="BP", pvalueCutoff = 0.01, conditional = FALSE, testDirection = "over")

gc()
     
hypBP <- hyperGTest(BP)

# Get the p-values of the test
gGhyp.pv <- pvalues(hypBP)
gGhyp.odds<-oddsRatios(hypBP)
gGhyp.counts<-geneCounts(hypBP)

sigGO.ID <- names(gGhyp.pv[gGhyp.pv < 0.01])

#Test the number of counts

gGhyp.counts<-as.data.frame(gGhyp.counts)
gGhyp.counts$GOterms<-rownames(gGhyp.counts)
gGhyp.counts<-gGhyp.counts[rownames(gGhyp.counts) %in% sigGO.ID,]

#Here only show the significant GO terms of BP (Molecular Function)
#For other categories, just follow the same procedure.
sigGO.Term <- getGOTerm(sigGO.ID)[["BP"]]

results_GO<-cbind(as.data.frame(gGhyp.pv[gGhyp.pv < 0.01]), as.data.frame(sigGO.Term), gGhyp.counts)

write.csv(results_GO, "GO_enrichment_BP_FetalvsAdult.csv")


#########
  #9
########
head(proms.island)
heat <- as.matrix(proms.island[,c(1:6)])
head(heat) 
library(gplots)


adult <- grep("^Adult",colnames(proms.island),perl=T)
fetal<- grep("^Fetal",colnames(proms.island),perl=T)

spcol <- c(rep("blue1",length(adult)),rep("red",length(fetal)))

pdf("heatmap.diffMeth.FetalvsAdult.pdf")

heatmap.2(heat,
             main="Diffmeth CpG's ranked",
             labRow=NA,
             trace="none",
             na.rm=T,
             col=greenred,
             ColSideColors=spcol,
             distfun=function(x) dist(x,method="manhattan"),
             dendrogram = "column")
dev.off()
