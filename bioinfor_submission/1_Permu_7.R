##------------------------------
## For my statistical genetics class, I wanted to see whether using a linear
## mixed effects model would be an improvement over a traditional ordinary least
## squares approach. The comparison between OLS and LMER had fewer significant
## found with LMER than with OLS. But still showed a comparable number of SNPs
## 23 for LMER compared to 30 for OLS, and 21 were the same SNPs in both tests.
##------------------------------

## Reads arguments from the command line
## this will also pass through arguments defined in the job submission file
## error logs are also created through the submission file
args=(commandArgs(TRUE))
for(i in 1:length(args)){
 eval(parse(text=args[[i]]))
}
##

## Defaults -- jobID, genotypes, and trait are defined from the submission file
if(!exists("jobID")){ jobID=104 } # for save file and defining chunk
  # this is defined as 1 through 125 in the submission file
if(!exists("genotypes")){ genotypes='calls' } # genotypes used
if(!exists("trait")){ trait='bhmd' } # trait - bone-heel mineral density
if(!exists("chunkSize")){ chunkSize=5000 } # number of SNPs in each chunk
if(!exists("bufferSize")){ bufferSize=1000 } # buffer to add to create chunk


## Load libraries
library(BGData)
library(data.table)

## Output Folder
dir.create("trait2")
setwd("trait2")
workingpath <- getwd()

## Load phenotypes
Y <- fread(paste0(
  "/mnt/research/quantgen/projects/UKB/PIPELINE500/output/adjusted_phenotypes/",
  trait, "/whites_unrelated_0.05/", genotypes, "/trait.tsv"),
  colClasses = list(character = c("FID", "IID")), data.table = FALSE,
  showProgress = FALSE)
rownames(Y) <- Y$IID

## Load genotypes
load.BGData(paste0(
  "/mnt/research/quantgen/projects/UKB/PIPELINE500/output/BGData/whites/",
  genotypes, "/BGData.RData"))
X <- DATA@geno
rm(DATA)

## Exclude missing values
y <- Y[, paste0("adjusted_", trait)] #bone-heel mineral density data
names(y) <- Y$IID
Y <- Y[!is.na(y),] # remove NA values in y from Y matrix
y <- y[!is.na(y)] # same for the y array


## Extract chunk with flanking regions (buffer)
nSNPs <- ncol(X)
iniSNP <- (jobID - 1) * chunkSize + 1 # define initial SNP in core (+1 so not 0)
endSNP <- min(nSNPs, iniSNP+chunkSize-1) # define final SNP in core
core <- iniSNP:endSNP # define core (chunk w/o flanking regions)
chunk <- seq(max(1, iniSNP-bufferSize), min(nSNPs, endSNP+bufferSize))
  # add flanking region to create chunk
isCore <- chunk %in% core # bolean to check chunk is in core
W <- X[names(y), chunk] # from the genotype matrix, extract indiduals in y
  # and SNPs in chunk

## GWAS for linear mixed-effects model with center as a random variable
GWAS0=GWAS(y~sex+age+(1|center),method='lmer',
  data=BGData(geno=W,pheno=data.frame(y=y, sex=Y$sex, age=Y$age, center=Y$center)))

Pvalue0=GWAS0[isCore,4] # save p-values if in the core SNPs


save(Pvalue0,file=paste0('chunk_',sprintf("%03d",jobID),".RData"))

quit(save='no')
