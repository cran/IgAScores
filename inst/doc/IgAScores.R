## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----dummy_data---------------------------------------------------------------
#load in IgAScores
library(IgAScores)

#dataframes with counts for the bacterial taxa in the IgA+ and IgA- fractions and presort sample, as would be produced by 16S rRNA appraoches such as DADA2
igapos <- data.frame(Sample1=c(100,0,1,2,10),Sample2=c(110,0,11,42,50),Sample3=c(140,60,10,3,0))
iganeg <- data.frame(Sample1=c(200,0,40,20,4),Sample2=c(10,30,110,2,5),Sample3=c(30,20,0,123,20))
presort <- data.frame(Sample1=c(150,10,50,30,5),Sample2=c(100,30,115,20,10),Sample3=c(30,20,10,100,25))

taxnames <- c("Taxon1","Taxon2","Taxon3","Taxon4","Taxon5")
rownames(igapos) <- taxnames
rownames(iganeg) <- taxnames
rownames(presort) <- taxnames

#convert the counts to relative abundances using the included helper function
igapos <- relabund(igapos)
iganeg <- relabund(iganeg)
presort <- relabund(presort)

#iga+ and iga- fraction sizes per sample (fraction, if a percentage divide by 100)
possize <- c(Sample1=0.04,Sample2=0.05,Sample3=0.03)
negsize <- c(Sample1=0.54,Sample2=0.47,Sample3=0.33)

#set a pseudo count for handling zero values in some scoring methods
#this should be of a similar value of the minimum non-zero observed value (e.g. if minum values is 0.007 use 0.001)
pseudo <- 0.001


## ----reqs, echo=F-------------------------------------------------------------
pr <- c("Probability Ratio","probratio","IgA+ abundance, IgA- abundance, IgA+ size, IgA- size, pseudo")
pp <- c("IgA+ Probability","prob","IgA+ abundance, Presort abundance, IgA+ size")
pa <- c("Palm index","palm","IgA+ abundance, IgA- abundance, pseudo")
ka <- c("Kau index","kau","IgA+ abundance, IgA- abundance, pseudo")
comb <- rbind(pr,pp,pa,ka)
colnames(comb) <- c("Score","Method name","Inputs required")

knitr::kable(comb, row.names = F)

## ----probrat------------------------------------------------------------------
#default method is probratio
prscores <- igascores(posabunds = igapos, negabunds = iganeg, 
                      possizes = possize, negsizes = negsize, 
                      pseudo = pseudo)

print(prscores)


## ----pp-----------------------------------------------------------------------
ppscores <- igascores(posabunds = igapos, possizes = possize, presortabunds = presort, method="prob")

print(ppscores)

## ----kp-----------------------------------------------------------------------
kscores <- igascores(posabunds = igapos, negabunds = iganeg, pseudo=pseudo, method="kau")
print(kscores)

pscores <-  igascores(posabunds = igapos, negabunds = iganeg, pseudo=pseudo, method="palm")
print(pscores)

## ----singlefuns---------------------------------------------------------------
igaprobabilityratio(posabund = igapos[1,1], negabund = iganeg[1,1], possize = possize[1], negsize = negsize[1], pseudo = pseudo)
igaprobability(withinabund = igapos[1,1], presortabund = presort[1,1], gatesize = possize[1])
kauindex(posabund = igapos[1,1], negabund = iganeg[1,1], pseudo = pseudo)
palmindex(posabund = igapos[1,1], negabund = iganeg[1,1], pseudo = pseudo)


## ----sim----------------------------------------------------------------------
#run the simulation with defaults
simdata <- simulateigaseq()

summary(simdata)

