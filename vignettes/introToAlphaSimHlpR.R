## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load packages------------------------------------------------------------
# Make sure you have the right packages installed
library(tidyverse); library(magrittr); library(AlphaSimR); library(alphaSimHlpR)
mutate<-dplyr::mutate # fixes the fact the AlphaSimR has a mutate() function
## NOTE: if you want mutate() from AlhpaSimR you have to do AlphaSimR::mutate() now.

## ---- cols.print=13-----------------------------------------------------------
schemeParams<-tibble(stageName=c("SDN", "CET", "PYT", "AYT"),
       nEntries=c(200, 60, 20, 10), # JL used nCrosses*nProgeny for nEntries at SDN stage...
       nReps=c(1, 1, 2, 2), # Number of reps used in each stage
       nLocs=c(1, 1, 1, 2), # Number of locations used in each stage
       nChecks=c(0, 1, 1, 1), # Number of checks used in each stage
       # Checks are replicated the same as experimental entries
       # Checks are to provide connectivity across trial stages and years
       nCheckPlotsPerRep=c(0,10,10,5), # instead of using entryToCheckRatio, specify numb. check plots per rep
       phenoEval=c(F,T,T,T), # When does pheontyping occur
       # this would also facilitate adding a "conservation" stage
       # in which lines can be sent, with small cost to hang out 
       # and be potential selection candidates in the future
       newOffspringEnter=c(T,F,F,F), # How to specify where to put results of crossing
       candidateParents=c(T,F,F,F), # from which stages can we draw parents for crossing
       errVars=c(NA,146,82,40),
       varietiesExit=c(F,F,F,T)); # not totally sure what to do with this... 
       # when pipeline ends (almost always) we output varieties
schemeParams

## -----------------------------------------------------------------------------
initialPop<-initiatePop(nFounders=50,nChr=1,segSites=50,nQTL=5,nSNP=40,genVar=40,meanDD=0.8,varDD=0.05,
                        schemeParams=schemeParams)
initialPop # tidy output

## -----------------------------------------------------------------------------
initialTP<-runPhenoTrial(nReps=2,nLocs=2,nChecks=0,
                         entries = initialPop$founders[[1]],
                         SP = initialPop$SP[[1]],
                         errVars = 40,
                         trialName = "HistoricalData")
initialTP

## -----------------------------------------------------------------------------
K<-makeGRM(pop=initialPop$founders[[1]],
           SP=initialPop$SP[[1]],
           type="Add")
ini_gblup<-fitModel(phenoDF = initialTP, grm = K)

## -----------------------------------------------------------------------------
initialProgeny<-selectAndCross(candSelCrit=ini_gblup$blups[[1]],
                               candPop=initialPop$founders[[1]],
                               nParents=10, nProgeny=10, nCrosses=20,
                               SP=initialPop$SP[[1]])

## -----------------------------------------------------------------------------
ini_blup<-fitModel(phenoDF = initialTP)

## -----------------------------------------------------------------------------
initrials<-advanceVDP(stageName="CET",
                      entries=initialPop$founders[[1]],
                      schemeParams = schemeParams, 
                      candSelCrit = ini_gblup$blups[[1]])

## ----input params-------------------------------------------------------------
# DF to keep track of all germplasm
germplasm<-bind_rows(tibble(yearOfOrigin=0,Pop=list(initialPop$founders[[1]])),
                     tibble(yearOfOrigin=1,Pop=list(initialProgeny)))
# DF to keep track of field (phenotype) data 
## by adding to this each year, we accumulate training data
## and track germplasm through the pipeline
fielddata<-initialTP %>% 
      mutate(year=0)
# DF with trials to run in year 1
## basically the initial staging
## all years after will follow the breeding scheme schematic DF
## must correspond with stages in the scheme so 
## each trial can be logically advanced through the pipeline
## will overwrite this each year, so it represents the current years _fields_ 
initrials<-bind_rows(tibble(stageName="SDN", 
                            entries=list(initialProgeny)),
                     initrials)
# DF specifying pipeline
schemeParams
# sim params
SP<-initialPop$SP[[1]]
# checks
checks<-initialPop$checks[[1]]

nYears2run<-6 # years to run breedings cheme for
nParents <- 10 # Number of parents in the crossing nursery
nCrosses <- 20 # Number of crosses entering the pipeline
nProgeny <- 10 # Number of progeny per cross

## -----------------------------------------------------------------------------
testrun<-runBreedingProgram(nYears2run,schemeParams,
                            germplasm,fielddata,initrials,SP,checks,
                            nParents,nCrosses,nProgeny)

## -----------------------------------------------------------------------------
testrun[["germplasm"]][[1]] %>% 
  mutate(meanGV=map_dbl(Pop,meanG)) %>% 
  ggplot(.,aes(x=yearOfOrigin,y=meanGV)) + geom_point() + geom_line() + 
  theme_bw() + 
  geom_point(data=testrun[["released"]][[1]] %>% 
               reduce(.,rbind) %>% 
               unnest(BestClone) %>% 
               mutate(meanGV=map_dbl(entries,meanG)),
             aes(x=yearReleased,y=meanGV),
             color='red')

