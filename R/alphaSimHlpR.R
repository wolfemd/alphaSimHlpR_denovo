


#' Title
#'
#' @param nFounders
#' @param nChr
#' @param segSites
#' @param schemeParams
#' @param nQTL
#' @param nSNP
#' @param genVar
#' @param meanDD
#' @param varDD
#'
#' @return
#' @export
#'
#' @examples
initiatePop<-function(nFounders,nChr,segSites,schemeParams,nQTL,nSNP,genVar,meanDD,varDD){
      # This version will use runMacs
      # Future, should use runMacs2 for flexible interface
      founderPop <- runMacs(nInd=nFounders, nChr=nChr, segSites=segSites)
      # New global simulation parameters from founder haplotypes
      SP <- SimParam$new(founderPop)
      # Additive and dominance trait architecture
      SP$addTraitAD(nQtlPerChr=nQTL, var=genVar, meanDD=meanDD, varDD=varDD, useVarA=FALSE)
      # Observed SNPs per chromosome
      SP$addSnpChip(nSNP)
      founders <- newPop(founderPop, simParam=SP)
      # Choose checks from founders
      if(any(schemeParams$nChecks > 0)){
            # choosing checks at random
            # but might want to choose checks with highest EBV
            checks <- selectInd(founders, nInd=max(schemeParams$nChecks), use="rand", simParam=SP)
      } else { checks <- NULL }
      initList <- tibble(SP=list(SP), founders=list(founders), checks=list(checks))
      return(initList)
}

#' runPhenoTrial function
#'
#' @param nReps
#' @param nLocs
#' @param nChecks
#' @param entries
#' @param checks
#' @param SP
#' @param errVars
#' @param trialName
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
runPhenoTrial<-function(nReps,nLocs,nChecks,entries,checks=NULL,SP,errVars,trialName,...){
      # NOTE: according to AlphaSimR, setPheno(reps=_) results
      #### in entry means, which are extracted with pheno()
      if(!is.null(checks) & nChecks>0){
            # here, if there are checks for this stage
            # they are added to the entry vector...
            # the result being that when setPheno() is subsequently called
            # entries _and_ checks are replicated.
            entries<-c(entries, checks)
      }
      trialrecords<-setPheno(entries, varE=errVars, reps=nReps*nLocs, simParam=SP)
      trialrecords<-tibble(id=trialrecords@id,
                           mother=trialrecords@mother,
                           father=trialrecords@father,
                           trialName=trialName,
                           pheno=as.numeric(pheno(trialrecords)),  # entry means
                           errVar=errVars/(nReps*nLocs)) # entry error
      return(trialrecords)
}

#' Title
#'
#' @param M
#' @param type
#'
#' @return
#' @export
#'
#' @examples
kinship<-function(M,type){
      # Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers)
      # M is Nind x Mrow, numeric matrix, with row/columanes to indicate SNP/ind ID
      # Two types of dominance matrix, as described in
      # Vitezica et al. 2013. Genetics (Genotypic and Classical)
      # Both predict identically.
      # Difference in meaning of variance components.
      # type == "Add" should match A.mat() / Van Raden, Method 2
      M<-round(M)
      freq <- colMeans(M,na.rm=T)/2
      P <- matrix(rep(freq,nrow(M)),byrow=T,ncol=ncol(M))
      if(type=="Add"){
            Z <- M-2*P
            varD<-sum(2*freq*(1-freq))
            K <- tcrossprod(Z)/ varD
            return(K)
      }
      if(type=="Dom"){
            W<-M;
            W[which(W==1)]<-2*P[which(W==1)];
            W[which(W==2)]<-(4*P[which(W==2)]-2);
            W <- W-2*(P^2)
            varD<-sum((2*freq*(1-freq))^2)
            D <- tcrossprod(W) / varD
            return(D)
      }
}

#' Title
#'
#' @param pop
#' @param SP
#' @param type
#'
#' @return
#' @export
#'
#' @examples
makeGRM<-function(pop, SP, type){
      grm<-kinship(M=pullSnpGeno(pop=pop, simParam=SP),
                   type=type)
      return(grm)
}


#' Title
#'
#' @param phenoDF
#' @param grm
#'
#' @return
#' @export
#'
#' @examples
fitModel<-function(phenoDF, grm=NULL){
      require(sommer)
      if(!is.null(grm)){
            # Enable prediction of genos that don't have phenos
            phenoDF$id<-factor(phenoDF$id, levels=rownames(grm)) }
      randformula<-if_else(!is.null(grm),"~vs(id, Gu=grm)", "~vs(id)") %>% as.formula
      phenoDF<-phenoDF %>% dplyr::mutate(WTs=1/errVar) # Make weights
      fm <- mmer(pheno ~1,
                 random=randformula,
                 rcov= ~units,
                 weights=WTs,
                 data=phenoDF,verbose = F)
      varcomps<-summary(fm)$varcomp
      blups<-tibble(id=as.character(names(fm$U[[1]]$pheno)),
                    BLUP=as.numeric(fm$U[[1]]$pheno))
      out<-tibble(varcomps=list(varcomps),
                  blups=list(blups))
      return(out)
}

#' Title
#'
#' @param selcrits
#' @param nToSel
#' @param col
#'
#' @return
#' @export
#'
#' @examples
truncSel<-function(selcrits,nToSel,col="BLUP"){
      # selcrits[order(selcrits[[col]],decreasing = T),]$id[1:nToAdvance] # base version
      clonesToSelect<-selcrits %>%
            arrange(desc(!!sym(col))) %>%
            slice(1:nToSel) %$%
            id # vector of GIDs as output
      return(clonesToSelect)
}

#' Title
#'
#' @param candSelCrit
#' @param col
#' @param candPop
#' @param nParents
#' @param nProgeny
#' @param nCrosses
#' @param SP
#'
#' @return
#' @export
#'
#' @examples
selectAndCross<-function(candSelCrit,col="BLUP",candPop,nParents,nProgeny,nCrosses,SP){
      parents<-truncSel(selcrits = candSelCrit,
                        nToSel = nParents,col=col)
      offspring<-randCross(candPop[parents],nCrosses=nCrosses, nProgeny=nProgeny, ignoreGender=T, simParam=SP)
      return(offspring)
}

#' Title
#'
#' @param stageName
#' @param entries
#' @param schemeParams
#' @param candSelCrit
#'
#' @return
#' @export
#'
#' @examples
advanceVDP<-function(stageName,entries,schemeParams,candSelCrit){
      lastStage<-schemeParams$stageName[nrow(schemeParams)]
      if(stageName==lastStage){
            nextStage<-"released";
            nToAdvance<-1 # for now, just pick top clone at end of each pipeline run
      } else {
            nextStage<-schemeParams$stageName[which(schemeParams$stageName==stageName)+1]
            nToAdvance<-schemeParams$nEntries[schemeParams$stageName==nextStage]
      }
      clones2advance<-truncSel(selcrits = candSelCrit,
                               nToSel = nToAdvance)
      nextTrial<-tibble(stageName=nextStage,
                        entries=list(entries[clones2advance]))
      return(nextTrial)
}


#' Title
#'
#' @param nYears2run
#' @param schemeParams
#' @param germplasm
#' @param fielddata
#' @param initrials
#' @param SP
#' @param checks
#' @param nParents
#' @param nCrosses
#' @param nProgeny
#'
#' @return
#' @export
#'
#' @examples
runBreedingProgram<-function(nYears2run,schemeParams,
                             germplasm,fielddata,initrials,SP,checks,
                             nParents,nCrosses,nProgeny){
      trials<-initrials
      released<-list()
      for(i in 1:nYears2run){
            trials<-left_join(trials,schemeParams)
            # Phenotype field trials for this year
            ## Updates training data _before_ selection and advancement
            if(any(trials$phenoEval==TRUE)){
                  print(paste0("Phenotype field trials - Year ",i))
                  fieldDataCurrentYear<-trials %>%
                        # only trials to be pheno eval'd
                        filter(phenoEval==TRUE) %>%
                        # label each pheno trial with stageName
                        dplyr::mutate(trialName=stageName) %>%
                        # use runPhenoTrial() function on each trial
                        pmap(.,runPhenoTrial,SP=SP) %$%
                        # combine trial data.frames
                        reduce(.,rbind) %>%
                        # add the year
                        mutate(year=i)
                  # add current year's data to "historical" record (==fielddata data.frame)
                  fielddata<-rbind(fielddata,fieldDataCurrentYear)
            }
            # Population improvement
            if(any(trials$candidateParents==TRUE)){
                  print(paste0("Selection and Crossing - Year ",i))

                  # pool of selection candidates
                  candIDs<-trials %>%
                        filter(candidateParents==TRUE) %$%
                        entries %>%
                        mergePops(.) %$% # combine the pop objects for any trial containing sel candidates
                        unique(.@id) # grab the unique ids
                  # training samples
                  tpIDs<-unique(fielddata$id) # anyone phenotyped is in the TP
                  # Make kinship matrix
                  K<-makeGRM(pop=mergePops(germplasm$Pop)[union(tpIDs,candIDs)],
                             SP=SP,
                             type="Add")
                  # Get GEBVs
                  gblup<-fitModel(phenoDF = fielddata %>% filter(id %in% tpIDs),
                                  grm = K)
                  # Select and cross
                  newOffspring<-selectAndCross(candSelCrit=gblup$blups[[1]] %>% filter(id %in% candIDs),
                                               candPop=mergePops(germplasm$Pop),
                                               nParents=nParents, nProgeny=nProgeny, nCrosses=nCrosses,
                                               SP=SP)
                  # Add new offspring to the germplasm record
                  germplasm<-bind_rows(germplasm,
                                       tibble(yearOfOrigin=i+1,Pop=list(newOffspring)))
            }
            # Product advancement
            print(paste0("Advance the VDP - Year ",i))
            ## Get BLUPs using all available field data
            ## as selection criteria
            blup<-fitModel(phenoDF = fielddata)
            newtrials<-trials %>%
                  ## if a trial has newOffspring (seedlings) and is not phenotyped
                  ## Use GEBVs else use BLUPs
                  mutate(candSelCrit=ifelse(phenoEval==FALSE & newOffspringEnter==TRUE,
                                            list(gblup$blups[[1]]),
                                            list(blup$blups[[1]])),
                         candSelCrit=map2(entries,candSelCrit,~filter(.y,id %in% .x@id))) %>%
                  dplyr::select(stageName,entries,candSelCrit) %>%
                  pmap(.,advanceVDP,schemeParams=schemeParams) %>%
                  reduce(.,rbind)
            ## separate off the "released" (best) clone at end of pipeline
            if(any(newtrials$stageName=="released")){
                  products<-newtrials %>% filter(stageName=="released")
                  newtrials<-newtrials %>% filter(stageName!="released")
                  released[[i]]<-tibble(yearReleased=i,BestClone=list(products))
            }
            ## "Plant" trials for next year
            ## starting with the newOffspring seedling nursery
            newtrials<-bind_rows(tibble(stageName="SDN",entries=list(newOffspring)),
                                 newtrials)
            trials<-newtrials
      }
      out<-list(trials=list(trials), # existing trials at end of simulation
                fielddata=list(fielddata), # accumulated field data records
                germplasm=list(germplasm), # accumulated germplasm [pop-class objects]
                released=list(released)) # accumulated best clones from end of VDP
      return(out)
}