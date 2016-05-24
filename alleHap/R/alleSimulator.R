#' @title Simulation of genetic data (alleles) and non-genetic data (family identifiers)
#' @description Data simulation can be performed taking into account many different factors such as number of families to generate, number of markers (allele pairs), number of different alleles per marker, type of alleles (numeric or character), number of different haplotypes in the population, probability of parent/offspring missing genotypes, proportion of missing genotypes per individual, probability of being affected by disease and recombination rate.
#' @param nFams Number of families to generate (integer: 1..1000+)
#' @param nChildren Number of children of each family (integer: 1..7 or NULL)
#' @param nMarkers Number of markers or allele pairs to generate (integer: 1..1000+)
#' @param numAllperMrk Number of different alleles per marker (vector or NULL)
#' @param chrAlleles Should alleles be expressed as characters A,C,G,T ? (boolean: FALSE, TRUE)
#' @param nHaplos Number of different haplotypes in the population (numeric)
#' @param missParProb Probability of parents' missing genotype (numeric: 0..1)
#' @param missOffProb Probability of offspring' missing genotype (numeric: 0..1)
#' @param ungenotPars Proportion of ungenotyped parents (numeric: 0..1)
#' @param ungenotOffs Proportion of ungenotyped offspring (numeric: 0..1)
#' @param phenProb Phenotype probability, e.g. being affected by disease (numeric: 0..1)
#' @param recombRate Recombination rate (numeric: 0..1)
#' @param invisibleOutput Data are not shown by default.
#' @param dataSummary A summary of the data is shown by default.
#' @return Families' genotypes and haplotypes.
#' @import stats
#' @export alleSimulator
#' @references Medina-Rodriguez, N. Santana A. et al. (2014) alleHap: an efficient algorithm to reconstruct zero-recombinant haplotypes from parent-offspring pedigrees. BMC Bioinformatics, 15, A6 (S-3).
#' @examples
#' 
#' ## Generation of 10 simulated families with 2 children per family and 20 markers
#' simulatedFams <- alleSimulator(10,2,20)  # List with simulated alleles and haplotypes
#' simulatedFams[[1]]                       # Alleles (genotypes) of the simulated families
#' simulatedFams[[2]]                       # Haplotypes of the simulated families
#' 
alleSimulator=function(nFams=2,nChildren=NULL,nMarkers=3,numAllperMrk=NULL,chrAlleles=TRUE,nHaplos=1200,
                       missParProb=0,missOffProb=0,ungenotPars=0,ungenotOffs=0,phenProb=0.2,recombRate=0,
                       invisibleOutput=TRUE,dataSummary=TRUE){  
  
  ############################################ I. INTERNAL FUNCTIONS #########################################
{  
  labelMrk=function(){                                          # Creation of the 'A','C','G','T' character labels
    labProb=sample(1:3,1)                                       # Probability label (among 1 and 3)              
    if (labProb==1) c('C','T')                                  # 'C' and 'T' labels 
    else if (labProb==2) c('G','A')                             # 'G' and 'A' labels 
    else sample(c('A','C','G','T'),2)                           # Two letters among 'A','C','G','T' labels
  }
  simHapSelection=function(n,numAllperMrk){                     # Selection of n different haplotypes among the total number of possible haplotypes
    hapLength=length(numAllperMrk)                              # Length of the haplotype 
    simHaps=matrix(0,nrow=n,ncol=hapLength)                     # simHaps matrix initialization
    for (k in 1:hapLength)                                      # Scanning of the columns of the simHaps matrix
      simHaps[,k]=sample(0:(numAllperMrk[k]-1),n,replace=TRUE)  # Sampling of each column of the matrix simHaps
    return(simHaps)                                             # Sampled simHaps matrix
  }
  simOffspring=function(nOffs,father,mother,recombRate){        # Generation of n children by selecting randomly an haplotype from each parent
    simChild=function(father,mother){                           # Randomly haplotype asignation to a child from parents
      rcb=rbinom(2,1,recombRate)                                # Randomly decision of recombination in parents
      if (rcb[1]) father=simRecombHap(father)                   # Simulation of the recombination of father's haplotypes
      if (rcb[2]) mother=simRecombHap(mother)                   # Simulation of the recombination of mother's haplotypes
      pos=sample(1:2,2,replace=TRUE)                            # Randomly decision of which haplotype from each parent is inherited by the child
      offspring=rbind(father[pos[1],],mother[pos[2],])          # Haplotypes in offspring
      haps=apply(offspring,1,paste,collapse="-")                # Convertion of offspring's haplotypes in a character vector
      markers=as.vector(apply(offspring,2,sort))                # Sorting of the alleles in alphanumeric order
      child=c(markers,sum(rcb),haps)                            # Vector with alleles, haplotypes and number of recombinant haplotypes
      return(child)
    }
    offspring=t(replicate(nOffs,simChild(father,mother)))       # Concatenation of the previous generated children
    return(offspring)
  }
  simOneFamily=function(famID,nCh,popHaplos,phenProb,recombRate){  # Simulation of one family from a population containing the haplotypes 'popHaplos'     
    if (is.null(nCh)|nCh==0|nCh>4) {
      if (is.null(nCh)|nCh==0) {
        dnof=c(0.1,0.45,0.25,0.15,0.044,0.005,0.001)               # Distribution of the number of offspring per family    
        nOffs=sample(1:length(dnof),1,prob=dnof)                   # Generation of the number of offspring
      } else {
        dnof=c(0.1,0.5,0.20,0.15,rep(0.05/(nCh-4),nCh-4))          # Distribution of the number of offspring per family    
        nOffs=4
      }
    } else nOffs=nCh 
    nHaplos=nrow(popHaplos)                                        # Number of haplotypes in population
    nMkrs=ncol(popHaplos)                                          # Number of markers in population
    if (nOffs==1) conf=sample(1:7,1)                               # Checks if the number of children is 1 and samples the configuration number
    else if (nOffs==2) conf=sample(3:7,1)                          # Checks if the number of children is 2 and samples the configuration number
    else if (nOffs==3) conf=sample(5:7,1)                          # Checks if the number of children is 3 and samples the configuration number
    else if (nOffs==4) conf=sample(6:7,1)                          # Checks if the number of children is 4 and samples the configuration number
    nucConf=c(1,1,2,2,3,4,4)                                       # Number of unique children per configuration
    nuhConf=c(1,2,2,3,2,3,4)                                       # Number of unique haplotypes per configuration
    parsConfs=list(pConf1=rbind(c(1,1),c(1,1)),
                   pConf2=rbind(c(1,1),c(2,2)),
                   pConf3=rbind(c(1,1),c(1,2)),
                   pConf4=rbind(c(1,1),c(2,3)),
                   pConf5=rbind(c(1,2),c(1,2)),
                   pConf6=rbind(c(1,2),c(1,3)),
                   pConf7=rbind(c(1,2),c(3,4)))                    # Parents configurations
    offsConfs=list(oConf1=t(c(1,1)),
                   oConf2=t(c(2,2)),
                   oConf3=rbind(c(1,1),c(1,2)),
                   oConf4=rbind(c(1,2),c(1,3)),
                   oConf5=rbind(c(1,1),c(2,2),c(1,2)),
                   oConf6=rbind(c(1,1),c(1,2),c(1,3),c(2,3)),
                   oConf7=rbind(c(1,3),c(1,4),c(2,3),c(2,4)))      # Offspring configurations
    parsIndex=parsConfs[[conf]]                                    # Parents index
    offsIndex=offsConfs[[conf]][sample(1:nucConf[conf],nOffs),,drop=FALSE]  # Offspring index
    H=popHaplos[sample(1:nHaplos,nuhConf[conf]),]                  # Family haplotypes from population
    familyHaps=NULL
    if (conf>1) familyHaps <- H[t(rbind(parsIndex,offsIndex)),]    # Family haplotypes when configuration 1 occurs
    else for (i in 1:((2+nCh)*2)) familyHaps=rbind(familyHaps,H)   # Family haplotypes for the rest of configurations
    dimnames(familyHaps)=NULL                                      # Initialization of the haplotype names
    oddHaps=familyHaps[seq(1,nrow(familyHaps),2),]                 # Odd haplotypes
    evenHaps=familyHaps[seq(2,nrow(familyHaps),2),]                # Even haplotypes
    oddHapsColl=apply(oddHaps,1,paste,collapse="-")                # Odd haplotypes collapse
    evenHapsColl=apply(evenHaps,1,paste,collapse="-")              # Even haplotypes collapse
    allsMat=NULL; for (i in 1:nrow(oddHaps)) {
      allsRow=NULL; for (j in 1:ncol(oddHaps)) 
        allsRow=cbind(allsRow,t(c(oddHaps[i,j],evenHaps[i,j])))
      allsMat=rbind(allsMat,allsRow)
    }                                                              # Matrix of collapsed haplotypes
    mkrsHaps=cbind(allsMat,0,oddHapsColl,evenHapsColl)             # Array with parents and offspring markers and haplotypes
    nf=nOffs+2                                                     # Number of family members
    indID=1:nf                                                     # Each individual is assigned an identification index indid
    sex=c(1,2,sample(1:2,nOffs,replace=TRUE))                      # Random assignation of sex to offsprings
    phen=rbinom(nf,1,phenProb)                                     # Random assignation of affection status to each subject
    patID=c(0,0,rep(1,nOffs))                                      # Paternal ID
    matID=c(0,0,rep(2,nOffs))                                      # Maternal ID
    fam=cbind(famID=famID,indID=indID,patID=patID,matID=matID,
              sex=sex,phen=phen,mkrsHaps,deparse.level=0)          # Data.frame composed by: famid,indid,patID,matID,sex,phenotype,markers (ncol: 5..nmarkers*2)
    return(fam)
  }
  simRecombHap=function(haps){                           # Simulation of the haplotype recombination
    hl=ncol(haps)                                        # Number of alleles in haplotypes
    breakPoint=sample(2:hl,1)                            # Position in which recombination takes place
    newHaps=haps                                         # Initialize the new recombinated haplotype as the original haplotype
    newHaps[1,breakPoint:hl]=haps[2,breakPoint:hl]       # Recombination
    newHaps[2,breakPoint:hl]=haps[1,breakPoint:hl]       # Recombination
    return(newHaps)                                      # Returns the recombinated haplotypes
  }
}
  ########################################### II. ALLELES PER MARKER #########################################
{
  ### Simulation of the number of alleles per marker if they are not supplied by user
  if (is.null(numAllperMrk)){
    if (!chrAlleles)                                                # Alleles are not character type
      numAllperMrk <- sample(2:9,nMarkers,replace=TRUE)             # Alleles range per marker
    else numAllperMrk <- rep(2,nMarkers)                            # Repeats 2 times if alleles are character type
  }
  else if (length(numAllperMrk)>=1){                                # If the number of alleles per marker is specified 
    chrAlleles=FALSE                                                # using a numeric value or a set of values (distributed)
    if (length(numAllperMrk)==1)                                    # If the no. of alleles per marker is a numeric value
      numAllperMrk <- rep(numAllperMrk,nMarkers)                    # Vector of numAllperMrk according to nMarkers
    else nMarkers <- length(numAllperMrk)                           # Number of markers according to numAllperMrk
  }
}  
  ####################################### III. HAPLOTYPES IN POPULATION ######################################
{  
  ### Generation the haplotypes which may exist in the population
  popHaplos <- simHapSelection(nHaplos,numAllperMrk)
  ### Adaptation of corresponding haplotypes when alleles are character type 
  if (chrAlleles){                                                  # If alleles are character type
    chrMrkLabel <- replicate(nMarkers,labelMrk())                   # 'A','C','G','T' character labels
    ph <- popHaplos                                                 # ph initialization
    for (j in 1:nMarkers){                                          # Scanning according to the number markers
      chl <- chrMrkLabel[,j]                                        # Character label of the marker
      popHaplos[,j] <- chl[ph[,j]+1]                                # overwriting of previous labels in the corresponding column of popHaplos 
    }    
  }
}  
  ########################################### IV. DATA CONCATENATION #########################################
{  
  ### Concatenation of the previous simulated data 
  families <- sapply(1:nFams,simOneFamily,nChildren,popHaplos,
                     phenProb,recombRate,simplify=FALSE)            # Simulation of one family
  families <- do.call(rbind, families)                              # Concatenation of the simulated families
  families <- data.frame(families)                                  # Saving of previous families as data.frame
}  
  ############################################## V. DATA LABELLING ###########################################
{  
  ### Labelling of previous generated data
  i1 <- gl(nMarkers,2)                                              # Grade levels to name the first allele column
  i2 <- gl(2,1,2*nMarkers)                                          # Grade levels to name the second allele column
  markerCols <- 6+1:(2*nMarkers)                                    # Markers' columns
  hapCols <- ncol(families)-(1:0)                                   # Haplotypes' columns 
  names(families) <- c("famID","indID","patID","matID","sex","phenot",paste("Mk",i1,"_",i2,sep=""),
                       "recombNr","Paternal_Hap","Maternal_Hap")    # Names of family columns 
}  
  ############################################# VI. DATA CONVERSION ##########################################
{  
  ### Conversion of the previous generated data into the most suitable type (integer and/or character)
  toInteger=function(f) as.numeric(levels(f))[f]                    # Integer conversion of factor levels
  toCharacter=function(f) as.character(f)                           # Character conversion of factor levels
  for (k in c(1:6,(ncol(families)-2)))                              # Scanning of the first 4 columns and recombNr column
    families[[k]] <- toInteger(families[[k]])                       # Integer conversion of famid, indid, sex, phenotype and recombNr
  if (chrAlleles) for (k in markerCols)                             # If alleles are of character type
    families[[k]] <- toCharacter(families[[k]])                     # Character conversion of the markers' columns
  else for (k in markerCols)                                        # If alleles are of numeric type
    families[[k]] <- toInteger(families[[k]])                       # Integer conversion of the markers' columns
  for (k in hapCols) families[[k]] <- as.character(families[[k]])   # Haplotype columns converted to character
}
  ######################################### VII. MISSING DATA GENERATION #####################################
{  
  ### Insertion of missing values in the previous generated data
  indivIDs <- families[,2]                                          # IDs of the families' individuals 
  parIDs <- which(indivIDs<3); offIDs <- which(indivIDs>2)          # IDs of parents and offspring
  ugp <- rbinom(length(parIDs),1,ungenotPars)                       # Binomial variable with 1 being the row completely ungenotyped
  ugo <- rbinom(length(offIDs),1,ungenotOffs)                       # Same with the offspring
  completeLostPrs <- parIDs[ugp==1]                                 # ID of parents completely ungenotyped
  completeLostOfs <- offIDs[ugo==1]                                 # ID of offsprings completely ungenotyped
  restPrts <- parIDs[ugp==0]                                        # ID of rest of parents
  restOffs <- offIDs[ugo==0]                                        # ID of rest of children
  lostRows <- matrix(1,nrow=sum(c(ugp,ugo)),ncol=nMarkers)          # Matrix for completely ungenotyped markers (parents and children)
  lostPrsMkrs <- matrix(rbinom(length(restPrts)*nMarkers,1,missParProb),ncol=nMarkers)  # Matrix for ungenotyped markers in parents
  lostOfsMkrs <- matrix(rbinom(length(restOffs)*nMarkers,1,missOffProb),ncol=nMarkers)  # Matrix for ungenotyped markers in offsprings
  lost <- rbind(lostRows,lostPrsMkrs,lostOfsMkrs)                   # Matrix with 1's being lost markers
  lost <- lost[order(c(completeLostPrs,completeLostOfs,restPrts,restOffs)),]  # The same matrix ordered
  lost <- lost[,gl(nMarkers,2)]                                     # The same matrix with columns duplicated (each marker comprises two columns)
  families[,markerCols][lost==1] <- NA                              # Final statement of lost markers in the positions specified by matrix lost
}  
  ############################################ VIII. FUNCTION OUTPUT #########################################
{  
  families[,1] <- paste("FAM0",families[,1],sep="")                    # Names of the families
  datasetAlls <- families[-((length(families)-2):length(families))]    # Data.frame composed by: famid,indid,sex,phenotype,markers
  if (dataSummary==TRUE) datasetAlls <- alleLoader(datasetAlls)        # 
  datasetHaps <- families[-(7:(6+nMarkers*2+1))]                       # Data.frame composed by: famid,indid,sex,phenotype,recombNr,haplotypes
  output <- list(datasetAlls=datasetAlls,datasetHaps=datasetHaps)      # 
  if (invisibleOutput)  return(invisible(output)) else return(output)  # 
}
}