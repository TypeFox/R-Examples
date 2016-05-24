#' @title Data loading of nuclear families (in .ped format)
#' @description The data to be loaded must be structured in .ped format and families must comprise by parent-offspring pedigrees.
#' @param data Data to be loaded.
#' @param invisibleOutput Data are not shown by default.
#' @param dataSummary A summary of the data is shown by default.
#' @return Loaded dataset.
#' @import tools
#' @import utils
#' @export alleLoader
#' @references Medina-Rodriguez, N. Santana A. et al. (2014) alleHap: an efficient algorithm to reconstruct zero-recombinant haplotypes from parent-offspring pedigrees. BMC Bioinformatics, 15, A6 (S-3).
#' @examples
#' 
#' ## Loading of a dataset in .ped format with alphabetical alleles (A,C,G,T)
#' example1 <- file.path(find.package("alleHap"), "examples", "example1.ped")
#' example1Alls <- alleLoader(example1)
#' head(example1Alls)
#' 
#' ## Loading of a dataset in .ped format with numerical alleles
#' example2 <- file.path(find.package("alleHap"), "examples", "example2.ped")
#' example2Alls <- alleLoader(example2)
#' head(example2Alls)
#' 
alleLoader=function(data, invisibleOutput=TRUE, dataSummary=TRUE){
  ################################### I. EXTENSION CHECK & DATA READ ####################################
{  
  if (is.character(data)) {
    filetype <- file_ext(data)
    if (filetype=="ped") datasetAlls <- read.table(data)
    else {
      cat("===========================================")
      cat("\n===== alleHap package: version",
      sessionInfo("alleHap")$otherPkgs$alleHap$Version,"======")
      cat("\n===========================================")
      cat("\n\nYour data could not be loaded!!!") 
      return(stop("The file must have a .ped extension"))
    }
  } else if (is.data.frame(data)|is.matrix(data)) {
    datasetAlls <- data
  } else {
    cat("===========================================")
    cat("\n===== alleHap package: version",
    sessionInfo("alleHap")$otherPkgs$alleHap$Version,"======")
    cat("\n===========================================")
    cat("\n\nYour data could not be loaded!!!")
    return(stop("Your data must be in data.frame or matrix format"))
  }
}
  ########################################## II. DATA CHECK #############################################
{
  ### Data identifiers
  famIDs <- datasetAlls[,1]                                   # Family identifiers
  indIDs <- datasetAlls[,2]                                   # Individual identifiers
  patIDs <- datasetAlls[,3]                                   # Paternal identifiers
  matIDs <- datasetAlls[,4]                                   # Maternal identifiers
  sexIDs <- datasetAlls[,5]                                   # Sex identifiers
  phenot <- datasetAlls[,6]                                   # Phenotype values
  famAlleles <- datasetAlls[,7:ncol(datasetAlls)]             # Family Alleles
  
  ### Data Counting
  nFams <- length(unique(famIDs))                             # Number of families
  nIndivs <- length(indIDs)                                   # Number of individuals
  nPars <- length(which(indIDs<3))                            # Number of parents (founders)
  nOffs <- length(which(indIDs>2))                            # Number of children (offspring)
  nMals <- length(which(sexIDs%in%c(1,"male")))               # Number of males
  nFems <- length(which(sexIDs%in%c(2,"female")))             # Number of females
  nMrks <- (ncol(datasetAlls)-6)/2                            # Number of markers
  
  ### Data Ranges
  if (length(unique(famIDs))>1) {                             # If there is more than one family
    minFamID <- min(famIDs); maxFamID <- max(famIDs)          # Range of Family IDs
  }
  minIndID <- min(indIDs); maxIndID <- max(indIDs)            # Range of Individual IDs
  levsPatID <- levels(as.factor(patIDs))                      # Paternal IDs
  levsMatID <- levels(as.factor(matIDs))                      # Maternal IDs
  levsSex <- levels(as.factor(sexIDs))                        # Sex values
  if (nlevels(as.factor(phenot))<=5) {                        # Phenotype values
    levsPhen <- levels(as.factor(phenot))                     # Phenotype levels
  } else { 
    minPhen <- min(phenot); maxPhen <- max(phenot)            # Min and max phenotype values 
  }
}
  ###################################### III. MISSING DATA CHECK ########################################
{  
  noPars <- length(which(datasetAlls[indIDs<3,2]%in%c(0,-9))|
                   which(is.na(datasetAlls[indIDs<3,2])))                  # Unknown parents
  noOffs <- length(which(datasetAlls[indIDs>2,2]%in%c(0,-9))|
                   which(is.na(datasetAlls[indIDs>2,2])))                  # Unknown offspring
  noPatID <- length(which(patIDs==0))                                      # Unknown Paternal IDs
  noMatID <- length(which(matIDs==0))                                      # Unknown Maternal IDs
  noSex <- length(which(!sexIDs%in%c(1,2,"male","female")))                # Unknown Sex
  noPhen <- length(which(phenot==-9))                                      # Unknown Phenotypes
  naAlls <- length(which(is.na(famAlleles)))                               # Unknown Genotypes
  pos09Alls <- which(famAlleles%in%c(0,-9))                                # Positions of 0 and -9 values
  na09Alls <- length(pos09Alls)                                            # Number of 0 and -9 values
  missAlls <- naAlls + na09Alls                                            # Unknown Alleles
  if (na09Alls>0) datasetAlls[,round(pos09Alls/nrow(datasetAlls))+7] <- NA # Genotype missing values
  mkrNA <- function(j) if (any(is.na(datasetAlls[,j]))) 1 else 0
  odColMrk <- seq(7,ncol(datasetAlls),by=2)                                # Odd columns
  evColMrk <- seq(8,ncol(datasetAlls),by=2)                                # Even columns
  NAmkrs <- length(which(colSums(rbind(sapply(odColMrk,mkrNA),
                                         sapply(evColMrk,mkrNA)))>0))      # Markers containing missing values
}
  ######################################## IV. FUNCTION OUTPUT ##########################################
{
  if (dataSummary==TRUE) {
    cat("===========================================")
    cat("\n===== alleHap package: version",
    sessionInfo("alleHap")$otherPkgs$alleHap$Version,"======")
    cat("\n===========================================")
    cat("\n\nData have been successfully loaded from: \n")
    if (is.character(data)) cat(data) else cat(getwd())
    cat("\n\n===== DATA COUNTING ======")
    cat(paste("\nNumber of families:",nFams))
    cat(paste("\nNumber of individuals:",nIndivs))
    cat(paste("\nNumber of founders:",nPars))
    cat(paste("\nNumber of children:",nOffs))
    cat(paste("\nNumber of males:",nMals))
    cat(paste("\nNumber of females:",nFems))
    cat(paste("\nNumber of markers:",nMrks))
    cat("\n===========================")
    cat("\n\n======== DATA RANGES =========")
    if (length(unique(famIDs))==1) cat(paste("\nFamily ID: ",unique(famIDs),sep=""))
    else cat(paste("\nFamily IDs: [",minFamID,",...,",maxFamID,"]",sep=""))
    cat(paste("\nIndividual IDs: [",minIndID,",...,",maxIndID,"]",sep=""))
    cat(paste("\nPaternal IDs: [",paste(levsPatID,collapse=","),"]",sep=""))
    cat(paste("\nMaternal IDs: [",paste(levsMatID,collapse=","),"]",sep=""))
    cat(paste("\nSex values: [",paste(levsSex,collapse=","),"]",sep=""))
    if (nlevels(as.factor(phenot))<=5) 
      cat(paste("\nPhenotype values: [",paste(levsPhen,collapse=","),"]",sep=""))
    else cat(paste("\nRange of Phenotypes: [",minPhen,",...,",maxPhen,"]",sep=""))
    cat("\n==============================")
    cat("\n\n===== MISSING DATA =====")
    cat(paste("\nMissing founders:",noPars))
    cat(paste("\nMissing children:",noOffs))
    cat(paste("\nMissing paternal IDs:",noPatID))
    cat(paste("\nMissing maternal IDs:",noMatID))
    cat(paste("\nMissing sex:",noSex))
    cat(paste("\nMissing phenotypes:",noPhen))
    cat(paste("\nMissing alleles:",missAlls))
    cat(paste("\nNumber of NA markers:",NAmkrs))
    cat("\n========================\n")
  }
  if (invisibleOutput) return(invisible(datasetAlls)) else return(datasetAlls)
}
}
