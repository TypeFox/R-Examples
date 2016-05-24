#' @title Imputation of missing alleles from a dataset composed by families.
#' @description By analyzing all possible combinations of a parent-offspring pedigree in which parental and/or offspring genotypes may be missing; as long as one child was genotyped, in certain cases it is possible an unequivocal imputation of the missing genotypes both in parents and children.
#' @param data Data containing the families' identifiers and the corresponding genetic data (or the path of the PED file).
#' @param invisibleOutput Data are not shown by default.
#' @param dataSummary A summary of the data is shown by default.
#' @return Imputed markers, Homozygosity (HMZ) matrix, marker messages and number of unique alleles per marker.
#' @import utils
#' @export alleImputer
#' @references Medina-Rodriguez, N. Santana A. et al. (2014) alleHap: an efficient algorithm to reconstruct zero-recombinant haplotypes from parent-offspring pedigrees. BMC Bioinformatics, 15, A6 (S-3).
#' @examples
#' 
#' ## Imputation of families containing parental missing data
#' simulatedFams <- alleSimulator(10,4,6,missParProb=0.2,dataSummary=FALSE) 
#' famsAlls <- simulatedFams[[1]]       # Alleles (genotypes) of the simulated families
#' alleImputer(famsAlls)                # Imputed alleles (genotypes)
#' 
#' ## Imputation of families containing offspring missing data
#' datasetAlls <- alleSimulator(10,4,6,missOffProb=0.2,dataSummary=FALSE)
#' famsAlls <- simulatedFams[[1]]       # Alleles (genotypes) of the simulated families
#' alleImputer(famsAlls)                # Imputed alleles (genotypes)
#' 
#' ## Imputation of a family marker containing missing values in one parent and one child
#' infoFam <- data.frame(famID="FAM03",indID=1:5,patID=c(0,0,1,1,1),
#'                       matID=c(0,0,2,2,2),sex=c(1,2,1,2,1),phenot=0)
#' mkr <- rbind(father=c(NA,NA),mother=c(1,3),child1=c(1,1),child2=c(2,3),child3=c(NA,NA))
#' colnames(mkr) <- c("Mkr1_1","Mkr1_2")
#' famMkr <- cbind(infoFam,mkr)         # Alleles (genotypes) of the family
#' alleImputer(famMkr)                  # Imputed alleles (genotypes)
#' 
alleImputer=function(data, invisibleOutput=TRUE, dataSummary=TRUE){
  
  ############################################### I. INTERNAL FUNCTIONS ##########################################
  {
    # Imputation of a marker (main function)
    mkrImputer=function(mkr){  
      #=== INTERNAL FUNCTIONS ===#
      classSet=function(mkr){                              # Convertion of the marker into a factor type  
        mkr <- apply(mkr,2,as.character)                   # Contertion to character
        dfm <- dim(mkr)                                    # Dimensions
        rnfm <- rownames(mkr)                              # Row names
        mkr <- factor(mkr,exclude=c(-9,NA,"NA","<NA>"))    # Factor
        dim(mkr) <- dfm                                    # Reassign dimensions
        rownames(mkr) <- rnfm                              # Reassign row names
        lfm <- levels(mkr)                                 # Levels
        return(list(mkr=mkr,dfm=dfm,lfm=lfm,rnfm=rnfm))    # List containing the converted marker, dimensions, levels, etc.
      }
      commAlls=function(x,y){                              # Checks the common alleles
        if (all(is.na(x))|all(is.na(y))) return(TRUE)       
        else if(length(intersect(x,y))>0) return(TRUE)        
        else return(FALSE)                                  
      }
      #=== QUALITY CONTROL ===#
      {  
        {
          # Verification of the unique non-missing alleles per marker
          lcs <- classSet(mkr)                                                 # Convertion of the marker alleles to factor type
          fmkr <- lcs$mkr;  lfm <- lcs$lfm; dfm <- lcs$dfm; rnfm <- lcs$rnfm   # Alleles, levels, dimensions and names of the marker 
          famSize <- nrow(fmkr)                                                # Family size (number of individuals of the family)
          homz <- as.integer(fmkr[,1]==fmkr[,2])                               # Identification of the homozygous allele pairs 
          fmkr <- cbind(fmkr,homz)                                             # Family marker with the homozygous identificator attached  
          alNr <- unique(as.vector(fmkr[,1:2]))                                # Number of unique alleles per marker
          alNr <- length(alNr[!is.na(alNr)])                                   # Number of unique alleles per marker excluding the missing values
          pph <- as.matrix(expand.grid(1:2,3:famSize,KEEP.OUT.ATTRS=TRUE),2)   # Indexation of the family members
          compat <- all(apply(pph,1, function(ind) 
                        commAlls(fmkr[ind[1],1:2],fmkr[ind[2],1:2])))          # Compatibility vector
          parents <- fmkr[1:2,]                                                # Parents marker plus homz column
          children <- fmkr[3:famSize,,drop=FALSE]                              # Children marker plus homz column
          childrenNoNA <- children[!is.na(children[,3]),,drop=FALSE]           # Children marker excluding individuals with missing values
          parAlls <- unique(as.vector(parents[,1:2]))                          # Alleles in parents
          childAllsNoNA <- unique(as.vector(childrenNoNA[,1:2]))               # Alleles in children excluding individuals with missing values
          unidAlls <- childAllsNoNA[(!childAllsNoNA%in%parAlls)]               # Children alleles which have been non-identified in parents
          whichHMZpars <- which(homz[1:2]==1)                                  # Homozygous parents
          whichHTZpars <- which(homz[1:2]==0)                                  # Heterozygous parents
          whichHMZchildNoNA <- which(homz[-(1:2)]==1)                          # Which are the homozygous children (excluding NA individuals)
          HMZchildrenNoNA <- children[whichHMZchildNoNA,1:2,drop=FALSE]        # Homozygous children (excluding NA individuals)
          uniqHMZchildNoNA <- unique(t(apply(HMZchildrenNoNA,1,sort)))         # Unique homozygous children (excluding NA individuals)
          nuHMZc <- nrow(uniqHMZchildNoNA)                                     # Number of unique homozygous children (excluding NA individuals)
          HMZallsNoNA <- if (length(whichHMZchildNoNA)>=1) 
            unique(as.vector(children[whichHMZchildNoNA,1:2])) else NULL       # Alleles in the homozygous children
          uniqueChildrenNoNA <- if (nrow(childrenNoNA)==1) childrenNoNA
          else unique(t(apply(childrenNoNA[,1:2],1,sort)))                     # Unique children (excluding NA individuals)
          nuc <- nrow(uniqueChildrenNoNA)                                      # Number of unique children (excluding NA individuals)
          whichHTZchildNoNA <- which(childrenNoNA[,3]==0)                      # Which are the heterozygous children (excluding NA individuals)
          HTZchildrenNoNA <- childrenNoNA[whichHTZchildNoNA,1:2,drop=FALSE]    # Heterozygous children (excluding NA individuals)
          uniqHTZchildNoNA <- unique(t(apply(HTZchildrenNoNA,1,sort)))         # Unique heterozygous children (excluding NA individuals)
          nuHTZc <- nrow(uniqHTZchildNoNA)                                     # Number of unique heterozygous children (excluding NA individuals)
          uHTZsum <- apply(uniqHTZchildNoNA,2,sum)                                 # Sum of unique heterozygous children alleles (per columns)
          missPars <- which(is.na(homz[1:2]))                                  # Missing genotype in parents 
          message <- NULL                                                      # Initialization of the marker message
        }
        ### Q.C. a) If each child has one common allele (at least) with each parent
        if (!compat) { 
          message <- "a"; fmkr[,] <- NA 
        } 
        #-- Considering parents without NA values
        if (length(missPars)==0) { 
          ### Q.C. b) If the number of unique alleles is 4 as maximum.
          if (alNr>4) { 
            message <- "b"; fmkr[,] <- NA 
          } 
          ### Q.C. c) Checks if a child has alleles not present in parents.
          else if (length(unidAlls)>0) {
            message <- "c"; fmkr[,] <- NA 
          }
          ### Q.C. d) There cannot be more than two unique homozygous children.
          else if (length(HMZallsNoNA)>2) { 
            message <- "d"; fmkr[,] <- NA
          }
        } 
        #-- Considering parents with NA values
        else {
          ### Q.C. e) If there are three or more unique heterozygous children and if they share the same allele.
          if (nuHTZc>=3&uHTZsum[1]==nuHTZc) {
            message <- "e"; fmkr[,] <- NA
          }
          ### Q.C. f) If one parent has missing genotype
          else if (length(missPars)==1) {
            if (length(whichHTZpars)>0&nuHMZc>2) {
              message <- "f1"; fmkr[,] <- NA
            }
            else if (length(whichHTZpars)>0&alNr==4&nuHMZc>=1) {
              message <- "f2"; fmkr[,] <- NA
            }
            else if (length(whichHMZpars)>0&nuc>2) {
              message <- "f3"; fmkr[,] <- NA
            }
          }
          ### Q.C. f) If both parents have missing genotype
          else if (length(missPars)==2) {
            if (nuc>4) {
              message <- "g1"; fmkr[,] <- NA
            }
            else if (alNr==4&nuHMZc>0) {
              message <- "g2"; fmkr[,] <- NA
            }
          }
        
        }
      }  
      #=== IMPUTATION PHASE ===#
      {  
        ### I.P. a) Children Imputation: Missing alleles in children are imputed only if a parent is homozygous 
        if (is.null(message)) {                                              # If there are irregular children  
          whichNAchildren <- which(is.na(children[,3]))                      # Children with missing values (NA values)      
          whichHOMZpar <- which(parents[,3]==1)                              # Parents with homozygous alleles
          if (length(whichNAchildren)>=1&length(whichHOMZpar)>=1){           # If there Na values in children and at least 
            for (k in whichNAchildren) {                                     # a parent is homozygous
              for (j in whichHOMZpar) children[k,j] <- parents[j,1]          # Imputation of the parent allele into the child
              children[k,3] <- as.integer(children[k,1]==children[k,2])      # Updating the homz value if child alleles are equal
            }
          }
        } 
        ### I.P. b) Parent Imputation
        else {
          #-- The alleles in any homozygous child are imputed in parents
          whichNApar <- which(is.na(parents[,3]))                              # Parents with missing values
          whichHMZchild <- which(children[,3]==1)                              # Homozygous children
          aleHOM <- if (length(whichHMZchild)>=1)                              # Unique homozygous alleles in children
            unique(as.vector(children[whichHMZchild,1:2])) else NULL            
          if (length(whichNApar)>=1&length(whichHMZchild)>=1&!is.null(aleHOM)){ 
            if (length(aleHOM)==2)                                             # If there are tow unique homozygous alleles
              for (k in whichNApar) parents[k,] <- c(aleHOM,FALSE)             # Imputation in both parents
            else for (k in whichNApar) parents[k,1] <- aleHOM                  # Imputation in one parent
          }
          #-- If one parent has missing alleles and the other not, and the heterozygous children are not present
          #   in the complete parent, those alleles are imputed to the other parent.
          whichNApar <- which(is.na(parents[,3]))                              # Parent with missing values (NA values) 
          if (length(whichNApar)==1){
            completePar <- 3-whichNApar                                        # Complete parent
            whichHETchild <- which(children[,3]==FALSE)                        # Children with heterozygous alleles
            aleHET <- unique(as.vector(children[whichHETchild,1:2]))           # Heterozygous alleles in children
            ale2par <- aleHET[!aleHET %in% parents[completePar,1:2]]           # Allele/s to impute in parent
            if (length(ale2par)>=1){                                           # If there are one or more alleles to impute
              lwn <- which(is.na(parents[whichNApar,1:2]))                     # Positions of the NA alleles in parents
              if (length(lwn)==1)                                              # If there is one NA allele in parents
                ale2par <- unique(c(parents[whichNApar,3-lwn],ale2par))        # Update of the parent allele
              if (length(ale2par)>2){                                          # If there are more than two alleles to impute
                message <- "h"                                                 # 
                parents[,] <- NA;  children[,] <- NA                           # Update to NA
              }
              else if (length(ale2par)==2)                                     # If there are two alleles to impute in parent
                parents[whichNApar,] <- c(ale2par,FALSE)                       # Update the alleles to impute in the NA parent and put homz to FALSE
              else parents[whichNApar,1] <- ale2par                            # Update the alleles to impute in the NA parent 
            }
          }
        }
      }  
      #=== FUNCTION OUTPUT ===#
      {  
        if (is.null(message)) {                                            # If there is no message
          mk <- rbind(parents,children)                                    # Children and parents imputed alleles are joined
          marker <- factor(mk[,1:2],labels=lfm)                            # Re-assignation of the labels containing the new imputed values
          dim(marker) <- dfm; rownames(marker) <- rnfm                     # Update of dimensions and row names
          return(list(marker=marker,allelesNumber=alNr))                   # Return the imputed marker, messages and the number of alleles
        } else return(list(marker=fmkr,message=message,allelesNumber=NA))  # Return the imputed marker and the number of alleles
      } 
    }
    # Imputation of a family
    famImputer=function(family,nMkrs,mkCols){
      imputedMkrs <- family[,-(1:6)]                                       # Extraction of the family alleles
      allelesNumber <- vector("numeric",nMkrs)                             # Initialization of allelesNumber 
      mkrIncidences <- vector("character",nMkrs)                           # Initialization of mkrIncidences
      for (k in 1:nMkrs) {                                                 # Scanning of each marker
        fkCols <- mkCols[,k]                                               # Pair of columns of each marker
        mkr <- imputedMkrs[,fkCols]                                        # Alleles of each marker
        fmrk <- mkrImputer(mkr)                                            # Imputation of each marker
        imputedMkrs[,fkCols] <- fmrk$marker[,1:2]                          # Alleles of each imputed marker
        allelesNumber[k] <- fmrk$allelesNumber                             # Number of unique alleles of each imputed marker
        if (exists("message",fmrk)) mkrIncidences[k] <- fmrk$message       # Saving the messages of each imputed marker
      }
      fam <- list(imputedMkrs=cbind(family[,1:6],imputedMkrs),
                  mkrIncidences=mkrIncidences) # 
      return(fam)                                                          # 
    }
    # Imputation of multiple families
    famsImputer=function(datasetAlls){
      fams <- list(imputedMkrs=NULL,mkrIncidences=NULL)                            # Initialization of the list of families
      idFams <- unique(datasetAlls[,1])                                            # Number of identification of the families
      nMkrs <- (ncol(datasetAlls)-6)/2                                             # Number of markers 
      mkCols <- matrix(1:(nMkrs*2),nrow=2)                                         # Columns of the biallelic markers
      for (f in idFams) {                                                          # Scanning of all family identifiers
        family <- datasetAlls[datasetAlls[,1]==f,]                                 # Family data
        family <- family[order(family[,3]),]                                       # Family sorting according to the third column (paternal ID)
        fam <- famImputer(family,nMkrs,mkCols)                                     # 
        fams$imputedMkrs <- rbind(fams$imputedMkrs,fam$imputedMkrs)                # 
        fams$mkrIncidences <- rbind(fams$mkrIncidences,
                                    c(as.character(family[1,1]),fam$mkrIncidences))# 
      }
      colnames(fams$mkrIncidences) <- c("Family",paste("Mkr",1:nMkrs,sep=""))
      fams$mkrIncidencesMessages <- matrix(c("Some children have no common alleles with a parent",
                                             "More alleles than possible in this marker",
                                             "Some children have alleles not present in parents",
                                             "Some homozygous children are not compatible in this marker",
                                             "Three or more unique heterozygous children share the same allele",
                                             "Heterozygous parent and more than two unique homozygous children",
                                             "Heterozygous parent, four unique alleles and more than one unique homozygous children",
                                             "Homozygous parent and more than two unique children",
                                             "More than four unique children geneotypes in the family",
                                             "Homozygous genotypes and four unique alleles in children",
                                             "A marker would require a parent to have three different alleles"),11)
      rownames(fams$mkrIncidencesMessages) <- c("a","b","c","d","e","f1","f2","f3","g1","g2","h")
      colnames(fams$mkrIncidencesMessages) <- "Incidence"
      fams$mkrIncidences <- as.data.frame(fams$mkrIncidences)
      return(fams)
    }
  }
  ################################################# II. DATA LOADING #############################################
  {
    datasetAlls <- alleLoader(data,invisibleOutput, dataSummary)             # 
  }
  ################################################# III. IMPUTATION  #############################################
  {
    imputationTime <- system.time(fams <- famsImputer(datasetAlls))[[3]]     #
  }
  ################################################# IV. DATA SUMMARY #############################################
  {
    nMissAlls <- length(which(is.na(datasetAlls[,-(1:6)])))                  # Number of missing alleles
    nImpMissAlls <- length(which(is.na(fams$imputedMkrs[,-(1:6)])))          # Number of missing alleles after imputation
    imputationRate <- nImpMissAlls/nMissAlls                                 # Imputation rate
    fams$imputationSummary <- data.frame(nMissAlls,nImpMissAlls,
                                         imputationRate,imputationTime)      # Data summary
  }
  #################################################  V. DATA STORING  ############################################
  {
    if (is.character(data)) {                                                # 
      baseName <- file_path_sans_ext(data)                                   # 
      outName <- paste(baseName,"imputed.ped",sep="_")                       # 
      write.table(fams$imputedMkrs,sep=" ",quote=FALSE,file=outName)         #
      if (any(fams$mkrIncidences[,-1]!="")) {                                # 
        outNameInc <- paste(baseName,"mkrIncidences.txt",sep="_")            #
        write.table(fams$mkrIncidences,sep=" ",quote=FALSE,file=outNameInc)  #  
      }
    }
  }
  ############################################### VI. FUNCTION OUTPUT  ###########################################
  {
    if (dataSummary==TRUE) {                                                 # 
      
      ### Imputation/incidences message printing
      if (all(fams$mkrIncidences[,-1]=="")) {                                # 
        cat("\nAlleles have been successfully imputed!!!")                   # 
        fams$mkrIncidences <- NULL                                           # 
        fams$mkrIncidencesMessages <- NULL                                   # 
      } else {                                                               #  
        cat("\nAlleles have been imputed, but some incidences in markers were detected.") 
      }
      
      ### Data summary printing
      cat("\n\n===== IMPUTATION SUMMARY =====")                              # 
      cat(paste("\nNumber of missing alleles:",nMissAlls))                   #
      cat(paste("\nNumber of imputed alleles:",nImpMissAlls))                # 
      cat(paste("\nImputation rate:",round(imputationRate,2)))               # 
      cat(paste("\nImputation time:",round(imputationTime,2)))               # 
      cat("\n==============================\n")
      
      ### Storage path
      if (is.character(data)) {                                              # 
        cat("\nImputed data have been stored in: \n")
        cat(paste("\n",outName,"\n",sep=""))                                 # 
        cat(paste(baseName,"_mkrIncidences.txt",sep=""))                     # 
      } #else cat(getwd())                                                    # 
      
    } 
    
    ### Cleaning all empty marker incidences
    if (all(fams$mkrIncidences[,-1]=="")) {                                  # 
      fams$mkrIncidences <- NULL                                             # 
      fams$mkrIncidencesMessages <- NULL                                     # 
    }    
    
    ### Returning (or not) the output
    if (invisibleOutput) return(invisible(fams)) else return(fams)           #     
  }

}
