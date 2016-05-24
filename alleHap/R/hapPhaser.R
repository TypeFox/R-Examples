#' @title Phasing of a dataset composed by several families.
#' @description By analyzing all possible combinations of a parent-offspring pedigree in which parents may be missing (missParProb>0), as long as one child was genotyped, it is possible an unequivocal reconstruction of many parental haplotypes. When neither parent was genotyped (missParProb==1), also it is possible to reconstruct at least two parental haplotypes in certain cases. Regarding offspring haplotypes, if both parents are completely genotyped (missParProb==0), in majority of cases partial offspring haplotypes may be successfully phased (missOffProb>0).
#' @param data Data containing non-genetic and genetic information of families (or PED file path).
#' @param NAsymbol Icon which will be placed in the NA values of the haplotypes.
#' @param alleSep Icon which will be used as separator of the haplotype alleles.
#' @param invisibleOutput Data are not shown by default.
#' @param dataSummary A summary of the data is shown by default.
#' @return Phased markers and haplotypes for each loaded family.
#' @import abind
#' @import stats
#' @import utils
#' @export hapPhaser
#' @references Medina-Rodriguez, N. Santana A. et al. (2014) alleHap: an efficient algorithm to reconstruct zero-recombinant haplotypes from parent-offspring pedigrees. BMC Bioinformatics, 15, A6 (S-3).
#' @examples
#' 
#' ## Haplotype reconstruction for 3 families without missing data.
#' simulatedFams <- alleSimulator(3,3,6,dataSummary=FALSE)  
#' (famsAlls <- simulatedFams[[1]])      # Alleles (genotypes) of the simulated families
#' phasedFams <- hapPhaser(famsAlls)     # List containing the phased families
#' phasedFams$phasedMkrs                 # Phased markers
#' phasedFams$haplotypes                 # phased haplotypes
#' 
#' ## Haplotype reconstruction of a family containing missing data in a parent. 
#' infoFam <- data.frame(famID="FAM002",indID=1:6,patID=c(0,0,1,1,1,1),
#'                      matID=c(0,0,2,2,2,2),sex=c(1,2,1,2,1,2),phenot=0)
#' Mkrs <- rbind(c(1,4,2,5,3,6),rep(NA,6),c(1,7,2,3,3,2),
#'               c(4,7,5,3,6,2),c(1,1,2,2,3,3),c(1,4,2,5,3,6))
#' colnames(Mkrs) <- c("Mk1_1","Mk1_2","Mk2_1","Mk2_2","Mk3_1","Mk3_2")
#' (family1 <- cbind(infoFam,Mkrs))
#' phasedFam1 <- hapPhaser(family1)      # List containing the phased family
#' phasedFam1$phasedMkrs                 # Phased markers
#' phasedFam1$haplotypes                 # Phased haplotypes
#'
hapPhaser=function(data,NAsymbol="?",alleSep="",invisibleOutput=TRUE,dataSummary=TRUE){

  ############################################# I. INTERNAL FUNCTIONS ###########################################
  { 
    # Haplotype phasing per family (main function)
    famPhaser=function(fam){
      
      # Internal Functions
      as.familyMatrix=function(familyArray){
        outArray <- data.frame(t(array(familyArray,dim=c(2*dim(familyArray)[2],dim(familyArray)[3]))))
        nm <- ncol(outArray)/2
        names(outArray) <- paste("Mk",gl(nm,2),"_",gl(2,1,nm*2),sep="")
        return(outArray)
      }
      childrenHapsFromParents=function(fam){
        for (child in 3:fam$famSize){
          for (p in 1:2){
            idh=fam$idHaps[child,p]  # Parent's haplotype inherited by child      
            if (idh!=0){
              parentHap=fam$phasedAlls[idh,,p]
              parentHap[which(fam$IDS[idh,,p]==0)]=NA
              uh=updateHap(parentHap,subject=list(phasedAlls=fam$phasedAlls[,,child],IDS=fam$IDS[,,child]),
                           hapIndex=p)
              if (!is.null(uh)){
                fam$phasedAlls[,,child]=uh$phasedAlls
                fam$IDS[,,child]=uh$IDS
                if (fam$idHaps[child,3-p]==0)
                  fam$idHaps[child,3-p]=inheritedHap(fam$phasedAlls[3-p,,child],fam$phasedAlls[,,3-p],
                                                     fam$IDS[3-p,,child],fam$IDS[,,3-p])            
              }
            }     
          }
        }
        return(fam)
      }
      inheritedHap=function(hap,parentHaps,hapIDS,pHapsIDS){   # Generation of idHaps matrix
        # This function determines which one of parentHaps is (or can be) equal to 'hap'
        # it returns 0 if both parentHaps are compatible with 'hap' and 3 if there is no compatibility
        if(is.null(dim(parentHaps))){ 
          dim(hap) <- dim(hapIDS) <- c(1,1)
          dim(parentHaps) <- dim(pHapsIDS) <- c(2,1)
        } 
        canBeSameHap=vector("logical",2)
        for (i in 1:2){
          pos=which(hapIDS+pHapsIDS[i,]==2)   
          if (length(pos)==0) canBeSameHap[i]=TRUE
          else canBeSameHap[i]=all(hap[pos]==parentHaps[i,pos])  
          if (length(pos)==ncol(pHapsIDS)&canBeSameHap[i]) return(i)
        }
        inhrtdHap=which(canBeSameHap)
        if (length(inhrtdHap)==2) inhrtdHap=0
        else if (length(inhrtdHap)==0) inhrtdHap=3  # Irregular inheritance
        return(inhrtdHap)
      }
      parentsHapsFromChildren=function(fam){ 
        phasedAlleles=function(subject,hapIndex){  # Find alleles correctly phased in the hap "hapIndex" of subject
          ## Alleles not phased are coded as NA
          phsdalls=fam$phasedAlls[hapIndex,,subject,drop=FALSE]
          phsdalls[fam$IDS[hapIndex,,subject]==0]=NA
          return(phsdalls)
        }
        # Update of parental haplotypes from the children information
        updatedParents <- NULL
        for (parent in 1:2){
          incompleteHaps <- which(apply(fam$IDS[,,parent],1,sum)<fam$nMkrs)
          if (length(incompleteHaps)>0){
            # Detects which children have its parent's hap unassigned (i.e. it is not known if the inherited
            # parent's hap is the first or the second haplotype of the parent)
            revCh <- 2+which(fam$idHaps[-(1:2),parent]==0)
            nrIdAlls <- apply(fam$IDS[parent,,revCh,drop=FALSE],3,sum)
            # Detect if any of this children has all its alleles phased (i.e. all IDS values equal to 1).
            # This means that a parent haplotype is fully identified:
            fullIdCh <- 2+which(nrIdAlls==fam$nMkrs)
            # If there is one such child, the haplotype is identified in parent:
            if (length(fullIdCh)>0){
              parentHap <- fam$phasedAlls[parent,,fullIdCh[1]]
              uh <- updateHap(haplo=parentHap,
                              subject=list(phasedAlls=fam$phasedAlls[,,parent],IDS=fam$IDS[,,parent]),hapIndex=1)
              if (!is.null(uh)){
                updatedParents <- c(updatedParents,parent)
                fam$phasedAlls[,,parent] <- uh$phasedAlls
                fam$IDS[,,parent] <- uh$IDS
                for (ch in 3:fam$famSize) {
                  fam$idHaps[ch,parent] <- inheritedHap(fam$phasedAlls[parent,,ch],fam$phasedAlls[,,parent],
                                                        fam$IDS[parent,,ch],fam$IDS[,,parent]) 
                }
                incompleteHaps <- which(apply(fam$IDS[,,parent],1,sum)<fam$nMkrs)
              }
            } 
            # Now if there is still any incomplete haplotype in parent, the algorithm searches in the
            # children with that haplotype and try to phase (or impute) the unidentified alleles:
            for (hap in incompleteHaps){
              parentHap <- phasedAlleles(parent,hap)
              idCh <- 2+which(fam$idHaps[-(1:2),parent]==hap)
              if (length(idCh)>0) for (i in idCh) parentHap <- rbind(parentHap,phasedAlleles(i,parent))
              parentHap <- apply(parentHap,2,function(x) if (all(is.na(x))) NA else unique(x[!is.na(x)]))
              uh <- updateHap(parentHap,subject=list(phasedAlls=fam$phasedAlls[,,parent],IDS=fam$IDS[,,parent]),
                              hapIndex=hap)
              if (!is.null(uh)){
                updatedParents <- c(updatedParents,parent)
                fam$phasedAlls[,,parent] <- uh$phasedAlls
                fam$IDS[,,parent] <- uh$IDS
              }
            }
          }
        }
        
        # If any child has haplotypes with undefined origin, check if now origin can be univocally determined 
        if (!is.null(updatedParents)){
          updatedParents <- unique(updatedParents)
          undefChHaps <- which(fam$idHaps[-(1:2),]==0,arr.ind=TRUE)
          if (length(undefChHaps)>0) {
            undefChHaps=matrix(undefChHaps,ncol=2)
            for (i in 1:nrow(undefChHaps)){
              ch <- undefChHaps[i,1]+2
              fromParent <- undefChHaps[i,2]
              up <- which(fromParent %in% updatedParents)
              for (j in up) fam$idHaps[ch,j] <- inheritedHap(fam$phasedAlls[j,,ch],
                                                             fam$phasedAlls[,,j],fam$IDS[j,,ch],fam$IDS[,,j])
            }
          }
        }
        return(fam)
      }  
      updateHap=function(haplo,subject,hapIndex){  # Update the subject's phasedAlls and IDS with 'haplo' in position 'hapIndex'
        haploIDS=ifelse(is.na(haplo),0,1)
        changes=which(haploIDS>subject$IDS[hapIndex,])
        if (length(changes)==0) return(NULL)
        if (length(changes)>0) for (j in changes){
          # If subject$IDS[hapIndex,j] is equal to 0 is because the allele in that position in the marker 
          # is NA, or because the pair of alleles in that marker are still without phasing
          subject$IDS[hapIndex,j]=1
          # If the other allele in the marker is already phased, is because the allele in this position 
          # is NA, and thus the value of haplo[j] can be directly assigned to it
          if (subject$IDS[3-hapIndex,j]==1) subject$phasedAlls[hapIndex,j]=haplo[j]
          else{
            wna=which(is.na(subject$phasedAlls[,j]))
            if (length(wna)==0){
              if (subject$phasedAlls[hapIndex,j]!=haplo[j]) subject$phasedAlls[,j]=subject$phasedAlls[2:1,j]
              subject$IDS[3-hapIndex,j]=1
            } else if(length(wna)==1){
              if (wna==hapIndex){
                if (subject$phasedAlls[3-hapIndex,j]==haplo[j]) subject$phasedAlls[,j]=subject$phasedAlls[2:1,j]
                else{
                  subject$phasedAlls[hapIndex,j]==haplo[j]
                  subject$IDS[3-hapIndex,j]=1
                } 
              } else{ # Case wna!=hapIndex
                if (subject$phasedAlls[hapIndex,j]!=haplo[j]){
                  subject$phasedAlls[3-hapIndex,j]=subject$phasedAlls[hapIndex,j]
                  subject$IDS[3-hapIndex,j]=1
                  subject$phasedAlls[hapIndex,j]=haplo[j]
                }                 
              }
            } else { # Case in which wna=2
              subject$phasedAlls[hapIndex,j]=haplo[j]
            }
          }          
        }
        return(subject)
      }
      writeFamilyHaps=function(fam){  # Writting of the phased haplotypes of a family
        writeHap=function(hap,ids){  # Writting of one haplotype
          hap[ids==0] <- NAsymbol
          return(paste(hap,collapse=alleSep))
        }
        haps <- matrix(,nrow=fam$famSize,ncol=2)
        for (i in 1:fam$famSize)
          for (j in 1:2)
            haps[i,j] <- writeHap(fam$phasedAlls[j,,i],fam$IDS[j,,i])
        haps <- matrix(haps,ncol=2)
        colnames(haps) <- c("hap1","hap2")
        return(haps)
      }
      
      # Function core
      nmiss0 <- length(which(is.na(fam$phasedAlls)))       # Counting of the number of missing values
      repeat{                                              # Find haplotypes in a recursive way 
        fam <- parentsHapsFromChildren(fam)                # Update of parental haplotypes from children information
        fam <- childrenHapsFromParents(fam)                # Updates children haplotypes from parents's information
        nmiss1 <- length(which(is.na(fam$phasedAlls)))     # Counting of the number of NAs after an iteration
        if (nmiss1==nmiss0) break else nmiss0 <- nmiss1    # If the number of NAs has been reduced
        # If TRUE updates nmiss0 to the new value and proceed to a new iteration, and if FALSE stop the loop
      }
      
      # Writting of the phased data 
      fam$haplotypes <- writeFamilyHaps(fam)
      fam$phasedMkrs <- as.familyMatrix(fam$phasedAlls)
      
      return(fam)
    }
    # Haplotype phasing for multiple families
    famsPhaser=function(fams){  
      
      #=== INTERNAL FUNCTIONS ===#
      as.familyArray=function(fam){
        outArray <- array(t(fam$imputedMkrs),dim=c(2,fam$nMkrs,fam$famSize), 
                          dimnames=list(haplo=c("hap1","hap2"),marker=paste("Mkr",1:fam$nMkrs,sep=""),
                                        subject=c("father","mother",paste("child",1:(fam$famSize-2),sep=""))))
        return(outArray)
      }
      as.familyMatrix=function(familyArray){
        outArray <- data.frame(t(array(familyArray,dim=c(2*dim(familyArray)[2],dim(familyArray)[3]))))
        nm <- ncol(outArray)/2
        names(outArray) <- paste("Mk",gl(nm,2),"_",gl(2,1,nm*2),sep="")
        return(outArray)
      }
      excludeCase5Mkrs=function(fam,case5){
        fam$nMkrs=fam$nMkrs-length(case5$mrk)
        fam$markers=fam$markers[,-case5$mrkCols]
        fam$imputedMkrs=fam$imputedMkrs[,-case5$mrkCols]
        fam$allelesNumber=fam$allelesNumber[-case5$mrk]
        fam$markerIncidences=fam$markerIncidences[-case5$mrk]
        fam$phasedAlls=fam$phasedAlls[,-case5$mrk,,drop=FALSE]
        return(fam)
      }
      includeCase5Mkrs=function(fam,fam0,case5){
        fam0$imputedMkrs[,-case5$mrkCols]=fam$imputedMkrs
        fam0$markerIncidences[-case5$mrk]=fam$markerIncidences
        fam0$phasedAlls[,-case5$mrk,]=fam$phasedAlls
        fam0$HMZ=matrix(0,nrow=fam0$famSize,ncol=fam0$nMkrs)
        fam0$HMZ[,-case5$mrk]=fam$HMZ
        fam0$HMZ[,case5$mrk]=t(apply(fam0$phasedAlls[,case5$mrk,,drop=FALSE],3,function(x) apply(x,2,function(a) as.integer(a[1]==a[2]))))
        fam0$IDS=array(0,dim=c(2,fam0$nMkrs,fam0$famSize))
        fam0$IDS[,-case5$mrk,]=fam$IDS
        fam0$idHaps=fam$idHaps
        fam0$haplotypes=writeFamilyHaps(fam0)
        fam0$phasedMkrs=as.familyMatrix(fam0$phasedAlls)
        nf0=names(fam0)
        nf=names(fam)
        nv=nf[which(!nf%in%nf0)]
        for (k in nv) fam0[[k]]=fam[[k]]
        return(fam0)
      }
      initializeIDS=function(fam){  
        ### Alleles arrangement 
        phasedAlls <- as.familyArray(fam)  
        ### HMZ Initialization
        HMZ <- t(apply(phasedAlls,3,function(x) apply(x,2,function(a) as.integer(a[1]==a[2]))))
        if (nrow(HMZ)==1) HMZ=t(HMZ)
        ### IDS Initialization
        IDS <- with(fam,array(,dim=c(2,fam$nMkrs,fam$famSize),
                              dimnames=list(haplo=c("hap1","hap2"),marker=paste("Mkr",1:fam$nMkrs,sep=""),
                                            subject=c("father","mother",paste("child",1:(fam$famSize-2),sep="")))))
        IDS[1,,] <- IDS[2,,] <- t(HMZ)
        IDS[is.na(IDS)] <- 0
        ## IDS in children
        for (i in 3:fam$famSize) {
          for (j in 1:fam$nMkrs){
            NAchild <- which(is.na(phasedAlls[,j,i]))
            NAfather <- which(is.na(phasedAlls[,j,1]))
            NAmother <- which(is.na(phasedAlls[,j,2]))
            # A child has only one NA value if the other value has been imputed from an homozygous parent
            if (length(NAchild)==1) {  
              hompar=which(HMZ[1:2,j]==1)
              chAl=phasedAlls[3-NAchild,j,i]
              proc=which(c(chAl %in% phasedAlls[,j,1],chAl %in% phasedAlls[,j,2]))
              if (length(proc)==1){
                if (NAchild==proc){
                  phasedAlls[,j,i]=phasedAlls[2:1,j,i]
                  NAchild=3-NAchild
                } 
                IDS[proc,j,i]=1
              } else if (length(proc)==2){
                if (length(hompar)==1){
                  if (NAchild==hompar) {
                    phasedAlls[,j,i]=phasedAlls[2:1,j,i]
                    NAchild=3-NAchild
                  }
                  IDS[hompar,j,i]=1
                }
              }
              if (length(hompar)==1){
                parAl=phasedAlls[1,j,hompar]
                if (chAl!=parAl){
                  phasedAlls[NAchild,j,i]=parAl
                  IDS[NAchild,j,i]=1
                }
              } else if (length(hompar)==2){
                phasedAlls[,j,i]=phasedAlls[1,j,1:2]
                IDS[,j,i]=1
              }
            } 
            else if (length(NAchild)==0){  # Child without missing alleles
              if (length(NAfather)+length(NAmother)==0&HMZ[i,j]==0){ # Parents without missing alleles
                fromWhichParent <- list(fromFather=which(phasedAlls[,j,i]%in%phasedAlls[,j,1]),
                                        fromMother=which(phasedAlls[,j,i]%in%phasedAlls[,j,2]))
                ia <- lapply(fromWhichParent,length)
                if (ia[[1]]+ia[[2]]<4){
                  ia1 <- which(ia==1)[1]   
                  if (fromWhichParent[[ia1]]!=ia1) phasedAlls[,j,i] <- phasedAlls[2:1,j,i]
                  IDS[,j,i] <- 1        
                }
              } else if (length(NAfather)>0&length(NAmother)==0){ # Missing alleles only in father
                fromFather <- which(!phasedAlls[,j,i]%in%phasedAlls[,j,2])  
                if (length(fromFather)==1){ 
                  if (fromFather==2) phasedAlls[,j,i] <- phasedAlls[2:1,j,i]
                  IDS[,j,i] <- 1
                }
              } else if (length(NAfather)==0&length(NAmother)>0){ # Missing alleles only in mother
                fromMother <- which(!phasedAlls[,j,i]%in%phasedAlls[,j,1])
                if (length(fromMother)==1){ 
                  if (fromMother==1) phasedAlls[,j,i] <- phasedAlls[2:1,j,i]
                  IDS[,j,i] <- 1
                }
              }
            }
          }
        }
        ## IDS in parents
        for (i in 1:2){    # IDS in parents is initialized from the child with more IDS==1
          if (fam$famSize==3) keyCh<-1
          else{
            nrIDS1=if (fam$nMkrs==1) IDS[i,,-(1:2)] else colSums(IDS[i,,-(1:2)],na.rm=TRUE)  
            keyCh<-which.max(nrIDS1)   
          }
          keyMks=which(IDS[i,,keyCh+2]==1)
          for (j in keyMks){
            pap=which(match(phasedAlls[,j,i],phasedAlls[i,j,keyCh+2])==1)
            if (!is.na(pap[1])) if (pap[1]==2) phasedAlls[,j,i]=phasedAlls[2:1,j,i]
            IDS[which(!is.na(phasedAlls[,j,i])),j,i]=1
          }
        }
        idHaps <- array(0,dim=c(fam$famSize,2)) # Which parents allele has been inherited by child
        for (i in 3:fam$famSize)
          for (j in 1:2)
            idHaps[i,j] <- inheritedHap(phasedAlls[j,,i],phasedAlls[,,j],IDS[j,,i],IDS[,,j])
        # Adding of new elements to fam
        fam$phasedAlls <- phasedAlls
        fam$HMZ <- HMZ
        fam$IDS <- IDS
        fam$idHaps <- idHaps
        return(fam)
      }
      inheritedHap=function(hap,parentHaps,hapIDS,pHapsIDS){   # Generation of idHaps matrix
        # This function determines which one of parentHaps is (or can be) equal to 'hap'
        # it returns 0 if both parentHaps are compatible with 'hap' and 3 if there is no compatibility
        if(is.null(dim(parentHaps))){ 
          dim(hap) <- dim(hapIDS) <- c(1,1)
          dim(parentHaps) <- dim(pHapsIDS) <- c(2,1)
        } 
        canBeSameHap=vector("logical",2)
        for (i in 1:2){
          pos=which(hapIDS+pHapsIDS[i,]==2)   
          if (length(pos)==0) canBeSameHap[i]=TRUE
          else canBeSameHap[i]=all(hap[pos]==parentHaps[i,pos])  
          if (length(pos)==ncol(pHapsIDS)&canBeSameHap[i]) return(i)
        }
        inhrtdHap=which(canBeSameHap)
        if (length(inhrtdHap)==2) inhrtdHap=0
        else if (length(inhrtdHap)==0) inhrtdHap=3  # Irregular inheritance
        return(inhrtdHap)
      }
      phaseHapsFromTwoChildren=function(fam){
        # Internal function to identify in which case are the parent-offspring haplotypes
        identifyCase=function(fam,phasedMarkers,childrenPair){
          writeHaps=function(hap,ids){  # Writting of one haplotype
            hap[ids==0] <- NAsymbol
            return(paste(hap,collapse=alleSep))
          }
          group=fam$phasedAlls[,phasedMarkers,c(1,2,childrenPair+2),drop=FALSE]
          haps=matrix(NA,nrow=4,ncol=2)
          for (i in 1:4) for (j in 1:2) haps[i,j]=writeHaps(group[j,,i],1)
          haps=matrix(haps,ncol=2)
          colnames(haps)=c("hap1","hap2")
          commChHap=intersect(haps[3,],haps[4,])  # Do the children have a haplotype in common?
          ncch=length(commChHap)  # How many different haplotypes have these two children?
          nChHaps=length(unique(as.vector(haps[3:4,])))
          commParHap=intersect(haps[1,],haps[2,]) # Do the parents have a haplotype in common?
          ncph=length(commParHap)
          nParHaps=length(unique(as.vector(haps[1:2,]))) # How many different haplotypes have the parents?
          case=0
          if (ncph==0){
            if (nChHaps==2){
              if (haps[1,1]==haps[1,2]) case=3
              else if(haps[2,1]==haps[2,2]) case=4
            } else if (nChHaps==3){
              if (nParHaps==4){
                if (commChHap%in%haps[1,]) case=1
                else case=2
              }
            }
          } else if (ncph==1){
            uch=unique(as.vector(haps[3:4,]))
            fcomm=length(intersect(uch,haps[1,]))
            mcomm=length(intersect(uch,haps[2,]))
            if (fcomm==2&mcomm==1&haps[2,1]!=haps[2,2]) case=6
            else if (mcomm==2&fcomm==1&haps[1,1]!=haps[1,2]) case=5
          }
          return(case)
        }
        
        # Identification of such markers containing IDS=1 and missing parental values 
        phasedParentsMarkers=which(apply(fam$IDS[,,1:2],2,sum)==4)
        if (is.null(phasedParentsMarkers)) return(fam)
        
        # Markers with two heterozygous children (at least) and both parents with missing data
        fam$HMZ=t(apply(fam$phasedAlls,3,function(x) apply(x,2,function(a) as.integer(a[1]==a[2]))))
        twoHetChMkr=which(apply(fam$HMZ[-(1:2),],2,function(x) length(which(x==0)))>1)
        naParents=which(apply(fam$HMZ[1:2,],2,function(x) if (any(is.na(x))) TRUE else FALSE))
        mk2ch=intersect(twoHetChMkr,naParents)  # Markers where the phaseHapsFromTwoChildren can be applied
        
        # Search for those sibs pairs containing one common allele and two alleles from each other
        imputedMarkers=NULL
        for (mk in mk2ch){
          chPairs=combn(which(fam$HMZ[-(1:2),mk]==0),2)
          selec=which(apply(chPairs,2,function(cp){
            alleles=fam$phasedAlls[,mk,cp+2,drop=FALSE]
            (length(unique(as.vector(alleles)))==3)&(length(intersect(alleles[,,1],alleles[,,2]))==1)
          }))
          if (length(selec)>0){
            chp=chPairs[,selec,drop=FALSE] # Heterozygous children pairs in this marker with a common allele
            for (s in 1:length(selec)){
              phasedMarkers=phasedParentsMarkers[which(apply(fam$IDS[,phasedParentsMarkers,chp[,s]],2,sum)==4)]
              case=identifyCase(fam,phasedMarkers,chp[,s])
              if (case>0){
                chAlle=fam$phasedAlls[,mk,chp[,s]+2,drop=FALSE]
                commAllele=intersect(chAlle[,,1],chAlle[,,2])
                if (case%in%c(1,4,5)) ref=2 else ref=1
                for (i in 1:2) if (chAlle[ref,,i]==commAllele) chAlle[,,i]=chAlle[2:1,,i]
                fam$phasedAlls[,mk,chp[,s]+2]=chAlle
                fam$IDS[,mk,chp[,s]+2]=1
                fam <- famPhaser(fam) 
                imputedMarkers=c(imputedMarkers,mk)
                break
              }
            }
          }
        }
        return(fam)
      }
      solveHapsWithMissingParents=function(fam){  
        phaseFamilyFromChildren=function(fam){  
          # Internal functions
          canBeEqual=function(hapPair1,hapPair2){
            # Compares two haplotype pairs
            for (j in 1:ncol(hapPair1)){
              equalMrk=(all(hapPair1[,j]==hapPair2[,j],na.rm=TRUE))|
                (all(hapPair1[,j]==hapPair2[,j][2:1],na.rm=TRUE))
              if (!equalMrk) return(FALSE)
            }
            return(TRUE)
          }  
          compareNA=function(v1,v2){
            same <- (v1 == v2) | (is.na(v1) & is.na(v2))
            same[is.na(same)] <- FALSE
            return(same)
          }
          fillChildrenNAs=function(children,parentsHaps){
            identify <- function(child, possibleSibs){
              nmk=dim(child)[2]
              sibComp=logical(4)
              for (k in 1:4){
                mkComp=logical(nmk)
                for (j in 1:nmk){
                  mkComp[j]=all(child[,j]==possibleSibs[,j,k],na.rm=TRUE)|all(child[2:1,j]==possibleSibs[,j,k],na.rm=TRUE)
                }
                sibComp[k]=all(mkComp)
              }
              wsbc=which(sibComp)
              if (length(wsbc)==0) return(list(child=NULL,idsChild=NULL))  ## Si ocurre esto el hijo es imposible y el programa debe parar
              else if (length(wsbc)==1){
                child=possibleSibs[,,wsbc]
                idsChild=array(1,dim=dim(child))
              }
              else{
                coincidences=possibleSibs[,,wsbc[1]]
                for (k in 2:length(wsbc)){
                  idtc=(coincidences==possibleSibs[,,wsbc[k]])
                  coincidences[!idtc]=NA
                }
                idsChild=array(0,dim=dim(child))
                idsChild[!is.na(coincidences)]=1
                for (j in 1:nmk){
                  addAl=child[which(!child[,j]%in%coincidences[,j]),j]
                  addAl=addAl[!is.na(addAl)]
                  if (length(addAl)>0){
                    coincGap=which(is.na(coincidences[,j]))
                    if (length(coincGap)<length(addAl)) return(list(child=NULL,idsChild=NULL))
                    else coincidences[coincGap,j]=addAl
                  }
                }
                child=coincidences
              }
              return(list(child=child,idsChild=idsChild))  
            }     
            possibleSibs=abind(rbind(parentsHaps[1,,1],parentsHaps[1,,2]),
                               rbind(parentsHaps[1,,1],parentsHaps[2,,2]),
                               rbind(parentsHaps[2,,1],parentsHaps[1,,2]),
                               rbind(parentsHaps[2,,1],parentsHaps[2,,2]),along=3)
            ids=array(0,dim=dim(children))
            for (k in 1:dim(children)[3]){
              ich=identify(children[,,k],possibleSibs)
              if (is.null(ich$child)) return(list(children=NULL,ids=NULL))
              children[,,k]=ich$child
              ids[,,k]=ich$idsChild
            }
            return(list(children=children,ids=ids))
          }
          find2Haps=function(mat){
            canBeEqual=function(v1,v2) all(v1==v2,na.rm=TRUE)
            nr=nrow(mat)
            differents=NULL
            for (i in 1:(nr-1)){
              for(j in 2:nr)
                if (!canBeEqual(mat[i,],mat[j,])){
                  differents=c(i,j)
                  break
                }
              if (!is.null(differents)) break
            }
            if (is.null(differents)){
              nna=apply(mat,1,function(x) length(which(is.na(x))))
              return(rbind(mat[which.min(nna),],rep(NA,ncol(mat))))
            } else{
              hap=rbind(mat[i,],mat[j,])
              hta=(1:nr)[-c(i,j)]
              repeat{
                natot0=length(which(is.na(hap)))
                for (k in hta){
                  if(canBeEqual(mat[k,],hap[1,])&!canBeEqual(mat[k,],hap[2,])){
                    h1na=intersect(which(is.na(hap[1,])),which(!is.na(mat[k,])))
                    hap[1,h1na]=mat[k,h1na]
                    hta=hta[-which(hta==k)]
                  } 
                  else if (!canBeEqual(mat[k,],hap[1,])&canBeEqual(mat[k,],hap[2,])){
                    h2na=intersect(which(is.na(hap[2,])),which(!is.na(mat[k,])))
                    hap[2,h2na]=mat[k,h2na]
                    hta=hta[-which(hta==k)]
                  }
                }
                natot1=length(which(is.na(hap)))
                if (natot0==natot1) break
              }
              return(hap)
            }
          }
          findHapsfrom3or4Children=function(children,parentsInfo=NULL){
            ### Internal Functions
            commonHap=function(twoChildren){
              # Finds the possible common haplotypes between two children (stored in the array 'twoChildren'). 
              commonAllelles <- vector("list")  
              nMkrs <- dim(twoChildren)[2]
              for (j in 1:nMkrs) commonAllelles[[j]] <- intersect(twoChildren[,j,1],twoChildren[,j,2])          
              commonHap <- expand.grid(commonAllelles,KEEP.OUT.ATTRS=FALSE)
              names(commonHap) <- paste("Mkr",1:ncol(commonHap),sep="")
              return(commonHap)
            }
            complementHap=function(hap,child){
              # Determines the complementary haplotype for a given one. 
              # This function supose that "hap" is placed in "child", and it finds the other haplotype.
              comp <- child[1,]
              flip <- (comp==hap)
              comp[flip] <- child[2,flip]
              return(comp)
            }
            findParHaps=function(chldHaps){  
              # Determination of the parental haplotypes from the haplotypes found in children
              foundHaps=t(apply(chldHaps,3,function(x) apply(x,1,paste,collapse="")))
              nHaps=length(unique(as.vector(foundHaps))) # Number of unique haplotypes
              hmzCh=sum(foundHaps[,1]==foundHaps[,2])    # Number of homozygous children
              if(nHaps==3&hmzCh==0)  
                return(list(phasedParents=NULL, phasedChildren=chldHaps, 
                            hapIncidences="Multiple compatible parental haplotypes",alignedParents=FALSE))
              hapIncidences=""
              parentOrdered=all(apply(foundHaps,2,function(x) length(unique(x)))<3)
              if (!parentOrdered){
                i=1
                repeat{
                  foundHaps[i,]=foundHaps[i,2:1]
                  parentOrdered=all(apply(foundHaps,2,function(x) length(unique(x)))<3)
                  if (!parentOrdered){
                    foundHaps[i,]=foundHaps[i,2:1]
                    i=i+1
                    if (i==4) {
                      parents=array(NA,dim=c(2,nMrk,2))
                      hapIncidences="Irregular inheritance detected"
                      break
                    }
                  } 
                  else{
                    chldHaps[,,i]=chldHaps[2:1,,i]
                    break
                  }
                }
              }
              if (hapIncidences==""){
                parents=abind(t(unique(chldHaps[1,,],MARGIN=2)), t(unique(chldHaps[2,,],MARGIN=2)),along=3)
                parents=abind(sortByRows(parents[,,1]),sortByRows(parents[,,2]),along=3)
                if (is.null(parentsInfo)) alignedParents=FALSE else {
                  eqPar=matrix(NA,2,2)
                  for (i in 1:2)
                    for (j in 1:2)
                      eqPar[i,j]=canBeEqual(parents[,,i],parentsInfo[,,j])
                  if (any(rowSums(eqPar)==0)|any(colSums(eqPar)==0)) {
                    alignedParents=FALSE
                    parents[,,1]<-parents[,,2]<-NA
                    hapIncidences="Parental information is not compatible with haplotypes found in children"
                  } else{
                    weq=which(!eqPar,arr.ind=TRUE)
                    if (length(weq)==0){
                      alignedParents=FALSE
                    } else{
                      alignedParents=TRUE
                      if (weq[1,1]==weq[1,2]){
                        parents=parents[,,2:1]
                        chldHaps=chldHaps[2:1,,]
                      } 
                    }
                  }
                }      
              }    
              return(list(phasedParents=parents, phasedChildren=chldHaps, 
                          hapIncidences=hapIncidences, alignedParents=alignedParents))
            }  
            hapInChild=function(hap,child){
              # Determines if the haplotype 'hap' is located in a given child. In such a case returns the children haplotypes
              # including 'hap' and its complementary haplotype.
              # If a haplotype is not found, the function returns "NULL".
              if (is.null(dim(child))) child=matrix(child,ncol=1)
              isin=(t(child)==hap)
              if (any(rowSums(isin)==0)) return(NULL)
              s=which(!isin[,1])
              child[,s]=child[2:1,s]
              return(matrix(child,2))
            }
            identifyFourthChild=function(famHaps,fourthChild){
              # Identification of haplotypes in the fourth child
              c4possibleHaps=list(t(famHaps$phasedParents[1,,1:2]),t(famHaps$phasedParents[2,,1:2]),
                                  rbind(famHaps$phasedParents[1,,1],famHaps$phasedParents[2,,2]),
                                  rbind(famHaps$phasedParents[2,,1],famHaps$phasedParents[1,,2]))
              c4possibleHapsSorted=lapply(c4possibleHaps,function(ph) apply(ph,2,sort))
              fourthChild=apply(fourthChild,2,sort)
              c4Haps=which(sapply(c4possibleHapsSorted,function(x) all(x==fourthChild | x==fourthChild[2:1,])))
              if (length(c4Haps)>0) famHaps$phasedChildren=abind(famHaps$phasedChildren, c4possibleHaps[[c4Haps]],along=3) 
              else{
                famHaps$phasedChildren=abind(famHaps$phasedChildren, fourthChild,along=3) 
                famHaps$hapIncidences="Haplotypes in one child are not compatible with the haplotypes found in the rest of the offspring"
              } 
              return(famHaps)
            }
            sortByRows=function(arr) {       # Sorting of a matrix row by row
              arr[do.call(order, lapply(1:NCOL(arr), function(i) arr[, i])), ]
            }
            
            ### Initialization of variables
            nMrk <- dim(children)[2]
            nchld <- dim(children)[3]
            if (nchld==4) {
              fourthChild <- children[,,4]
              children <- children[,,1:3]
            }
            if (nchld<3){
              return(list(phasedParents=array(NA,dim=c(2,nMrk,2)),phasedChildren=children,alignedParents=FALSE,
                          hapIncidences="Less than two children detected",
                          childrenIDS=array(0,dim=c(2,nMrk,nchld)), parentsIDS=array(0,dim=c(2,nMrk,2))))
            }
            
            ### Haplotype finding from three unique children without missing values nor recombination.
            commonHaps <- vector("list")
            for (i in 1:3) commonHaps[[i]] <- commonHap(children[,,-i])
            nch <- sapply(commonHaps,nrow)
            refH=which(nch>0)
            chHaps=list()
            for (ref in refH) {    
              chld=(1:3)[-ref]
              # The complemmentary haplotypes are added in each child
              possibleHaps <- array(rbind(t(commonHaps[[ref]]), apply(commonHaps[[ref]],1,complementHap,children[,,chld[1]]),
                                          apply(commonHaps[[ref]],1,complementHap,children[,,chld[2]])),dim=c(nMrk,3,nch[ref]))    
              for (ih in 1:nch[ref]){
                identifiedHaps <- t(possibleHaps[,,ih])
                uHaps <- unique(identifiedHaps)
                for (j in 1:nrow(uHaps)) {
                  hapsChRef <- hapInChild(uHaps[j,],children[,,ref])
                  if (!is.null(hapsChRef)){
                    nHaps <- nrow(unique(rbind(uHaps,hapsChRef)))
                    if (nHaps<=4){              
                      chh <- rbind(identifiedHaps[c(1,2,1,3),],hapsChRef)
                      chh <- chh[order(2*c(sort(rep((1:3)[-ref],2)),ref,ref)-c(1,0)),]            
                      nh=nrow(unique(chh)) # Number of different haplotypes in the three children
                      ndh=nrow(unique(chh[duplicated(chh),])) # Num. of haplotypes which are repeated
                      if (!(nh==4&ndh==1)) 
                        chHaps[[length(chHaps)+1]]<-array(cbind(sortByRows(chh[1:2,]),sortByRows(chh[3:4,]),
                                                                sortByRows(chh [5:6,])),dim=c(2,nMrk,3))           
                    }
                  }
                }
              }
            }
            chHaps=unique(chHaps)
            
            ### Identification and reconstruction of the posible parental haplotypes
            famHaps=lapply(chHaps,findParHaps)
            noIncidences=which(sapply(famHaps,function(x) nchar(x$hapIncidences))==0)
            # If there is only one posible haplotype structure in children
            famHaps <- if (length(noIncidences)==0) famHaps[1] else famHaps[noIncidences]
            # If there is a fourth child (diferent from the rest of the offspring)
            if (nchld>3){
              famHaps=lapply(famHaps,identifyFourthChild,fourthChild)
              noIncidences=which(sapply(famHaps,function(x) nchar(x$hapIncidences))==0)
              famHaps= if (length(noIncidences)==0) famHaps[1] else famHaps[noIncidences]
            }
            
            ### IDS calculation
            if(length(famHaps)>1) return(NULL) else famHaps=famHaps[[1]]
            ids <- if (famHaps$hapIncidences=="") 1 else 0
            famHaps$childrenIDS <- array(ids,dim=c(2,nMrk,nchld))
            famHaps$parentsIDS <- array(ids,dim=c(2,nMrk,2))
            
            ### Function output
            return(famHaps)
          }
          isCase5=function(mrk){
            # Indicates if a marker is in case 5
            globAlNr=length(unique(as.vector(mrk[!is.na(mrk)])))  # Global number of alleles in marker
            indivAlNr=apply(mrk,3,function(x) if(any(is.na(x))) 2 else length(unique(x)))
            return((globAlNr==2) & all(indivAlNr==2))
          }
          possibleSibs=function(parentsHaps){
            # Possible siblings
            return(abind(rbind(parentsHaps[1,,1],parentsHaps[1,,2]),
                         rbind(parentsHaps[1,,1],parentsHaps[2,,2]),
                         rbind(parentsHaps[2,,1],parentsHaps[1,,2]),
                         rbind(parentsHaps[2,,1],parentsHaps[2,,2]),along=3))
          }
          selChMk=function(ndx,children){
            alleMCh=alleMat[ndx,]
            noNAMks=which(colSums(alleMCh)==0)
            if (length(noNAMks)>0){
              c5=which(sapply(noNAMks,function(k) isCase5(children[,k,ndx,drop=FALSE])))
              if (length(c5)>0) noNAMks=noNAMks[-c5]
            }
            if (length(noNAMks)>1){
              sch=array(apply(children[,noNAMks,ndx,drop=FALSE],3, function(M) apply(M,2,sort)),dim=c(2,length(noNAMks),3))
              idupCh=duplicated(sch,MARGIN=3)
              uch=ndx[!idupCh]  # Unique complete Children
              if (length(uch)>=3) return(list(children=uch,markers=noNAMks))
            }
          }
          tryParentAlign=function(fam,selectedMarkers,partialFam,children,idsCh){
            f1=f2=fam
            alignmentMks=fam$alignmentMks
            f1$phasedAlls[,selectedMarkers,]=abind(partialFam$phasedParents,children,along=3)
            f1$IDS[,selectedMarkers,]=abind(partialFam$parentsIDS,idsCh,along=3)
            f2$phasedAlls[,selectedMarkers,]=abind(partialFam$phasedParents[,,2:1],children[2:1,,],along=3)
            f2$IDS[,selectedMarkers,]=abind(partialFam$parentsIDS[,,2:1],idsCh[2:1,,],along=3)
            f1$phasedAlls[,alignmentMks$alinMk,]=alignmentMks$phasedAlls
            f1$IDS[,alignmentMks$alinMk,]=alignmentMks$IDS
            f2$phasedAlls[,alignmentMks$alinMk,]=alignmentMks$phasedAlls
            f2$IDS[,alignmentMks$alinMk,]=alignmentMks$IDS
            f1=updateParentsIDS(f1)
            f2=updateParentsIDS(f2)
            compat1=!any(f1$idHaps==3)
            compat2=!any(f2$idHaps==3)
            if (compat1&!compat2){
              fam=famPhaser(f1)
            } else if (compat2&!compat1){
              fam=famPhaser(f2)
            } else if (compat1&compat2){
              for (i in alignmentMks$alinMk){
                cal=intersect(fam$phasedAlls[,i,1],fam$phasedAlls[,i,2])
                if (length(cal)==1) fam$phasedAlls[,i,1]<-fam$phasedAlls[,i,2]<-c(cal,NA)
                else if (length(cal)==0) fam$phasedAlls[,i,1]<-fam$phasedAlls[,i,2]<-c(NA,NA)
                fam$IDS[,i,]=0
              }
              fam$alignedParents=FALSE
            } else {
              fam$hapIncidences <- "Irregular inheritance detected"
            }
            return(fam)
          }
          updateParentsIDS=function(fam){
            # IDS in parents is initialized from the child with more IDS==1
            for (i in 1:2){
              fam$IDS[,,i]=0
              nrIDS1=if (fam$nMkrs==1) fam$IDS[i,,-(1:2)] else colSums(fam$IDS[i,,-(1:2)],na.rm=TRUE)  
              keyCh<-which.max(nrIDS1)  
              keyMks=which(fam$IDS[i,,keyCh+2]==1)
              for (j in keyMks){
                chAl=fam$phasedAlls[i,j,keyCh+2]
                newa=which(!chAl%in%fam$phasedAlls[,j,i])
                if(length(newa)>0){
                  pos=which(is.na(fam$phasedAlls[,j,i]))[1:length(newa)]
                  fam$phasedAlls[pos,j,i]=chAl[newa]
                }
                pap=which(match(fam$phasedAlls[,j,i],fam$phasedAlls[i,j,keyCh+2])==1)
                if (!is.na(pap[1])) if (pap[1]==2) fam$phasedAlls[,j,i]=fam$phasedAlls[2:1,j,i]
                fam$IDS[which(!is.na(fam$phasedAlls[,j,i])),j,i]=1
              }
            }
            # Which parents allele has been inherited by child 
            fam$idHaps<-array(0,dim=c(fam$famSize,2))
            for (i in 3:fam$famSize)
              for (j in 1:2)
                fam$idHaps[i,j]=inheritedHap(fam$phasedAlls[j,,i],fam$phasedAlls[,,j],fam$IDS[j,,i],fam$IDS[,,j])
            return(fam)
          }
          
          if (!exists("alignedParents",fam)){
            if (all(is.na(fam$phasedAlls[,,1:2]))) fam$alignedParents=FALSE
            else{
              idsM1=ifelse(apply(fam$IDS[,,-(1:2)],2,sum)>0,TRUE,FALSE)  # Markers with ids=1 in some sib
              c5=apply(fam$phasedAlls[,,-(1:2)],2,function(x){
                x=array(x,dim=c(2,1,fam$famSize-2))
                isCase5(x)}
              )
              # Determination of which haps come from father and which from mother (which markers can be used for alignment of parents)
              alinMk=which(!apply(fam$phasedAlls[,,1:2],2,function(M) all(compareNA(M[,1],M[,2])))&idsM1&!c5)  
              fam$alignedParents=if (length(alinMk)>0) TRUE else FALSE
              fam$alignmentMks=if (fam$alignedParents) list(alinMk=alinMk,phasedAlls=fam$phasedAlls[,alinMk,,drop=FALSE],
                                                            IDS=fam$IDS[,alinMk,,drop=FALSE]) else NULL
            }
          }
          nChildren=fam$famSize-2
          alleMat=t(apply(fam$phasedAlls[,,-c(1,2)],3,function(x) apply(x,2,function(x) any(is.na(x))))+0)
          # The sum by rows indicates the number of missing markers on each child      
          
          combChildren=combn(nChildren,3) # Search for marker combinations for three full children (at least)
          scm=apply(combChildren,2,selChMk,fam$phasedAlls[,,-c(1,2)])  
          if (!is.null(scm)) scm=scm[-which(sapply(scm,is.null))]
          nmk=sapply(scm,function(cm){length(cm$imputedMkrs)})
          if (all(nmk<=1)) return(fam)
          
          sl=order(nmk,decreasing = TRUE)
          sl=sl[nmk[sl]>1]
          selSet=NULL
          selMks=NULL
          for (k in sl){
            nwm=which(!scm[[k]]$imputedMkrs%in%selMks)
            if (length(nwm)>0){
              selSet=c(selSet,k)
              selMks=c(selMks,scm[[k]]$imputedMkrs[nwm])
            }
          }
          for (s in selSet) if(any(is.na(fam$phasedAlls[,scm[s][[1]]$imputedMkrs,]))){
            selected=scm[s][[1]]
            idsCh=fam$IDS[,selected$imputedMkrs,-(1:2)]
            children=fam$phasedAlls[,selected$imputedMkrs,-c(1,2)]
            partialFam=findHapsfrom3or4Children(fam$phasedAlls[,selected$imputedMkrs,selected$children+2],
                                                parentsInfo=fam$phasedAlls[,selected$imputedMkrs,1:2])
            if (!is.null(partialFam$phasedParents)){
              children[,,selected$children]=partialFam$phasedChildren
              idsCh[,,selected$children]=partialFam$childrenIDS
              if(length(selected$children)<nChildren){
                unidentifiedChildren=children[,,-selected$children,drop=FALSE]
                fillCh=fillChildrenNAs(unidentifiedChildren,partialFam$phasedParents)
                children[,,-selected$children]=fillCh$children
                idsCh[,,-selected$children]=fillCh$ids
              }
              
              if (partialFam$alignedParents){
                fam$phasedAlls[,selected$imputedMkrs,1:2]=partialFam$phasedParents
                fam$phasedAlls[,selected$imputedMkrs,-(1:2)]=children
                fam$IDS[,selected$imputedMkrs,-(1:2)]=idsCh
                fam=updateParentsIDS(fam)
              } else{
                if (!is.null(fam$alignmentMks)) 
                  fam=tryParentAlign(fam,selected$imputedMkrs,partialFam,children,idsCh)
                if (is.null(fam$alignmentMks)|!fam$alignedParents){
                  oldCh=fam$phasedAlls[,,-(1:2)]
                  oldCh[fam$IDS[,,-(1:2)]==0]=NA
                  newCh=children
                  newCh[idsCh==0]=NA
                  newCh1<-newCh2<-oldCh
                  newCh1[,selected$imputedMkrs,]=newCh
                  newCh2[,selected$imputedMkrs,]=newCh[2:1,,]
                  ids1old=t(apply(fam$IDS[,selected$imputedMkrs,-(1:2)],3,function(M) ifelse(colSums(M)>0,1,0)))
                  ids1new=t(apply(idsCh,3,function(M) ifelse(colSums(M)>0,1,0)))
                  ids1comm=ifelse(ids1old+ids1new==2,1,0)
                  ids1comm=which(ids1comm==1,arr.ind=TRUE)
                  het=apply(ids1comm,1,function(v) children[1,v[2],v[1]]!=children[2,v[2],v[1]])
                  hetMk=which(het)
                  if (length(hetMk)>0){
                    hcm=ids1comm[hetMk[1],]
                    child0=fam$phasedAlls[,selected$imputedMkrs,hcm[1]+2]
                    gmk0=child0[,hcm[2]]
                    gmk=children[,hcm[2],hcm[1]]
                    if (all(gmk==gmk0)) {
                      fam$phasedAlls[,selected$imputedMkrs,-(1:2)]=children
                      fam$IDS[,selected$imputedMkrs,-(1:2)]=idsCh
                      fam$phasedAlls[,selected$imputedMkrs,(1:2)]=partialFam$phasedParents
                    } else{
                      fam$phasedAlls[,selected$imputedMkrs,-(1:2)]=children[2:1,,]
                      fam$IDS[,selected$imputedMkrs,-(1:2)]=idsCh[2:1,,]
                      fam$phasedAlls[,selected$imputedMkrs,(1:2)]=partialFam$phasedParents[,,2:1]
                    }
                    fam=updateParentsIDS(fam)
                  } 
                  else if (all(compareNA(newCh1,newCh2[2:1,,]))){ 
                    fam$phasedAlls[,selected$imputedMkrs,1:2]=partialFam$phasedParents
                    fam$phasedAlls[,selected$imputedMkrs,-(1:2)]=children
                    fam$IDS[,selected$imputedMkrs,-(1:2)]=idsCh
                    fam$alignedParents=FALSE  
                    # This means that parents will be assigned haplotypes but it is not known which 
                    # haps really correspond to father and which correspond to mother.
                    fam=updateParentsIDS(fam)
                  } 
                  else{
                    cpr=compareNA(newCh1,newCh2[2:1,,])
                    fam$phasedAlls[,selected$imputedMkrs,-(1:2)]=children
                    fam$IDS[,selected$imputedMkrs,-(1:2)]=idsCh
                    fam$IDS[,,-(1:2)]=(cpr+0)*fam$IDS[,,-(1:2)]
                    fch=fam$phasedAlls[,,-(1:2)]
                    ids0=fam$IDS[,,-(1:2)]==0
                    fch[ids0]=NA
                    posParHaps=list(unique(t(fch[1,,])),unique(t(fch[2,,]))) # Possible parental haplotypes
                    phaps=lapply(posParHaps,find2Haps)
                    for (i in 1:2) fam$phasedAlls[,,i]=phaps[[i]]
                    fam=updateParentsIDS(fam)
                  }
                }
              }
              fam=famPhaser(fam)
            }
          }
          if (!fam$alignedParents){
            fam$unAlignedParentsHaps=fam$phasedAlls[,,1:2]
            dimnames(fam$unAlignedParentsHaps)[[3]]=c("Parent1","Parent2")
            cnf=list(fam$phasedAlls[,,1]==fam$phasedAlls[,,2],
                     fam$phasedAlls[,,1]==fam$phasedAlls[2:1,,2])
            fid=fam$IDS[,,1]*fam$IDS[,,2]
            cnf=lapply(cnf,function(id) {id[is.na(id)]=FALSE; id*fid})
            idsC=cnf[[which.max(sapply(cnf,sum))]]
            phParAlls=fam$phasedAlls[,,1]
            phParAlls[!idsC]=NA
            fam$phasedAlls[,,1:2]=phParAlls
            fam$IDS[,,1:2]=idsC
            fam$haplotypes=writeFamilyHaps(fam)
          }
          return(fam)
        }
        na0=sum(is.na(fam$phasedAlls))
        repeat{
          fam= tryCatch(phaseFamilyFromChildren(fam),error=function(e) return(fam))
          na1=sum(is.na(fam$phasedAlls))
          if (na1==na0) break
          na0=na1
        }
        return(fam)
      }
      whichCase5Mkrs=function(fam){ # Identification of those markers in case 5 from fam$imputedMkrs matrix
        isCase5=function(mrkPos){
          mrk=fam$imputedMkrs[,mrkPos]
          globAlNr=length(unique(as.vector(mrk[!is.na(mrk)])))  # Global number of alleles in marker
          indivAlNr=apply(mrk,1,function(x) if(any(is.na(x))) 2 else length(unique(x)))
          if ((globAlNr==2) & all(indivAlNr==2)) return(TRUE) else return(FALSE)
        }
        mrksPos=matrix(1:ncol(fam$imputedMkrs),nrow=2)
        case5Mks=which(apply(mrksPos,2,isCase5))
        return(list(mrk=case5Mks,mrkCols=mrksPos[,case5Mks]))
      }
      writeFamilyHaps=function(fam){  # Writting of the phased haplotypes of a family
        writeHap=function(hap,ids){  # Writting of one haplotype
          hap[ids==0] <- NAsymbol
          return(paste(hap,collapse=alleSep))
        }
        haps <- matrix(,nrow=fam$famSize,ncol=2)
        for (i in 1:fam$famSize)
          for (j in 1:2)
            haps[i,j] <- writeHap(fam$phasedAlls[j,,i],fam$IDS[j,,i])
        haps <- matrix(haps,ncol=2)
        colnames(haps) <- c("hap1","hap2")
        return(haps)
      }
      
      #=== INITIALIZATION OF FAMS LIST ELEMENTS ===#
      idFams <- unique(fams$imputedMkrs[,1])                       # Family identifier
      nMkrs <- (ncol(fams$imputedMkrs)-6)/2                        # Number of markers 
      fams$phasedMkrs <- as.data.frame(setNames(replicate(nMkrs*2,character(0), simplify=F), 
                                                names(fams$imputedMkrs)[-(1:6)]))  # Phased Markers
      fams$haplotypes <- NULL                                      # Phased Haplotypes
      fams$Vdata <- NULL                                           # Attached information
      fams$hapIncidences <- NULL                                   # Initialization of the list of families
      
      #=== FAMILIES SCAN ===#
      for (f in idFams) {                                          # Scanning of each family
        family <- fams$imputedMkrs[fams$imputedMkrs[,1]==f,]       # Family data
        family <- family[order(family[,3]),]                       # Family sorting according to the third column (paternal ID)
        fam <- list(imputedMkrs=family[,-(1:6)],
                    famSize=nrow(family),nMkrs=nMkrs)              # 
        fam$phasedAlls <- as.familyArray(fam)                      # 
        case5Mks <- whichCase5Mkrs(fam)                            # 
        if (length(case5Mks$mrk)>0){                               #
          fam0 <- fam                                              #
          fam <- excludeCase5Mkrs(fam,case5Mks)                    #
        }
        if (nMkrs-length(case5Mks$mrk)<2) {                                      # It has to be two markers (at least) to start the phasing process
          names(fam0$imputedMkrs) <- names(fams$phasedMkrs)                      # 
          fams$phasedMkrs <- rbind(fams$phasedMkrs,fam0$imputedMkrs)             # 
          fam <- initializeIDS(fam0)                                             # IDS matrix creation
          fams$haplotypes <- rbind(fams$haplotypes,writeFamilyHaps(fam))         # 
          fam$hapIncidences <- ""                                                # 
        } else {
          fam <- initializeIDS(fam)                                              # IDS matrix creation
          if (any(fam$idHaps[-(1:2),]==3))  {                                    # At least one child has haplotypes non-inherited from a parent
            fam$hapIncidences <- "Irregular inheritance detected"                # Irregular inheritance due to recombination event, genotyping error or inheritance from non-declared parent
            names(fam$imputedMkrs) <- names(fams$phasedMkrs)                     # 
            fams$phasedMkrs <- rbind(fams$phasedMkrs,fam$imputedMkrs)            # 
            fams$haplotypes <- rbind(fams$haplotypes,matrix(NA,fam$famSize,2))   # 
            fam$IDS <- array(0,dim=c(2,fam$nMkrs,fam$famSize))                   # 
          } else {
            fam <- famPhaser(fam)                                                               # Haplotype phasing for each family
            missParsMkrs=function(fam) apply(fam$phasedAlls[,,1:2],2,function(x) any(is.na(x))) #
            if (any(missParsMkrs(fam))&fam$famSize>4) fam <- solveHapsWithMissingParents(fam)   # For three or more children
            if (any(missParsMkrs(fam))&fam$famSize>3) fam <- phaseHapsFromTwoChildren(fam)      # For some cases with two children
            if (length(case5Mks$mrk)>0) fam <- includeCase5Mkrs(fam,fam0,case5Mks)              #
            HMZupdate=function(fam) 
              apply(fam$phasedAlls,3,function(x) apply(x,2,function(a) as.integer(a[1]==a[2]))) #
            fam$HMZ <- t(HMZupdate(fam))                                                        #
            if (nrow(fam$HMZ)==1) fam$HMZ=t(fam$HMZ)                                            #
            names(fam$phasedMkrs) <- names(fams$phasedMkrs)                                     # 
            fams$phasedMkrs <- rbind(fams$phasedMkrs,fam$phasedMkrs)                            # 
            fams$haplotypes <- rbind(fams$haplotypes,fam$haplotypes)                            # 
            fam$hapIncidences <- ""                                                             # 
          }
        }
        idsCount <- apply(fam$IDS,3,sum)                                                      # 
        V1 <- fam$nMkrs*2-idsCount                                                            # 
        fullHapsInds <- which(idsCount==nMkrs*2)                                              # 
        partialHapsInds <- which(idsCount<nMkrs*2&idsCount>0)                                 # 
        V2 <- rep(0,fam$famSize)                                                              # 
        if (length(fullHapsInds)>0) V2[fullHapsInds] <- 2                                     # 
        if (length(partialHapsInds)>0) V2[partialHapsInds] <- 1                               # 
        V3 <- ifelse(apply(fam$idHaps,1,sum)>0,1,0)                                           # 
        fam$Vdata <- data.frame(nonAlls=V1,fullHaps=V2,IDSindiv=V3)                           # 
        rownames(fam$Vdata) <- NULL                                                           #  
        
        fams$Vdata <- rbind(fams$Vdata,fam$Vdata)                                             #
        fams$hapIncidences <- rbind(fams$hapIncidences,c(as.character(f),fam$hapIncidences))  #
      }
      
      #=== FAMILIES BOUND ===# 
      fams$phasedMkrs <- cbind(fams$imputedMkrs[,1:6],fams$phasedMkrs)         # 
      fams$haplotypes <- cbind(fams$haplotypes,fams$Vdata)                     # 
      fams$hapIncidences <- as.data.frame(fams$hapIncidences)                  # 
      colnames(fams$hapIncidences) <- c("Family","Incidence")                  #
      
      return(fams)
    }  
    # Data Summary 
    summarizeData=function(fams){
      nAlls <- (ncol(fams$imputedMkrs)-6)*nrow(fams$imputedMkrs)                       # Total number of alleles
      pNonPhasAlls <- sum(fams$haplotypes$nonAlls)/nAlls                               # Number of non-phased alleles
      pPhasAlls <- 1-pNonPhasAlls                                                      # Number of phased alleles
      nHaps <- nrow(fams$haplotypes)*2                                                 # Number of haplotypes
      pFullHaps <- sum(fams$haplotypes$fullHaps)/nHaps                                 # Number of phased haplotypes
      pFullNAhaps <- length(which(fams$haplotypes$fullHaps==0))/nHaps                  # Number of missing haplotypes
      pPartialHaps <- length(which(fams$haplotypes$fullHaps==1))/nHaps                 # No. partially phased haplotypes
      fams$phasingSummary <- data.frame(nAlls,pPhasAlls,pNonPhasAlls,nHaps,pFullHaps,    
                                        pFullNAhaps,pPartialHaps,phasingTime)          # Phasing summary
      return(fams)
    }
  }
  ################################################ II. IMPUTATION ###############################################
  {
    fams <- alleImputer(data,invisibleOutput,dataSummary)        # Preliminary imputation (marker by marker)
  }
  ################################################# III. PHASING ################################################
  {
    phasingTime <- system.time(fams <- famsPhaser(fams))[[3]]    # Phasing time    
  }
  ############################################### IV. DATA SUMMARY  #############################################
  {
    fams <- summarizeData(fams)                                  # Phased data summary
  }
  ###############################################  V. DATA STORING  #############################################
  {
    if (is.character(data)) {     
      baseName <- file_path_sans_ext(data)
      outName1 <- paste(baseName,"phased.ped",sep="_")
      outName2 <- paste(baseName,"haplotypes.txt",sep="_")
      write.table(fams$imputedMkrs,sep=" ",quote=FALSE,file=outName1)
      if (any(fams$hapIncidences[,-1]!="")) {
        outName3 <- paste(baseName,"hapIncidences.txt",sep="_")
        write.table(fams$hapIncidences,sep=" ",quote=FALSE,file=outName3)     
      }
    }    
  }
  ############################################## VI. FUNCTION OUTPUT  ###########################################
  {
    if (dataSummary==TRUE) {
      
      ### Phasing message printing
      if (all(fams$hapIncidences[,-1]=="")) {
        cat("\n\nHaplotypes have been successfully phased!!!")
        fams$hapIncidences <- NULL
      } else {
        cat("\n\nIncidences were detected. Some haplotypes could not be phased properly.")
      }
      
      ### Data summary printing
      cat("\n\n=========== PHASING SUMMARY ============")
      cat(paste("\nProportion of phased alleles:",round(fams$phasingSummary$pPhasAlls,4)))
      cat(paste("\nProportion of non-phased alleles:",round(fams$phasingSummary$pNonPhasAlls,4)))
      cat(paste("\nProportion of missing haplotypes:",round(fams$phasingSummary$pFullNAhaps,4)))
      cat(paste("\nProportion of partial haplotypes:",round(fams$phasingSummary$pPartialHaps,4)))
      cat(paste("\nProportion of full haplotypes:",round(fams$phasingSummary$pFullHaps,4)))
      cat(paste("\nPhasing time:",round(fams$phasingSummary$phasingTime,4)))
      cat("\n=======================================\n")
      
      ### Storage path
      if (is.character(data)) {
        cat("\nPhased data have been stored in: \n")
        cat(paste("\n",outName1,"\n",sep="")); cat(outName2)
        if (any(fams$hapIncidences[,2]!="")) cat(paste("\n",outName3,sep=""))
      } #else cat(getwd())
    }     
    
    ### Cleaning empty haplotype Incidences
    if (all(fams$hapIncidences[,-1]=="")) {
      fams$hapIncidences <- NULL
    }

    ### Returning (or not) the output
    if (invisibleOutput) return(invisible(fams)) else return(fams)              #         
  }

}