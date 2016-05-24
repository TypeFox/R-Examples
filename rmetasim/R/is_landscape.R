#Allan Strand 12/03/02
#
#This function is intended to check and see if a landscape is of the correct form
#
#
is.landscape <- function(Rland=NULL,verb=TRUE,exact=FALSE)
  {
    ok<-TRUE;
    if (!is.list(Rland))
      {
        if (verb) {print("Landscape not a list");}
        ok<-FALSE;
      }

    if (is.null(Rland$intparam))
      {
        if (verb) {print("intparam not found.");}
        ok<-FALSE;
      }

    if (is.null(Rland$switchparam))
      {
        if (verb) {print("switchparam not found.");}
        ok<-FALSE;
      }

    if (is.null(Rland$floatparam))
      {
        if (verb) {print("floatparam not found.");}
        ok<-FALSE;
      }

#check integer parameters
    
#check boolean parameters
    
#check float parameters
    
#check demography

    if (is.null(Rland$demography))
      {
        if (verb) {print("demography not found.");}
        ok<-FALSE;
      }
    else
      {
        #check if the matrices have the correct dimensions
        for (i in 1:length(Rland$demo$localdem))
          {
            if
            (
             (
              dim(Rland$demo$localdem[[i]]$LocalS)!=c(Rland$intparam$stages,Rland$intparam$stages) ||
              dim(Rland$demo$localdem[[i]]$LocalR)!=c(Rland$intparam$stages,Rland$intparam$stages) ||
              dim(Rland$demo$localdem[[i]]$LocalM)!=c(Rland$intparam$stages,Rland$intparam$stages)
              )
             )
              {
                if (verb) {
                  print("One or more of the local demography matrices is not of the correct dimensions")
                  print(paste("Rland$inparame$stages",Rland$intparam$stages))
                }
                ok <- FALSE
              }
            
            if (max(apply(Rland$demography$localdem[[i]]$LocalS,2,sum))>1)
              {
                if (verb) {print(paste("Local survival matrix",i,"has a column that sums to a number greater than one"))}
                ok <- FALSE
              }
          }
        for (i in 1:length(Rland$demo$epochs))
          {
            if (
                (length(Rland$demo$epochs[[i]]$Extinct)!=Rland$intparam$habitats) ||
                (length(Rland$demo$epochs[[i]]$Carry)!=Rland$intparam$habitats) ||
                (length(Rland$demo$epochs[[i]]$Localprob)!=Rland$intparam$numdemos) ||
#                (length(Rland$demo$epochs[[i]]$RndChooseProb)!=Rland$intparam$numepochs) ||
                (min(dim(Rland$demo$epochs[[i]]$S)==c(Rland$intparam$stages*Rland$intparam$habitats,Rland$intparam$stages*Rland$intparam$habitats))==0) ||
                (min(dim(Rland$demo$epochs[[i]]$R)==c(Rland$intparam$stages*Rland$intparam$habitats,Rland$intparam$stages*Rland$intparam$habitats))==0) ||
                (min(dim(Rland$demo$epochs[[i]]$M)==c(Rland$intparam$stages*Rland$intparam$habitats,Rland$intparam$stages*Rland$intparam$habitats))==0)
                )
              {
                if (verb) {print(paste("One or more of the epoch paramters is of incorrect dimension in epoch",i))}
                ok <- FALSE
              }

            if (Rland$switchparam$randdemo==1||length(Rland$demography$localdem)==1)#still needs to check when more than one localdem and not random placement of demographies
              {
                for (j in 1:length(Rland$demography$localdem))
                  {
                    strt <- seq(1,dim(Rland$demo$epochs[[i]]$S)[1],dim(Rland$demography$localdem[[j]]$LocalS)[1])
                    stp <- strt+(dim(Rland$demography$localdem[[j]]$LocalS)[1]-1)
                    slice <- cbind(strt,stp)
#                    print(slice)
                    tmpS <- Rland$demography$epochs[[i]]$S
                    if (Rland$intparam$habitats!=1)
                      for (l in 1:Rland$intparam$habitats)
                        {
                          tmpS[c(slice[l,1]:slice[l,2]),c(slice[l,1]:slice[l,2])] <- Rland$demography$localdem[[j]]$LocalS
                        }
                    if (max(apply(tmpS,2,sum))>1)
                      {
                        if (verb) {
                          tcls <- which(apply(tmpS,2,sum)>1)
                          print(paste("Columns",paste(tcls,collapse=", ")," in the landscape S matrix associated with localdem",j,"total more than 1"))}
                        ok <- FALSE
                      }
                    
                  }
              }
          }
      }
    
#check loci
    if (is.null(Rland$loci))
      {
        if (verb) {print("loci not found.");}
        ok<-FALSE;
      }
    else
      {
        #first see if the number of loci match the number given in intparam
        if (length(Rland$loci)!=Rland$intparam$locusnum)
          {
            if (verb) {print("conflict between size of loci object and size specified in $intparam");}
            ok <- FALSE 
          }
      }
    
#check individuals

    if (is.null(Rland$individuals))
      {
        if (verb) {print("no individuals section found.");}
        ok<-FALSE;
      }
    else
      {
        #check if popsize>0
        if (!(dim(Rland$individuals)[1]>0))
          {
            if (verb) {print("No individuals in this landscape: spontaneous generation is notallowed")}
            ok <- FALSE

          }
#check if the max number of classes corresponds to the number in demography        
        if (max(Rland$individuals[,1])>((Rland$intparam$habitats*Rland$intparam$stages)-1))
          {
            if (verb) {print("There are individuals in stages other than demography allows")}
            ok <- FALSE
          }
                                        #check if the number and ploidy of loci work out
        if (dim(Rland$individuals)[2]!=(landscape.democol() + sum(landscape.ploidy(Rland))))
          {
            if (verb) {print("The number of loci do not correspond to the number of columns in $individuals\nThis indicates a problem with locus number and/or ploidy")}
            ok <- FALSE
          }
                                        #check if birthdates are later than the current generation
        if (max(Rland$individuals[,3])>Rland$intparam$currentgen)
          {
            if (verb) {print("There are individuals born in the future")}
            ok <- FALSE
          }
      }
    return(ok)
  }

