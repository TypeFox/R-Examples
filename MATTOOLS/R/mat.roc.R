 

## The function is currently defined as
mat.roc <-function(inModern, modTaxa=c(), colClasses=NULL,numAnalogs=2,rocEvalSeq=seq(0,2,0.05),counts=F,aucmethod="trap")
{
if(modTaxa[1] >= modTaxa[2]) stop("The number specifying the first column containing the taxa must be smaller than the number specifying the last column in modTaxa=c(FirstTaxon,LastTaxon).") 

if (counts) {

inModern[,modTaxa]=inModern[,modTaxa]/rowSums(inModern[,modTaxa])
if(is.na(sum(rowSums(inModern[,modTaxa])))) stop("Your taxon dataset has missing values or other values that are causing problems.  Check for missing values.")

  }

#  1.  Get the number of biological zones and create a list with these zones named plus an overall zone set

        biozones=as.character(inModern[,colClasses])
        biozones=unique(biozones)
        numzones=length(biozones)
        outlist=vector(mode="list",length=numzones+1)
        names(outlist)=c(biozones,"Overall")
        midlist=vector(mode="list")

#  2. Create two vectors one to hold all inzone2inzone A+ analogs and one to hold all outzone2inzone A- non-analogs

        allzone2zone = vector("numeric")
        alloutzone2zone = vector("numeric")

#  2. for each zone

        for(i in 1:numzones){
                inzonesub=inModern[,colClasses]==biozones[i]
                inzone1 = inModern[inzonesub,]
                inzone2 = inModern[!inzonesub,]
                outzone2inzone = mattools.roc(inFossil=inzone2, inModern=inzone1, modTaxa,numAnalogs=numAnalogs)[1:numAnalogs-1,]
                inzone2inzone = mattools.roc(inFossil=inzone1, inModern=inzone1, modTaxa,numAnalogs=numAnalogs)[2:numAnalogs,]
                allzone2zone=c(allzone2zone,inzone2inzone)
                alloutzone2zone=c(alloutzone2zone,outzone2inzone)

        #  3. For current biological zone undertake ROC analysis and store result in appropriate list component

                currentroclist=mat.ROCcalc(truth=c(rep(0,length(outzone2inzone)),rep(1,length(inzone2inzone))), data=c(outzone2inzone,inzone2inzone), evalseq=rocEvalSeq,method=aucmethod)
                midlist=list(zonezone=inzone2inzone,outzone=outzone2inzone,ROC=currentroclist,within=length(inzone1[,1]),outside=length(inzone2[,1]))
                outlist[[i]]=midlist    
        }

#  4. For overall zonezone and outzone undertake final ROC analysis and put in the.

        currentroclist=mat.ROCcalc(truth=c(rep(0,length(alloutzone2zone)),rep(1,length(allzone2zone))), data=c(alloutzone2zone,allzone2zone), evalseq=rocEvalSeq,method=aucmethod)
        midlist=list(zonezone=allzone2zone,outzone=alloutzone2zone,ROC=currentroclist,within=length(allzone2zone),outside=length(alloutzone2zone))
        outlist[[numzones+1]]=midlist   

#  5. Create a summary table
  
        trecs=names(outlist)
        tcols=names(outlist[[1]][[3]])[c(6,8,9)]
        outmat=matrix(NA,nrow=length(trecs),ncol=6)
        
        for (i in 1:(numzones+1)) {
                crows=c(outlist[[i]][[3]]$optDissVal[1],outlist[[i]][[3]]$AUC,outlist[[i]][[3]]$SEAUC,length(outlist[[i]][[3]]$optDissVal),outlist[[i]][[4]],outlist[[i]][[5]])
                outmat[i,]=crows
                

        }
        
        rownames(outmat)=trecs
        colnames(outmat)=c(tcols,"nsol","nWithin","nOutside")
        outlist=c(outlist,roctable=NULL)
        outlist$roctable=outmat
        
        
print(outmat)
outlist

  }



