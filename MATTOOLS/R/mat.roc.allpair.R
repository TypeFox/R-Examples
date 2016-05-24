 
mat.roc.allpair<- function(inModern, modTaxa=c(), colClasses=NULL,numAnalogs=2,rocEvalSeq=seq(0,2,0.05),counts=F,aucmethod="trap")
{
if(modTaxa[1] >= modTaxa[2]) stop("The number specifying the first column containing the taxa must be smaller than the number specifying the last column in modTaxa=c(FirstTaxon,LastTaxon).") 

if (counts) {

inModern[,modTaxa]=inModern[,modTaxa]/rowSums(inModern[,modTaxa])

  }

#  1.  Get the number of biological zones and create a list with these zones named plus an overall zone set

        biozones=as.character(inModern[,colClasses])
        biozones=unique(biozones)
        numzones=length(biozones)
        outlist=vector(mode="list",length=(numzones*numzones))
        #names(outlist)=c(biozones,"Overall")
        midlist=vector(mode="list")
        masterzonevec=rep(biozones,each=numzones)
        trecs=vector("character")
#  2. Create two vectors one to hold all inzone2inzone A+ analogs and one to hold all outzone2inzone A- non-analogs

        allzone2zone = vector("numeric")
        alloutzone2zone = vector("numeric")
        #overalloutzone=vector("numeric")

#  2. for each zone

        counter=1

        for(k in 1:numzones){

                inzonesub=inModern[,colClasses]==biozones[k]
                inzone1 = inModern[inzonesub,]
                inzone2inzone = mattools.roc(inFossil=inzone1, inModern=inzone1, modTaxa, numAnalogs=numAnalogs)[2:numAnalogs,]

                for(i in (1:numzones)[-k]){
                        inzonesub=inModern[,colClasses]==biozones[i]            
                        inzone2 = inModern[inzonesub,]
                        outzone2inzone = mattools.roc(inFossil=inzone2, inModern=inzone1, modTaxa,numAnalogs=numAnalogs)[1:numAnalogs-1,]
                        #allzone2zone=c(allzone2zone,inzone2inzone)
                        alloutzone2zone=c(alloutzone2zone,outzone2inzone)

        #  3. For current biological zone undertake ROC analysis and store result in appropriate list component

                        currentroclist=mat.ROCcalc(truth=c(rep(0,length(outzone2inzone)),rep(1,length(inzone2inzone))), data=c(outzone2inzone,inzone2inzone), evalseq=rocEvalSeq,method=aucmethod)
                        midlist=list(zonezone=inzone2inzone,outzone=outzone2inzone,ROC=currentroclist)
                        outlist[[counter]]=midlist      
                        trecs=c(trecs,biozones[i])
                        counter=counter+1
                }

        currentroclist=mat.ROCcalc(truth=c(rep(0,length(alloutzone2zone)),rep(1,length(inzone2inzone))), data=c(alloutzone2zone,inzone2inzone),         evalseq=rocEvalSeq,method=aucmethod)
        midlist=list(zonezone=inzone2inzone,outzone=alloutzone2zone,ROC=currentroclist)
        outlist[[counter]]=midlist
        trecs=c(trecs,"Overall")
        alloutzone2zone = vector("numeric")     
        counter=counter+1
        }

#  5. Create a summary table
  
        tcols=names(outlist[[1]][[3]])[c(6,8,9)]
        outmat=matrix(NA,nrow=length(trecs),ncol=4)
        for (i in 1:((numzones*numzones))) {
                crows=c(outlist[[i]][[3]]$optDissVal[1],outlist[[i]][[3]]$AUC,outlist[[i]][[3]]$SEAUC,length(outlist[[i]][[3]]$optDissVal))
                outmat[i,]=crows
                

        }
        
        rownames(outmat)=paste(masterzonevec,trecs)
        colnames(outmat)=c(tcols,"nsol")
        outlist=c(outlist,roctable=NULL)
        outlist[[i+1]]=outmat

print(outmat)
names(outlist)=c(paste(masterzonevec,"_vs._",trecs,sep=""),"roctable")
outlist

  }



