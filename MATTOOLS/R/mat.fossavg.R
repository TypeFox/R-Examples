 
mat.fossavg <- function(dObj, modEnvCol, fossEnvCol=0, fossCols=0, cutoff=0, distance = -1, wmethod="none", numanalogs=length(dObj$position[,1]))
{
   
inModern=get(dObj$inModern)

longlat=dObj$llmod

if (numanalogs > length(dObj$position[,1])) stop(paste("Your dissimilarity object has less analogs than",numanalogs))

#DECLARE EMPTY VECTORS FOR OUTPUT
        outcorr=vector('numeric')
        outavgvec=vector('numeric')

#
#CREATE DISSIMILARITY MATRIX

        dimnames(dObj$sqdist)=NULL
        dimnames(dObj$position)=NULL

        sqmat=t(dObj$sqdist)
        distmat=t(dObj$distance)
        posmat=t(dObj$position)

        numcols=ncol(sqmat)
        numrows=nrow(sqmat)
        
#
#CREATE MATRIX OF RECONSTRUCTED ENVIRONMENTAL VARIABLES

        reconvec = inModern[,modEnvCol]
        reconmat = matrix(data=reconvec[t(dObj$position)],nrow=numrows,ncol=numcols)
        longvec = inModern[,longlat[1]]
        latvec = inModern[,longlat[2]]
        longmat = matrix(data=longvec[t(dObj$position)],nrow=numrows,ncol=numcols)
        latmat = matrix(data=latvec[t(dObj$position)],nrow=numrows,ncol=numcols)
        
#
##GET NAME OF ENV VARIABLE

        recname=names(inModern)[modEnvCol]

#
#IF CUTOFF>0 THEN REMOVE ALL RECONSTRUCTED VALUES > THE DISSIMILARITY THREASHOLD

        if (cutoff > 0) {

                reconmat[sqmat>cutoff]=NA
                sqmat[sqmat>cutoff]=NA

        }

        if (distance > -1 ){
        
                reconmat[distmat>distance]=NA
                sqmat[distmat>distance]=NA
                #avgvec=apply(reconmat[,1:numanalogs],1,mean,na.rm=T)
                #outrecon=cbind(x=dObj$x,y=dObj$y,recname=avgvec)

        }

#
#GET WEIGHTED AVERAGE VECTOR AND RECONMAT
        
        if(wmethod=="inv.dissim" & numanalogs > 1){
                
                weightmat=1/sqmat[,1:numanalogs]
                weightmat[is.na(reconmat[,1:numanalogs])]=NA
                #weightmat[reconmat[,1:numanalogs]==NA]=NA
                numerators=reconmat[,1:numanalogs]*weightmat
                sqavg=rowMeans(sqmat[,1:numanalogs],na.rm=T)
                dissweight=rowSums(numerators,na.rm=T)/rowSums(weightmat,na.rm=T)
                theSD=sqrt(rowSums(weightmat*(reconmat[,1:numanalogs]-dissweight)^2,na.rm=T) / (rowSums(weightmat,na.rm=T))  )
                outrecon=cbind(x=dObj$x,y=dObj$y,recname=dissweight,SDMINUS=dissweight-theSD,SDPLUS=dissweight+theSD,avgdissim=sqavg)           
        
        }
        else if (wmethod=="inv.dist" & numanalogs > 1){

                dmat=distmat
                dmat[distmat==0]=1
                weightmat=1/dmat[,1:numanalogs]
                #weightmat[reconmat[,1:numanalogs]==NA]=NA
                weightmat[is.na(reconmat[,1:numanalogs])]=NA
                sqavg=rowMeans(sqmat[,1:numanalogs],na.rm=T)
                numerators=reconmat[,1:numanalogs]*weightmat
                dissweight=rowSums(numerators,na.rm=T)/rowSums(weightmat,na.rm=T)
                theSD=sqrt(rowSums(weightmat*(reconmat[,1:numanalogs]-dissweight)^2,na.rm=T) / (rowSums(weightmat,na.rm=T))  )
                outrecon=cbind(x=dObj$x,y=dObj$y,recname=dissweight,SDMINUS=dissweight-theSD,SDPLUS=dissweight+theSD,avgdissim=sqavg)           

        }
        else if (wmethod=="inv.rank" & numanalogs > 1){

                weightmat=1/matrix(1:numanalogs,numrows,length(1:numanalogs),byrow=T)
                #weightmat[reconmat[,1:numanalogs]==NA]=NA
                weightmat[is.na(reconmat[,1:numanalogs])]=NA
                sqavg=rowMeans(sqmat[,1:numanalogs],na.rm=T)
                numerators=reconmat[,1:numanalogs]*weightmat
                dissweight=rowSums(numerators,na.rm=T)/rowSums(weightmat,na.rm=T)
                theSD=sqrt(rowSums(weightmat*(reconmat[,1:numanalogs]-dissweight)^2,na.rm=T) / (rowSums(weightmat,na.rm=T))  )
                outrecon=cbind(x=dObj$x,y=dObj$y,recname=dissweight,SDMINUS=dissweight-theSD,SDPLUS=dissweight+theSD,avgdissim=sqavg)   
                #outrecon=weightmat

        }
        else if (wmethod =="equal.wt" & numanalogs > 1){

                avgvec=apply(reconmat[,1:numanalogs],1,mean,na.rm=T)
                theSD=apply(reconmat[,1:numanalogs],1,sd,na.rm=T)
                sqavg=rowMeans(sqmat[,1:numanalogs],na.rm=T)

                outrecon=cbind(x=dObj$x,y=dObj$y,recname=avgvec,SDMINUS=avgvec-theSD,SDPLUS=avgvec+theSD,avgdissim=sqavg)

        }
        else if (wmethod =="single" & numanalogs == 1){

                outrecon=cbind(x=dObj$x,y=dObj$y,recname=reconmat[,1],SDMINUS=reconmat[,1],SDPLUS=reconmat[,1],dissim=sqmat[,1],rowinmod=posmat[,1],recx=longmat[,1],recy=latmat[,1])

        }
        else if (wmethod =="none" & numanalogs >1){

                outrecon=cbind(x=dObj$x,y=dObj$y,reconmat[,1:numanalogs],sqmat[,1:numanalogs],longmat[,1:numanalogs],latmat[,1:numanalogs],posmat[,1:numanalogs] )
                                namesvec=c("x","y",paste(names(inModern)[modEnvCol],1:numanalogs,sep=""),
                        paste("dVal",1:numanalogs,sep=""),paste("x",1:numanalogs,sep=""),paste("y",1:numanalogs,sep=""),paste("rowinmod",1:numanalogs,sep=""))

                colnames(outrecon)=namesvec     

        }

#
#OUTPUT THE CORRELATIONS AND JUMP VALUES FOR PLOTTING
if (length(fossCols)!=0){

        fossilout=get(dObj$inFossil)[,fossCols]

  }
##list(weightmat=weightmat,numerators=numerators,dissweight=dissweight,outrecon=outrecon)
outrecon=cbind(fossilout,outrecon)

attributes(outrecon)=c(attributes(outrecon),"numanalogs"=numanalogs,"weightmethod"=wmethod,"cutoff"=cutoff,"distance"=distance,"date"=date())
#list(a=outrecon,b=sqmat,c=reconmat)
#sqmat
#write.table(cbind(outrecon,sqmat,reconmat,weightmat,numerators),"c:\tst.txt",quote=F,row.names=F)
return(outrecon)

  }



