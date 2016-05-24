 

## The function is currently defined as
mat.jumprecon<-function(dObj, modEnvCol,fossEnvCol,cutoff=0, alpha=NULL,distance = -1, numanalogs=length(dObj$position[,1]))
{
   
inModern=get(dObj$inModern)
inFossil=get(dObj$inFossil)
longlat=dObj$llmod

if (numanalogs > length(dObj$position[,1])) stop(paste("Your dissimilarity object has less analogs than",numanalogs))

#DECLARE EMPTY VECTORS FOR OUTPUT
        outavgvec=vector('numeric')

#
#CREATE DISSIMILARITY MATRIX

        dimnames(dObj$sqdist)=NULL
        dimnames(dObj$position)=NULL

        sqmat=t(dObj$sqdist)
        numcols=ncol(sqmat)
        numrows=nrow(sqmat)
        distmat=t(dObj$distance)

#
#CREATE MATRIX OF RECONSTRUCTED ENVIRONMENTAL VARIABLES

        reconvec = inModern[,modEnvCol]
        reconmat = matrix(data=reconvec[t(dObj$position)],nrow=numrows,ncol=numcols)
        longvec = inModern[,longlat[1]]
        latvec = inModern[,longlat[2]]
        longmat = matrix(data=longvec[t(dObj$position)],nrow=numrows,ncol=numcols)
        latmat = matrix(data=latvec[t(dObj$position)],nrow=numrows,ncol=numcols)
        
        numeratorStartCol = 2
        denominatorStartCol = 1
        numeratorEndCol = numcols
        demonimatorEndCol = numcols-1

#
#IF CUTOFF>0 THEN REMOVE ALL RECONSTRUCTED VALUES > THE DISSIMILARITY THREASHOLD

        if (cutoff > 0) {

                reconmat[sqmat>=cutoff]=NA

        }

        if (distance > -1 ){
        reconmat[distmat>distance]=NA
        
        avgvec=apply(reconmat[,1:numanalogs],1,mean,na.rm=T)
        outrecon=cbind(x=dObj$x,y=dObj$y,avgvec)
        }


        #
        #MAKE THE JUMP MATRIX   
 
        jumpmat = ((sqmat[,numeratorStartCol: numeratorEndCol]/sqmat[,denominatorStartCol: demonimatorEndCol]) - 1) * 100
        jumpmat[is.infinite(jumpmat)]=0
        numrows=nrow(jumpmat)

        #
        #ITERATE OVER JUMPS FROM ZERO TO ONE-HUNDRED PERCENT

        
                
        tfmat=jumpmat<alpha
        tfmat=cbind(rep(TRUE,length(tfmat[,1])),tfmat)

        for(k in 1:numrows){

        jumppoint = match(FALSE,tfmat[k,])

                if(!is.na(jumppoint)) {
                tfmat[k,jumppoint:numcols]=FALSE                           
                }

        }
        
        errorplus=apply(reconmat,1,max)
        errorminus=apply(reconmat,1,min)
        
        reconmat[!tfmat]=NA
        sqmat[!tfmat]=NA
        #outrecon=cbind(dObj$long,dObj$lat,outavgvec)
        avgvec=rowMeans(reconmat,na.rm=T)
        avgdissim=rowMeans(sqmat,na.rm=T)
        
        #errorminus=avgvec-errorplus    
#errorplus=errorplus+avgvec
        outrecon=cbind(x=dObj$x,y=dObj$y,recon=avgvec,avgdissim=avgdissim,ep=errorplus,em=errorminus)

#
#OUTPUT THE CORRELATIONS AND JUMP VALUES FOR PLOTTING

#list(jumpmat,tfmat,sqmat,reconmat,cbind(inFossil[,fossEnvCol],outrecon))

cbind(inFossil[,fossEnvCol],outrecon)

  }



