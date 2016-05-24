 
## The function is currently defined as
mat.jump <-function(dObj, inModern, envColumn, cutoff=0, calib=T)
{

  #DECLARE EMPTY VECTORS FOR OUTPUT
  outcorr=vector('numeric')
  outavgvec=vector('numeric')

#
#CREATE DISSIMILARITY MATRIX

dimnames(dObj$sqdist)=NULL
dimnames(dObj$position)=NULL

sqmat=t(dObj$sqdist)
numcols=ncol(sqmat)
numrows=nrow(sqmat)

#
#CREATE MATRIX OF RECONSTRUCTED ENVIRONMENTAL VARIABLES

reconvec = inModern[,envColumn]
reconmat = matrix(data=reconvec[t(dObj$position)],nrow=numrows,ncol=numcols)

        if(calib) 
        {
                numeratorStartCol = 3
                numeratorEndCol = numcols
                denominatorStartCol = 2
                demonimatorEndCol = numcols-1
        }
        else
        {
                numeratorStartCol = 2
                denominatorStartCol = 1
                numeratorEndCol = numcols-1
                demonimatorEndCol = numcols
        }

#
#MAKE THE JUMP MATRIX   
 
jumpmat = ((sqmat[,numeratorStartCol: numeratorEndCol]/sqmat[,denominatorStartCol: demonimatorEndCol]) - 1) * 100
jumpmat[is.infinite(jumpmat)]=0

#
#IF CUTOFF>0 THEN REMOVE ALL RECONSTRUCTED VALUES > THE DISSIMILARITY THREASHOLD

        if (cutoff > 0) {

                reconmat[sqmat>=cutoff]=NA

        }

numrows=nrow(jumpmat)

#
#ITERATE OVER JUMPS FROM ZERO TO ONE-HUNDRED PERCENT

for(i in 0:100){
        
        tfmat=jumpmat>i 

        for(k in 1:numrows){

        jumppoint = match(TRUE,tfmat[k,])

                if(!is.na(jumppoint)) {
                                                        currow=c(reconmat[k,denominatorStartCol],reconmat[k,numeratorStartCol:numcols][1:jumppoint])
currow=currow[!is.na(currow)]
                     outavgvec=c(outavgvec,sum(currow)/length(currow))
                 }
                else
                {

                     currow=reconmat[k,denominatorStartCol:numcols]
currow=currow[!is.na(currow)]
                     outavgvec=c(outavgvec,sum(currow)/length(currow))

                }
        }
        
        tcor = cor(reconvec, outavgvec,"pairwise.complete.obs")^2
        outcorr = c(outcorr, tcor)
        outavgvec=vector('numeric')
        
  }

#
#OUTPUT THE CORRELATIONS AND JUMP VALUES FOR PLOTTING

tseq=0:100
maxalpha=max(outcorr)

list(alpha = tseq[which(outcorr==maxalpha)], alphacor=cbind(alpha=0:100,correlation=outcorr))

  }



