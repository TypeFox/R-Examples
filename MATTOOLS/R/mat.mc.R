 

## The function is currently defined as
mat.mc <- function (inModern,modTaxa=c(NULL,NULL),probs=c(0.05,0.025,0.01,0.001),freqint=seq(0, 2, 0.02),sampleSize=length(inModern[,1]),method="sawada",withReplace=T,counts=F) 
{

outvec = vector("numeric")

        if (counts) {

                inModern[,modTaxa]=inModern[,modTaxa]/rowSums(inModern[,modTaxa])

        }

if(method == "sawada") {

        set1=inModern[sample(1:length(inModern[,1]),sampleSize,withReplace),modTaxa]
        set2=inModern[sample(1:length(inModern[,1]),sampleSize,withReplace),modTaxa]
        set1=sqrt(set1)
        set2=sqrt(set2)
        set3=set1-set2
        set3=set3*set3
        sqvec=rowSums(set3)

  }
else if (method=="bartlein"){

        sqvec = mattools.roc(inModern,inModern,modTaxa,numAnalogs=length(inModern[,1]))

  }
      
#GET CUMULATIVE FREQUENCIES
    tlen=length(sqvec)
    for (i in freqint) {
        outvec = c(outvec, length(sqvec[sqvec <= i])/tlen)
    }

cumcurve = cbind(sqdist=freqint, problteq=outvec)
critvalue=approx(cumcurve[,2],cumcurve[,1],probs)

list(sqdist=sqvec,cumcurve=cumcurve,cutoffs=critvalue,method=method,samplesize=sampleSize,replacement=withReplace,probabilities=probs,wascounts=counts)

  }



