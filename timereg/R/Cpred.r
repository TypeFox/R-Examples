Cpred<-function(cum,xval,start.val=0,cum.startval=0,order=FALSE,comp="smaller",strict=TRUE) 
{
designX<-as.matrix(cum); 
cumtimes <- designX[,1]
px<-as.integer(dim(designX)[2]);
nx<-as.integer(dim(designX)[1]);
nval<-length(xval); 
pred<-rep(0,nval); 

###sout<-.C("Cpred", 
###as.double(cumtimes),as.integer(nx),as.integer(px),
###as.double(xval),as.integer(nval),as.integer(pred),
###as.integer(Tminus),PACKAGE="timereg")

### sindex from prodlim
xval.order <- sindex.prodlim(cumtimes,xval,comp=comp,strict=strict)
pred.begin <-  xval.order
pred.begin[xval.order==0] <- 1
###predcum <- as.matrix(designX[pred.begin,-1])
predcum <- designX[pred.begin,-1,drop=FALSE]
predcum[xval.order==0,] <- cum.startval

if (order==FALSE) return(cbind(xval,predcum)) else return(list(xval.order=xval.order,pred.begin=pred.begin))
}

## sindex fra prodlim, thanks to Thomas Gerds
sindex.prodlim <- function (jump.times, eval.times, comp = "smaller", strict = FALSE)
{   
    stopifnot(is.numeric(jump.times))
    stopifnot(is.numeric(eval.times))
        N <- length(jump.times)
    if (comp == "greater") {
       N - sindex.prodlim(jump.times = jump.times, eval.times = eval.times, comp = "smaller", strict = !strict)
    }
    else {
           neval <- length(eval.times)
           if (!(neval > 0 && N > 0))
	                stop("missing data")
           new.order <- order(eval.times)
           ind <- .C("sindex", index = integer(neval), as.double(sort(jump.times)),
                  as.double(eval.times[new.order]), as.integer(N),
                  as.integer(neval), as.integer(strict), PACKAGE = "timereg")$index
		   ind[order(new.order)]
    }
}
