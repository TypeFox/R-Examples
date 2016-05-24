TOS2D <-function (image, detrend = FALSE, nsamples = 100, theTS = avespecvar, verbose = TRUE,...)
{

dname <- deparse(substitute(image))
statname <- deparse(substitute(theTS))

x <- cropimage(image)

# detrend if necessary:
if(detrend){	
    	x <- medpolish(x)$resid    
}

oa<-list(...)
si <-"smooth"%in%names(oa)

if(si){
	xcddews <- cddews(x,...)
}
else{
	xcddews <- cddews(x,smooth=FALSE,...)
}
TS <- theTS(xcddews)

m <- rowMeans(xcddews$S)
m[m < 0] <- 0
p <- array(m,dim=dim(xcddews$S))
   
xavspec2 <- xcddews
xavspec2$S <- p

for (b in 1:nsamples) {
	if(verbose){
		if((b%%50)==0){
                	cat("bootstrap number",b,"\n")
        	}
	}
        xbs2 <- LS2Wsim(xavspec2)
	if(si){
        	xbs2.cddews <- cddews(xbs2,...)
	}
	else{
        	xbs2.cddews <- cddews(xbs2,smooth=FALSE,...)
	}
        TS <- c(TS, theTS(xbs2.cddews))
}

if (verbose){
	cat("\n")
}


p.value <-getpval(TS,verbose=verbose)

l<-list(data.name=dname,samples=TS,statistic=statname,p.value=p.value)

class(l)<-"TOS2D"

return(l)
    
}

