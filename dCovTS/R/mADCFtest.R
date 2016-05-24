mADCFtest <- function(x,type = c("truncated", "bartlett", "daniell", "QS",    
 "parzen"), p, b = 0, parallel = FALSE)
{
 type <- match.arg(type)
 data.name <- deparse(substitute(x))
 if (!is.matrix(x))stop('Only multivariate time series with dimension d>2')
 if (!is.numeric(x)) 
     stop("'x' must be numeric")
 if(!all(is.finite(x))) stop('Missing or infitive values')
 n <- as.integer(NROW(x))
 q <- as.integer(NCOL(x))
 MaxLag <- n-1
 t <- rep(0,MaxLag)
 for(k in 1:MaxLag){
  kern <- kernelFun(type,k/p)
  if (kern !=0){
   t[k] <- (n-k)*kern^2*sum(mADCF(x,lags=k,output=FALSE)^2)
  }
 }
 stat <- sum(t)
 if(!b==0){
 Tnstar <- TstarBoot2(x,type,p,b,parallel)
 pvalue <- sum(Tnstar>=stat)/(b+1)
 }
 p.value <- ifelse(b == 0,NA,pvalue)
  if(b==0){
   Tnstar <- NULL
  } else Tnstar <- Tnstar
 dataname <- paste(data.name,","," kernel type: ", type,", bandwidth=",p, ", boot replicates ", b, sep = "")
 names(stat) <- "Tnbar"
   e=list(method = paste("Multivariate test of independence based on distance correlation", sep = ""), 
        statistic = stat, p.value = p.value, replicates=Tnstar,data.name=dataname)
 class(e) <- "htest"
 return(e)
}

