"wawotest.default" <-
function(x,alter=c("greater","less"), ...){
 alter <- match.arg(alter) 
 N <- length(x)
 x <- na.omit(x)
 if( N !=length(x))
   cat(paste(N-length(x),"NA removed\n"))
 x <- scale(x)
 N <- length(x)
 x2 <- c(x[-1],x[1])
 S1 <- sum(x)
 S2 <- sum(x^2)
 S3 <- sum(x^3)
 S4 <- sum(x^4)
 a <- sum(x*x2)
 ea <- (S1^2-S2)/(N-1)
 va <- ((S2^2-S4)/(N-1))+
   ((S1^4-4*S1^2*S2+4*S1*S3+S2^2-2*S4)/
    ((N-1)*(N-2)))-((S1^2-S2)^2/(N-1)^2)
 za <- (a-ea)/sqrt(va)
 p <- pnorm(za)
 if(alter=="greater")
   p=1-p
 return(c(a=a,ea=ea,va=va,za=za,p=p)) 
}

