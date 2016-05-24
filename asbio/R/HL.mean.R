HL.mean <- function (x, conf = NULL, method = "exact") 
{
    d <- outer(x, x, "+")/2
    sus <- d[lower.tri(d, diag = TRUE)]
    pm <- median(sus)

if(!is.null(conf)){
    alpha <- 1 - conf
    if(method == "normal"){
        n <- length(x)
        c.alpha <- floor((n*(n+1))/4-qnorm(1-(alpha/2))*sqrt((n*(n+1)*(2*n+1))/24))
        M <- (n^2+n)/2 
        soso <- sort(sus)
        ci.L <- soso[c.alpha]
        ci.U <- soso[M + 1 - c.alpha]
    }
    
    if(method == "exact"){
        n <- length(x)
        M <- (n^2+n)/2 
        c.alpha <- M + 1 - (qsignrank(1-(alpha/2), n = n)+1)
        soso <- sort(sus)
        ci.L <- soso[c.alpha]
        ci.U <- soso[M + 1 - c.alpha]
        }
}

if(is.null(conf)) res <- pm
if(!is.null(conf)){ 
res <- list()
res$head<-paste(paste(as.character(conf*100),"%",sep=""),c("Confidence interval for population pseudomedian"))
res$ci<-c(pseudomedian=pm,lower=ci.L,upper=ci.U)
res$ends<-c("Estimate",paste(as.character(c((1-conf)/2,1-((1-conf)/2))*100),"%",sep=""))
class(res)<-"ci"
}
res
}