qqunif <- function(u,type="unif",logscale=TRUE,base=10,col=palette()[4],lcol=palette()[2],ci=FALSE,alpha=0.05,...)
{
  u <- u[!is.na(u)]
  n <- length(u)
  xlabel <- ifelse(logscale,paste("-log",base,"(Expected)",sep=""), "Expected")
  ylabel <- ifelse(logscale,paste("-log",base,"(Observed)",sep=""), "Observed")
  if (ci) c <- abs(qnorm(alpha/2))
  if (type=="exp") {
     n1 <- 1/(n:1)
     n2 <- n1^2
     lambda <- 1/log(base)
     m <- cumsum(n1)*lambda
     z <- qqplot(m,-log(u,base),xlab=xlabel,ylab=ylabel,col=col,...)
     if (ci)
     {
        v <- cumsum(n2)*lambda^2
        s <- sqrt(v)
        lcl <- m - c*s
        ucl <- m + c*s
        lines(m,lcl)
        lines(m,ucl)
     }
  } else if (type=="unif") {
     m <- (1:n)/(n+1)
     if (logscale) z <- qqplot(-log(m,base),-log(u,base),xlab=xlabel,ylab=ylabel,col=col,...)
     else z <- qqplot(m,u,xlab=xlabel,ylab=ylabel,col=col,...)
     if (ci)
     {
        v <- (1:n)*(n-(1:n)+1)/(n+1)^2/(n+2)
        s <- sqrt(v)
        lcl <- m - c*s
        ucl <- m + c*s
        lid <- (lcl>0)
        uid <- (ucl<=1)
        if (logscale)
        {
           a <- -log(m[lid],base)
           b <- -log(lcl[lid],base)
           c <- -log(m[uid],base)
           d <- -log(ucl[uid],base)
        } else {
           a <- m[lid]
           b <- lcl[lid]
           c <- m[uid]
           d <- ucl[uid]
        }
        lines(a,b)
        lines(c,d)
     }
  } else stop ("invalid type")
  abline(0,1,col=lcol)
# polygon(c(m,rev(m)),c(ucl,rev(lcl)),col="gray")
  invisible(z)
}

#09-11-2009 refine U(0,1) and exp(lambda)
#27-08-2008 first attempt of CI
