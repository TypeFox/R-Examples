##' Generic correlation method
##'
##' @title Generic method for extracting correlation coefficients of model object
##' @param x Object
##' @param \dots Additional arguments
##' @author Klaus K. Holst
##' @export
"correlation" <- function(x,...) UseMethod("correlation")

##' @export
correlation.lvmfit <- function(x,z=TRUE,level=0.05,adj=TRUE,...) {
  pp <- matrices2(Model(x), with(index(x),seq_len(npar+npar.mean+npar.ex)))$P
  pos <- pp[lower.tri(pp)][(index(x)$P0)[lower.tri(pp)]==1]
  if (length(pos)<1) return(NULL)
  pp0 <- pp
  pp0[index(x)$P0!=1] <- 0; pp0[lower.tri(pp0)] <- 0
  coords <- c()
  mynames <- vars(x)
  n <- nrow(pp0)
  res <- c()
  for (i in pos) {
    idx <- which(pp0==i)
    rowpos <- (idx-1)%%n + 1
    colpos <- ceiling(idx/n)
    coefpos <- c(i,pp0[rbind(c(rowpos,rowpos),c(colpos,colpos))])
    pval <- pp[rbind(c(rowpos,rowpos),c(colpos,colpos))]
    phi.v1.v2 <- numeric(3);
    newval <- coef(x)[coefpos]
    phi.v1.v2[coefpos!=0] <- newval
    phi.v1.v2[coefpos==0] <- pval[tail(coefpos==0,2)]
    f <- function(p) {
      p[1]/sqrt(p[2]*p[3])
    }
    rho <- f(phi.v1.v2)
    if (z) {
        zrho <- atanh(rho)
        var.z <- 1/(nrow(model.frame(x))-ifelse(adj,3,0)) ## n-k-3
        ci.z <- zrho + c(-1,1)*qnorm(1-level/2)*sqrt(var.z)
        ci.rho <- tanh(ci.z)
        zs <- 1/sqrt(var.z)*zrho
        p.z <- 2*(pnorm(-abs(zs))) # p-value using z-transform for H_0: rho=0.
        est <- c(rho, NA, ci.rho[1], ci.rho[2])
    } else {
      Sigma.phi.v1.v2 <- matrix(0,3,3)
      Sigma.phi.v1.v2[coefpos!=0,coefpos!=0] <- vcov(x)[coefpos,coefpos]
      nabla.f <- function(p) {
        c(1/sqrt(p[2]*p[3]), -f(p)/(2*p[2]), -f(p)/(2*p[3]))
      }
      rho.var <- t(nabla.f(phi.v1.v2))%*%Sigma.phi.v1.v2%*%nabla.f(phi.v1.v2)
      est <- c(rho, sqrt(rho.var),  rho + c(-1,1)*qnorm(1-level/2)*sqrt(rho.var))
    }
    res <- rbind(res,est)
    rownames(res)[nrow(res)] <- paste(mynames[c(rowpos,colpos)],collapse="~")
  }
  colnames(res) <- c("Correlation","Std.Err","lowerCI","upperCI")
  return(res)
}
