###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2014
#
par2Psi <- 
function(par, k, bscov, control) {
#
################################################################################
# 
  Psi <- switch(bscov,
    # IF UNSTRUCTURED, CROSSPRODUCT FROM CHOLESKY
    unstr = {
      L <- diag(0,k)
      L[lower.tri(L,diag=TRUE)] <- par
      tcrossprod(L)
    },
    # DIAGONAL: THE EXPONENTIAL OF THE PARAMETERS
    diag = diag(exp(par),k),
    # IDENTITY: THE EXPONENTIAL OF THE PARAMETER   
    id = diag(k)*exp(par[1]),
    # COMPOUND SYMMETRY
    cs = {
      R <- matrix((exp(par[2])-1/(k-1))/(exp(par[2])+1),k,k)
      R[row(R) == col(R)] <- 1
      exp(par[1]*2) * R
    },
    # HETEROSCEDASTIC COMPOUND SYMMETRY
    hcs = {
      R <- matrix((exp(par[k+1])-1/(k-1))/(exp(par[k+1])+1),k,k)
      R[row(R) == col(R)] <- 1
      D <- diag(sqrt(exp(par[seq(k)]*2)),k)
      D%*%R%*%D
    },
    # AUTOREGRESSIVE OF FIRST ORDER
    ar1 = {
      cor <- plogis(par[k+1])*2-1
      R <- cor^abs(outer(seq(k),seq(k),"-"))
      D <- diag(sqrt(exp(par[seq(k)])),k)
      D%*%R%*%D
    },
    # PROPORTIONAL
    prop = {
      Psifix <- if(!is.null(fix <- control$Psifix)) {
        if(is.matrix(fix)) fix else xpndMat(fix)
      } else diag(k)
      exp(par)*Psifix
    },
    # KNOWN CORRELATION
    cor = inputcov(sqrt(exp(par)),control$Psicor),
    # FIXED
    fixed = {
      Psifix <- if(!is.null(fix <- control$Psifix)) {
        if(is.matrix(fix)) fix else xpndMat(fix)
      } else diag(0,k)
      Psifix
    }
  )
#
  Psi
}
