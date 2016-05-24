###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2014
#
initpar <- 
function(Xlist, ylist, Slist, nalist, k, m, p, bscov, control) {
#
################################################################################
#
  # CHECK Psifix AND Psicor IF NEEDED 
  if(bscov%in%c("prop","fixed")) {
    if(!is.null(fix <- control$Psifix)) {
      Psifix <- if(is.matrix(fix)) fix else xpndMat(fix)
      Psifix <- checkPD(Psifix,k,"Psifix")
    }
  } else if(bscov=="cor") {
    Psicor <- inputcov(rep(1,k),control$Psicor)
    Psicor <- checkPD(Psicor,k,"Psicor")
  } 
#
  # CHECK AND RETURN STARTING VALUES (IF PROVIDED)
  if(!is.null(initPsi <- control$initPsi)) {
    if(is.vector(initPsi)) initPsi <- xpndMat(initPsi)
    initPsi <- checkPD(initPsi,k,"initPsi")
  # IF NOT PROVIDED AND NEEDED, DETERMINE Psi WITH IGLS
  } else if(!bscov%in%c("fixed")) {
    initPsi <- diag(0.001,k)
    for(i in seq(control$igls.iter)) 
      initPsi <- iter.igls(initPsi,Xlist,ylist,Slist,nalist,k,m)
  }
#
################################################################################
#
  # PARAMETERS OF Psi DEPENDENT ON STRUCTURE AND PARAMETERIZATION
  par <- switch(bscov,
    # IF UNSTRUCTURED, LOWER TRIANGULAR OF THE CHOLESKY DECOMPOSITION
    unstr = vechMat(t(chol(initPsi))),
    # DIAGONAL: THE LOGARITHM OF THE DIAGONAL ELEMENTS
    diag = log(diag(initPsi)),
    # IDENTITY: THE LOGARITHM OF THE MEAN OF THE DIAGONAL ELEMENTS
    id = log(mean(diag(initPsi))),
    # COMPOUND SYMMETRY (WITH POSITIVE-DEFINITENESS CONSTRAINT)
    cs = {
      if(k<2) stop("bscov='cs' only meaningful with more than 1 outcome")
      cor <- mean(cov2cor(initPsi)[row(initPsi)!=col(initPsi)])
      if(cor <= -1/(k-1)) cor <- -1/k
      c(log(mean(diag(initPsi)))/2, log((cor+1/(k-1))/(1-cor)))
    },
    # HETEROGENEOUS COMPOUND SYMMETRY (WITH POSITIVE-DEFINITENESS CONSTRAINT)
    hcs = {
      if(k<2) stop("bscov='hcs' only meaningful with more than 1 outcome")
      cor <- mean(cov2cor(initPsi)[row(initPsi)!=col(initPsi)])
      if(cor <= -1/(k-1)) cor <- -1/k
      c(log(diag(initPsi))/2, log((cor+1/(k-1))/(1-cor)))
    },
    # AUTOREGRESSIVE OF FIRST ORDER
    ar1 = {
      if(k<2) stop("bscov='hcs' only meaningful with more than 1 outcome")
      cor <- mean(cov2cor(initPsi)[row(initPsi)-col(initPsi)==1L])
      c(log(diag(initPsi)), (qlogis((cor+1)/2)))
    },
    # PROPORTIONAL: THE LOGARITHM OF THE MEAN RATIO OF THE DIAGONAL ELEMENTS
    prop = log(mean(diag(Psifix/initPsi))),
    # KNOWN CORRELATION: SAME AS DIAGONAL
    cor = log(diag(initPsi)),
    # FIXED
    fixed = NULL
  )
#
  par
}
