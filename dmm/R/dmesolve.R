dmesolve <-
function(mdf,fixform = Ymat ~ 1,components=c("VarE(I)","VarG(Ia)"),cohortform=NULL,posdef=T,gls=F,glsopt=list(maxiter=200,bdamp=0.8,stoptol=0.01), dmeopt="qr",ncomp.pcr="rank",relmat="inline",dmekeep=F,dmekeepfit=F) {
# dmesolve() - dyadic model equations solved by ols and (optionally) gls
  #
  # fixform is the model formula for fixed effects
  # df is dataframe, k is fixed effects, l is traits
  # k and l determined from dataframe during AOV step
  # cohortform is formula for constructing cohort environment coding from (multiple) 
  #            columns in df
  # v is no of individual var/cov components (incl residual)
  # m & n determined from dataframe
  # m = no of individuals in pedigree
  # n = no of individuals with data
  #
  # mdf can be a df or a list containig a df and a list of relmat's
  # df must include Id, SId, DId columns,
  #     and they must be numeric and class "NULL",
  #     and Id must be sorted in ascending order and without duplicates,
  #     and each DId and SId must occur as an Id,
  #     as prepared by mdf2() function
  # df must include a matrix of the Y variables which must be called "Ymat"
  #     as prepared by mdf2() function
  # df must include any X variates or codes mentioned in fixform
  #     and these must be of class "factor" if codes (but not if covariates)
  #     as prepared by mdf2() function
  # df must include a Sex column if using sex-linked component partitioning
  #     and this must be of class "factor"
  #
# make component tables
  ctable <- make.ctable()
# check components in ctable
  if(any(is.na(match(components,ctable$all)))){
    print(components)
    stop("Component(s) not recognized:\n")
  }
#
  if(is.null(mdf$rel)) {
    df <- mdf
    cat("Data file is a normal dataframe:\n")
    if(relmat == "withdf") {
      stop("dmm: cant have 'relmat=withdf' option for a normal dataframe\n")
    }
  }
  else {
    df <- mdf$df
    cat("Data file is a list containing a dataframe and a list of relationship matrices:\n")
  }
# cat("Model formula for fixed effects:\n")
# print(fixform)
  cat("Random effect partitioned into components: Residual:\n")
  v <- length(components)
  cat("No of individual variance components partitioned(v):",v,"\n")
# print(components)
# cat("Model formula for Cohort environment:\n")
# print(cohortform)
#
# AOV step:
# initial aov() of fixed effects only
#
  fixed.aov <- aov(fixform,df,x=T,y=T,qr=T)
  cat("OLS-b step:\n")
# cat("AOV step:\n")
# cat("Initial aov of fixed effects only\n")
# print(fixed.aov)
# cat("OLS summary(aov):\n")
# print(summary(fixed.aov))
# print(anova(fixed.aov))
# fixed.uaov <- drop1(fixed.aov)  # only works in univariate case
# cat("OLS coef(aov):\n")
# print(coef(fixed.aov))
# cat("fixed.aov$coef\n")
# print(fixed.aov$coef)
# cat("fixed.aov$effects\n")
# print(fixed.aov$effects)
# cat("fixed.aov$y:\n")
# print(fixed.aov$y)
# model.tables fails for mean only model and for multivariate case
# cat("Tables.aov:\n")
# print(model.tables(fixed.aov))
  k <- ncol(fixed.aov$x)
  l <- ncol(as.matrix(fixed.aov$y))
  cat("no of fixed effect df (k) = ",k,"\n")
  cat("no of traits (l) = ",l,"\n")
# cat("AOV substep completed:\n")
  aov.list <- list(aov=fixed.aov, mdf=substitute(mdf), fixform=fixform)
#

# cohort setup
# cat("Cohort setup:\n")
  if(any(!is.na(match(components,ctable$cohort)))){
    cohortlabs <- attributes(terms(cohortform))$term.labels
    celabs <- c("DId",cohortlabs)
#   cat("cohortlabs:\n")
#   print(cohortlabs)
#   cat("celabs:\n")
#   print(celabs)
    cohortparts <- match(cohortlabs,colnames(df))
    ceparts <- match(celabs, colnames(df))
  }
  else {
    ceparts <- NULL
    cohortparts <- NULL
    celabs <- "DId"
    cohortlabs = NULL
  }
#   cat("cohortparts:\n")
#   print(cohortparts)
#   cat("ceparts:\n")
#   print(ceparts)

#  setup am object
# am is an antemodel object containing k,l,m,n,r,v,Z1,..Zr,Rel,..Rel,X,Y
  cat("Setup antemodel matrices:\n")
  am <- am.zandrel(mdf,df,k,l,v,as.matrix(fixed.aov$x),as.matrix(fixed.aov$y),
                   cohortparts,components,relmat,ctable)
  if( l == 1) {
#   dimnames(am$y) <- list(NULL,"Ymat")   #list(NULL,dimnames(df[[2]])[dimnames(df)[[2]] == "Ymat"]])
    dimnames(am$y) <- list(NULL,as.character(terms(fixform)[[2]]))
  }
# above change made to cope with univariate case

  cat("no of individuals in pedigree (m) = ",am$m,"\n")
  cat("no of individuals with data and X codes (n) = ",am$n,"\n")
#  cat("Z matrix\n")
#  print(am$z)
  am.list <- list(am=am, components=components, v=v, cohortform=cohortform, cohortlabs=cohortlabs,
              celabs=celabs, cohortparts=cohortparts, ceparts=ceparts)
#  cat("Antimodel matrices completed:\n")

#
# OLS AOV analysis -> initial estimates of b - why is this repeated??
# cat("AOV followup:\n")
  x.qr <- fixed.aov$qr
  b <- qr.coef(x.qr,am$y)
  krank <- x.qr$rank
  if(krank < am$k) {
    stop("Rank of X ",krank," .ne. no of fixed effects ",am$k,"\n")
  }
# cat("OLS AOV Estimates of b:\n")
# print(b)  # b here is k x l
  cat("Rank of X:",x.qr$rank,"  No of Fixed Effects:",k,"\n")

# # Setup (Y-Xb) matrix - residuals from OLS AOV
#  cat("Setup vec(Y-Xb) matrix:\n")
   ymxb <- am$y - am$x %*% b
#  cat("dim of ymxb matrix:",dim(ymxb),"\n")
#  cat("Y-Xb matrix:\n")
#  print(ymxb)
#

# SS after OLS AOV fit
# cat("Calculate residual variance after OLS fit:\n:")
  ssa <- t(ymxb) %*% ymxb
# cat("ssa:\n")
# print(ssa)
  degf <- am$n - am$k
  vara <- ssa/degf
# cat("Adjusted variance of individuals after OLS (Y-Xb)'(Y-Xb):\n")
#   print(vara)
#   cat("Above adjusted variances are observed phenotypic variances for this sample:\n")

# Var of b estimates
# cat("Calculate var/cov matrix of b estimates:\n")
# vb <- kronecker(vara, solve(crossprod(qr.R(x.qr))), make.dimnames=T)
  # note crossprod works for multivariate case - these 2 code lines are equivalent
  vb <- kronecker(vara, solve(t(qr.R(x.qr)) %*% qr.R(x.qr)), make.dimnames=T)
  # vb is in lxl blocks each kxk
# cat("Variances of estimates of b:\n")
# print(vb)
# cat("SE's of estimates of b:\n")
  seb <- matrix(sqrt(diag(vb)), am$k, am$l, dimnames=dimnames(b))
# print(seb)
# cat("AOV followup completed:\n")

  ols.fixed.list <- list(b=b,seb=seb,vara=vara,totn=am$n,degf=degf)

#
#
# DME step:
#   individual components after OLS fit
  cat("DME substep:\n")
#  cat("Dyadic equations for (Y-Xb) with b from OLS:\n")
# setup evec - data vector matching emat
   evec <- kronecker(ymxb,ymxb,make.dimnames=T)
#
# cat("(ymxb)(ymxb)'\n")
# print(evec)
#


# setup emat - expectation matrix (also emat.qr and zaz = vmat)
    dyad.explist <- dyad.am.expect(am,gls,dmeopt) 

#
# calc siga from evec and emat
    degfd <- am$n * am$n - am$v
#
    if(dmekeep) {
      dme.exp.list <- list(dme.wmat=dyad.explist$emat,dme.psi=evec,dme.mean=dyad.explist$emat.mean,dme.var=dyad.explist$emat.var, dme.correl=dyad.explist$emat.cor)
    } 
    else{
      dme.exp.list <- list(dme.mean=dyad.explist$emat.mean,dme.var=dyad.explist$emat.var,dme.correl=dyad.explist$emat.cor)
    }

# fit options for DME's
  if(dmeopt == "qr") {
    cat("QR option on dyadic model equations:\n")
    siga <- qr.coef(dyad.explist$emat.qr,evec)  # siga is r x l^2
    vard <- crossprod(qr.resid(dyad.explist$emat.qr,evec))
    vard <- vard/degfd
#   cat("Residual var for DME (vard):\n")
#   print(vard)
    vsiga <- kronecker(vard, solve(crossprod(qr.R(dyad.explist$emat.qr))), make.dimnames=T)
    sesiga <- matrix(sqrt(diag(vsiga)), am$v, am$l * am$l, dimnames=dimnames(siga))
#   cat("vsiga:\n")
#   print(vsiga)
#   cat("sesiga:\n")
#   print(sesiga)
    if(dmekeepfit) {
      dme.fit.list <- list(dme.fit=dyad.explist$emat.qr,dmeopt=dmeopt)
    }
    else {
      dme.fit.list <- list(dmeopt=dmeopt)
    }
  }

  else if(dmeopt == "lm") {
   dme.lm <- lm(evec ~ -1 + ., as.data.frame(dyad.explist$emat),x=T,y=T,qr=T)
#  dme.lm <- lm(evec ~  ., as.data.frame(dyad.explist$emat))
   cat("LM option on dyadic model equations:\n")
#  print(dme.lm)
#  cat("Summary of lm() on dyadic model eqns:\n")
#  print(summary(dme.lm))
#  print(anova(dme.lm))
#   cat("Estimates by lm:\n")
#   print(summary(dme.lm)[[1]]$coef[,1]) # [[1]] is first trait only
#   cat("Std.Error of estimates by lm:\n")
#   print(summary(dme.lm)[[1]]$coef[,2])
     # extract siga and sesiga
#    siga <- matrix(0, am$v, am$l * am$l,dimnames=list(dimnames(dyad.explist$emat)[[2]], dimnames(evec)[[2]]))
     siga <- matrix(0, am$v, am$l * am$l,dimnames=list(colnames(dyad.explist$emat), colnames(evec)))
     sesiga <- matrix(0, am$v, am$l * am$l, dimnames=dimnames(siga))
     if(am$l == 1) {
       siga[ ,1] <- summary(dme.lm)$coef[,1]
       sesiga[ ,1] <- summary(dme.lm)$coef[,2]
     }
     else {
       for (l2 in 1 : (am$l * am$l)) {
         siga[ ,l2] <- summary(dme.lm)[[l2]]$coef[,1]
         sesiga[ ,l2] <- summary(dme.lm)[[l2]]$coef[,2]
       }
     }
#   cat("sesiga:\n")
#   print(sesiga)
#   vard <- crossprod(resid(dme.lm))
    vard <- crossprod(qr.resid(dme.lm$qr,evec))
    vard <- vard/degfd
#   cat("Residual var for DME (vard):\n")
#   print(vard)
    vsiga <- kronecker(vard,solve(crossprod(qr.R(dme.lm$qr))),make.dimnames=T)
#   cat("vsiga: \n")
#   print(vsiga)
#   plot(profile(dme.lm))
#  cat("Sensitivity:\n")
#  print(influence.measures(dme.lm))
    if(dmekeepfit) {
      dme.fit.list <- list(dme.fit=dme.lm,dmeopt=dmeopt)
    } 
    else {
      dme.fit.list <- list(dmeopt=dmeopt)
    }
  }

  else if(dmeopt == "lmrob") {
   if(am$l > 1) {
     stop("Lmrob option does not work in multivariate case:\n")
   }
   dme.lmrob <- lmrob(evec ~ -1 + ., as.data.frame(dyad.explist$emat),x=T,y=T,qr=T)
#  dme.lmrob <- lmrob(evec ~  ., as.data.frame(dyad.explist$emat))
   cat("LMROB option on dyadic model equations:\n")
#  print(dme.lmrob)
#  cat("Summary of lmrob() on dyadic model eqns:\n")
#  print(summary(dme.lmrob))
#  print(anova(dme.lmrob))
     # extract siga and sesiga
     siga <- matrix(0, am$v, am$l * am$l,dimnames=list(colnames(dyad.explist$emat), colnames(evec)))
     sesiga <- matrix(0, am$v, am$l * am$l, dimnames=dimnames(siga))
     if(am$l == 1) {
       siga[ ,1] <- summary(dme.lmrob)$coef[,1]
       sesiga[ ,1] <- summary(dme.lmrob)$coef[,2]
     }
     else {
       for (l2 in 1 : (am$l * am$l)) {
         siga[ ,l2] <- summary(dme.lmrob)[[l2]]$coef[,1]
         sesiga[ ,l2] <- summary(dme.lmrob)[[l2]]$coef[,2]
       }
     }
#   cat("sesiga:\n")
#   print(sesiga)
#   vard <- crossprod(resid(dme.lmrob))
    vard <- crossprod(qr.resid(dme.lmrob$qr,evec))
    vard <- vard/degfd
#   cat("Residual var for DME (vard):\n")
#   print(vard)
    vsiga <- kronecker(vard,solve(crossprod(qr.R(dme.lmrob$qr))),make.dimnames=T)
#   cat("vsiga: \n")
#   print(vsiga)
    if(dmekeepfit) {
      dme.fit.list <- list(dme.fit=dme.lmrob, dmeopt=dmeopt)
    }
    else {
      dme.fit.list <- list(dmeopt=dmeopt)
    }
  }
  
  else if(dmeopt == "pcr"){
   if(ncomp.pcr == "all") {
     myncomp <- am$v
   }
   else if(ncomp.pcr == "rank"){
     myncomp <- dyad.explist$emat.qr$rank
   }
   else if(is.numeric(ncomp.pcr)){
     myncomp <- min(am$v, ncomp.pcr)
   }
   else{
     stop("Invalid option ncomp.pcr: ",ncomp.pcr,"\n")
   }
   dme.pcr <- mvr(evec ~ -1 + ., ncomp=myncomp,  data=as.data.frame(dyad.explist$emat),method="svdpc",validation="CV",model=T,x=T,y=T,jackknife=T)
#  dme.pcr <- mvr(evec ~  ., ncomp=myncomp, data=as.data.frame(dyad.explist$emat),method="svdpc",validation="CV",model=T,x=T,y=T,jackknife=T)
   ncomp <- dme.pcr$ncomp
   cat("PCR option on dyadic model equations:\n")
   summary(dme.pcr)
#  print(dme.pcr)
#  print(attributes(dme.pcr))
#   cat("Estimates by pcr:\n")
#   print(coef(dme.pcr))
#   cat("Scores by pcr:\n")
#   print(scores(dme.pcr))
#   cat("Loadings by pcr:\n")
#   print(loadings(dme.pcr))
    vsiga <- var.jack(dme.pcr,ncomp=myncomp,covariance=T)[,,1]
#   cat("Vsiga:\n")
#   print(vsiga - nearPD(vsiga,ensureSymmetry=T)$mat)
#   vsiga <- nearPD(vsiga,ensureSymmetry=T)$mat

# extract siga and sesiga
     siga <- matrix(coef(dme.pcr)[,,1], am$v, am$l * am$l,dimnames=list(colnames(dyad.explist$emat), colnames(evec)))
     sesiga <- matrix(sqrt(diag(vsiga)), am$v, am$l * am$l, dimnames=dimnames(siga))
#     cat("Sesiga:\n")
#     print(sesiga)
# residuals and their co/variances
    residmat <- resid(dme.pcr)[,,ncomp]
    vard <- crossprod(residmat)
    vard <- vard/degfd
#   cat("Residual var for DME (vard):\n")
#   print(vard)
    if(dmekeepfit) {
      dme.fit.list <- list(dme.fit=dme.pcr,dmeopt=dmeopt)
    }
    else {
      dme.fit.list <- list(pcr.loadings=loadings(dme.pcr),dmeopt=dmeopt)
    }
  }

# cat("Partitioned variance  components from DME equations:\n")
# print(siga)

#
# check siga[1,], siga[2,], siga[3,], siga[5,] are positive definite
# adjust siga[,] if not pd
  if(posdef){
    siga <- siga.posdef(siga,am,ctable)
#   cat("Partitioned variance components made positive definite:\n")
#   print(siga)
  }
  cat("DME substep completed:\n")


#  sum variance components (by DME method) to get phenotypic variances
# cat("Summing variance component estimates to get phenotypic component:\n")
# cat("Sums of individual variance components(by DME method):\n")
    covp <- apply(siga,2,sum)
#  reorganize sums into an lxl matrix
   covp <- matrix(covp,am$l,am$l,dimnames=list(colnames(am$y),colnames(am$y)))
#  print(covp)
#  cat("Above sums of components estimate population phenotypic variance:\n")

   ols.random.list <- list(siga=siga, sesiga=sesiga, vard=vard, degfd=degfd)
#
# Var/Covariances to genetic parameters and SE's
#  cat("Components to genetic parameters and SE's:\n")
   ols.genpar.list <- comtopar(am$v,am$l,siga,vara,vsiga,sesiga,ctable)
#  print(ols.genpar.list)

   ols.list <- c(aov.list, ols.fixed.list, dme.exp.list, dme.fit.list, ols.random.list, ols.genpar.list)
#  outlist <- list(ols=ols.list)
   outlist <- ols.list
   cat("OLS-b step completed:\n")

#
# GLS analysis
# GLS on b step
  if(gls && posdef) {
   cat("\nGLS-b step:\n")
#  cat("Fixed effects iterated -> GLS b\n")
#  cat("Partitioned variance components from DME after GLS b \n")
   if(am$l > 1) {
     cat("Warning: Multivariate GLS is not same as multiple univariate GLS's\n")
   }

#  GLS iteration of b's
   gls.list <- gls.iter.b(am, b, siga, dyad.explist, glsopt, dmeopt, ctable, ncomp.pcr, dmekeepfit)
#  print(gls.list)

     if(gls.list$ok) {
#      cat("Components to genetic parameters and SE's:\n")
       gls.genpar.list <- comtopar(am$v,am$l,gls.list$siga, gls.list$msa, gls.list$vsiga,gls.list$sesiga,ctable)
       gls.list.out <- list(b=gls.list$b, seb=gls.list$seb, siga=gls.list$siga, 
        sesiga=gls.list$sesiga, vard=gls.list$vard, msr=gls.list$msr,
        msrdf=gls.list$msrdf, msa=gls.list$msa)
       gls.list.out <- c(gls.list.out,gls.genpar.list,gls.list$dme.fit.list)
       # setup return object
#      outlist <- c(list(ols=ols.list), list(gls=gls.list.out))
       outlist <- c(ols.list, list(gls=gls.list.out))
       cat("GLS-b step completed successfully:\n")
    }
    else if (!gls.list$ok) {
#      outlist <- c(list(ols=ols.list))
       outlist <- ols.list
       cat("GLS-b step abandoned:\n")
    }
  }
  return(outlist)
}
