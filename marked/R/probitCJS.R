#' Perform MCMC analysis of a CJS model
#' 
#' Takes design data list created with the function \link{make.design.data} for model "probitCJS" 
#' and draws a sample from the posterior distribution using a Gibbs sampler.
#' 
#' 
#' @param ddl list of dataframes for design data; created by call to
#' \code{\link{make.design.data}}
#' @param dml list of design matrices created by \code{\link{create.dm}} from
#' formula and design data
#' @param parameters A model specification list with a list for Phi and p containing a 
#'        formula and optionally a prior specification which is a named list. See 'Priors' section below.
#' @param design.parameters Specification of any grouping variables for design
#' data for each parameter
#' @param burnin number of iteration to initially discard for MCMC burnin
#' @param iter number of iteration to run the Gibbs sampler for following burnin
#' @param initial A named list (Phi,p). If null and imat is not null, uses cjs.initial to create initial values; otherwise assigns 0
#' @param imat A list of vectors and matrices constructed by \code{\link{process.ch}} from the capture history data
#' @section Prior distribution specification:
#'        The prior distributions used in \code{probitCJS} are multivatiate normal with mean mu a
#'        and variance tau*(X'X)^{-1} on the probit scale for fixed effects. The matrix X is 
#'        the design matrix based on the model specification (located in \code{parameters$Phi$formula} and
#'        \code{parameters$p$formula} respectively). Priors for random effect variance 
#'        components are inverse gamma with shape parameter 'a' and rate parameter 'b'. Currently, 
#'        the default values are mu=0 and tau=0.01 for phi and p parameters. For all randome effects
#'        deault values are a=2 and b=1.0E-4. In addition to the variance component each random effect
#'        can be specified with a known unscaled covariance matrix, Q,  if random effects are correlated.
#'        For example, to obtain a random walk model for a serially correlated effect see 
#'        Examples section below. To specify different hyper-parameters for the prior distributions, it must be done 
#'        in the \code{parameters} list. See the Examples section for changing the prior distributions. Note that a and b can be
#'        vectors and the Qs are specified via a list in order of the random effects specified in the 
#'        model statements.
#'         
#' @return A list with MCMC iterations and summarized output:
#' \item{beta.mcmc}{list with elements Phi and p which contain MCMC iterations for each beta parameter} 
#' \item{real.mcmc}{list with elements Phi and p which contain MCMC iterations for each real parameter} 
#' \item{beta}{list with elements Phi and p which contain summary of MCMC iterations for each beta parameter} 
#' \item{reals}{list with elements Phi and p which contain summary of MCMC iterations for each real parameter} 
#' @export
#' @import coda truncnorm Matrix
#' @author Devin Johnson
#' @examples
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' # Analysis of the female dipper data
#' data(dipper)
#' dipper=dipper[dipper$sex=="Female",]
#' # following example uses unrealistically low values for burnin and 
#' # iteration to reduce package testing time
#' fit1 = crm(dipper,model="probitCJS",model.parameters=list(Phi=list(formula=~time),
#'  p=list(formula=~time)), burnin=100, iter=1000)
#' fit1
#' # Real parameter summary
#' fit1$results$reals
#' 
#' # Changing prior distributions:
#' fit2 = crm(dipper,model="probitCJS",
#'    model.parameters=list(
#'      Phi=list(formula=~time, prior=list(mu=rep(0,6), tau=1/2.85^2)),
#'      p=list(formula=~time, prior=list(mu=rep(0,6), tau=1/2.85^2))
#'    ), burnin=100, iter=1000)
#' fit2
#' # Real parameter summary
#' fit2$results$reals
#' 
#' # Creating a Q matrix for random walk effect for 6 recapture occasions
#' A=1.0*(as.matrix(dist(1:6))==1)
#' Q = diag(apply(A,1,sum))-A
#' 
#' # Fit a RW survival process
#' fit3 = crm(dipper,model="probitCJS",
#'    model.parameters=list(
#'      Phi=list(
#'        formula=~(1|time), 
#'        prior=list(mu=0, tau=1/2.85^2, re=list(a=2, b=1.0E-4, Q=list(Q)))
#'        ),
#'      p=list(formula=~time, prior=list(mu=rep(0,6), tau=1/2.85^2))
#'    ), burnin=100, iter=1000)
#' fit3
#' # Real parameter summary (Not calculated for random effects yet)
#' fit3$results$reals
#' 
#' }
probitCJS = function(ddl,dml,parameters,design.parameters,burnin, iter, initial=NULL, imat=NULL){
  
  ### DEFINE SOME FUNCTIONS ###
  
  sample.z = function(id, mu.y, mu.z, yvec){
    .Call("sampleZ", ID=id, PVec=pnorm(mu.y), PhiVec=pnorm(mu.z), yVec=yvec, PACKAGE="marked")
  }
  make.ztilde.idx = function(id, zvec){
    .Call("makeZtildeIdx", ID=id, zvec=zvec, PACKAGE="marked")
  }
  
  #browser()
  
  ### DATA MANIPULATION ###
  ##  restrict design data so Time>=Cohort and recreate design matrices 
  restricted.dml = create.dml(ddl,parameters,design.parameters,restrict=TRUE) 
  is.phi.re=!is.null(restricted.dml$Phi$re)
  is.p.re=!is.null(restricted.dml$p$re)
  yvec = as.vector(t(ddl$ehmat[,-1]))
  yvec=yvec[ddl$p$Time>=ddl$p$Cohort]
  n=length(yvec)
  ddl$p = ddl$p[ddl$p$Time>=ddl$p$Cohort,]
  Xy = as.matrix(restricted.dml$p$fe)
  pn.p = colnames(Xy)
  Xz = as.matrix(restricted.dml$Phi$fe)
  pn.phi = colnames(Xz)
  id = ddl$p$id
  if(is.phi.re){
    n.phi.re=sapply(restricted.dml$Phi$re$re.list, function(x){dim(x)[2]})
    K.phi=do.call("cBind", restricted.dml$Phi$re$re.list)
    ind.phi.re=rep(1:length(n.phi.re), n.phi.re)
    if(length(n.phi.re)>1){
      m.phi.re=model.matrix(~factor(ind.phi.re)-1)
    } else{
      m.phi.re=matrix(1, nrow=n.phi.re[[1]])
    }
    colnames(m.phi.re)=names(n.phi.re)
  }
  if(is.p.re){
    n.p.re=sapply(restricted.dml$p$re$re.list, function(x){dim(x)[2]})
    K.p=do.call("cBind", restricted.dml$p$re$re.list)
    ind.p.re=rep(1:length(n.p.re), n.p.re)
    if(length(n.p.re)>1){
      m.p.re=model.matrix(~factor(ind.p.re)-1)
    } else m.p.re=matrix(1, nrow=n.p.re[[1]])
    colnames(m.p.re)=names(n.p.re)
  }
  ###  PRIOR DISTRIBUTIONS ###
  if(is.null(parameters$Phi$prior)){
    tau.b.z = 0.01
    mu.b.z = rep(0,ncol(Xz))
  }else{
    tau.b.z = parameters$Phi$prior$tau
    mu.b.z = parameters$Phi$prior$mu
  }		
  if(is.null(parameters$p$prior)){
    tau.b.y = 0.01
    mu.b.y = rep(0,ncol(Xy))	
  }else{
    tau.b.y = parameters$p$prior$tau
    mu.b.y = parameters$p$prior$mu
  }
  if(is.phi.re){
    #n.phi.re=length(dml$Phi$re)
    if(is.null(parameters$Phi$prior$re)){
      a.phi=rep(2,length(n.phi.re))
      b.phi=rep(1.0E-4,length(n.phi.re))
      Q.phi=lapply(n.phi.re, function(x){Diagonal(x)})
      rnks.phi=n.phi.re
    }
    else{
      a.phi=parameters$Phi$prior$re$a
      b.phi=parameters$Phi$prior$re$b
      Q.phi=parameters$Phi$prior$re$Q
      if(!is.list(Q.phi)) {Q.phi = list(Q.phi)}
      rnks.phi=unlist(lapply(Q.phi, function(x){sum(eigen(x)$values>nrow(x)*.Machine$double.eps)}))
    }
    Q.phi=.bdiag(Q.phi)
    dimnames(Q.phi)=list(colnames(K.phi), colnames(K.phi))
  }
  if(is.p.re){
    #n.p.re=length(dml$p$re)
    if(is.null(parameters$p$prior$re)){
      a.p=rep(2,length(n.p.re))
      b.p=rep(1.0E-4,length(n.p.re))
      Q.p=lapply(n.p.re, function(x){Diagonal(x)})
      rnks.p=n.p.re
    }else{
      a.p=parameters$p$prior$re$a
      b.p=parameters$p$prior$re$b
      Q.p=parameters$p$prior$re$Q
      if(!is.list(Q.p)) {Q.p = list(Q.p)}
      rnks.p=unlist(lapply(Q.p, function(x){sum(eigen(x)$values>nrow(x)*.Machine$double.eps)}))
    }
    Q.p=.bdiag(Q.p)
    dimnames(Q.p)=list(colnames(K.p), colnames(K.p))
  }
  
  ### Initial values
  if(is.null(initial)) {
    beta = cjs.initial(dml,imat=imat,link="probit")
  } else beta = set.initial(c("Phi","p"),dml,initial)$par
  beta.z = beta$Phi  
  beta.y = beta$p
  
  # browser()
  
  if(is.phi.re){
    alpha.phi = rep(0,ncol(K.phi))
    eta.phi = as.vector(K.phi%*%alpha.phi)
    tau.phi = a.phi/b.phi
    Tau.phi.mat = Matrix(diag(x=rep(tau.phi, n.phi.re)))
    Q.alpha.phi = suppressMessages(Tau.phi.mat%*%Q.phi)
  } else eta.phi = rep(0,nrow(Xz))
  if(is.p.re){
    alpha.p = rep(0,ncol(K.p))
    eta.p = as.vector(K.p%*%alpha.p)
    tau.p = a.p/b.p
    Tau.p.mat = Diagonal(x=rep(tau.p, n.p.re))
    Q.alpha.p = suppressMessages(Tau.p.mat%*%Q.p)
  } else eta.p = rep(0,nrow(Xy))
  
  ### STORAGE ###
  beta.z.stor = matrix(NA, iter, ncol(Xz))
  beta.y.stor = matrix(NA, iter, ncol(Xy))
  colnames(beta.z.stor) = pn.phi
  colnames(beta.y.stor) = pn.p
  if(is.phi.re){
    alpha.phi.stor = matrix(NA, iter, length(alpha.phi))
    colnames(alpha.phi.stor) = colnames(K.phi)
    eta.phi.stor = matrix(NA, iter, length(eta.phi))
    tau.phi.stor = matrix(NA, iter, length(n.phi.re))
  }
  if(is.p.re){
    alpha.p.stor = matrix(NA, iter, n.p.re)
    colnames(alpha.p.stor) = colnames(K.p)
    eta.p.stor = matrix(NA, iter, length(eta.p))
    tau.p.stor = matrix(NA, iter, length(n.p.re))
  }
  ### BEGIN MCMC ###
  message("probitCJS MCMC beginning...")
  message("p model = ", as.character(parameters$p$formula))
  message("phi model = ", as.character(parameters$Phi$formula))
  flush.console()
  
  # browser()
  
  tot.iter = burnin + iter
  st = Sys.time()
  for(m in 1:tot.iter){
    
    ### UPDATE Z ###
    zvec = sample.z(id=id, mu.y=as.vector(Xy%*%beta.y+eta.p), mu.z=as.vector(Xz%*%beta.z+eta.phi), yvec)
    
    ### UPDATE Z.TILDE ### 
    a = ifelse(zvec==0, -Inf, 0)
    b = ifelse(zvec==0, 0, Inf)
    z.tilde = truncnorm::rtruncnorm(n, a=a, b=b, mean=as.vector(Xz%*%beta.z+eta.phi), sd=1)
    
    ### BETA.Z UPDATE ###
    idx.z.tilde = make.ztilde.idx(id, zvec)
    Q.b.z = tau.b.z*crossprod(Xz[idx.z.tilde,])
    V.beta.z.inv = crossprod(Xz[idx.z.tilde,]) + Q.b.z
    m.beta.z = solve(V.beta.z.inv, crossprod(Xz[idx.z.tilde,],z.tilde[idx.z.tilde]-eta.phi[idx.z.tilde]) + crossprod(Q.b.z,mu.b.z))
    beta.z = m.beta.z + solve(chol(V.beta.z.inv), rnorm(ncol(Xz),0,1))
    if(m>burnin)beta.z.stor[m-burnin,] = beta.z
    
    ### UPDATE Y.TILDE ### 
    a = ifelse(yvec==0, -Inf, 0)
    b = ifelse(yvec==0, 0, Inf)
    y.tilde = truncnorm::rtruncnorm(n, a=a, b=b, mean=as.vector(Xy%*%beta.y+eta.p), sd=1)
    
    ### BETA.Y UPDATE ###
    Q.b.y = tau.b.y*crossprod(Xy[zvec==1,])
    V.beta.y.inv = crossprod(Xy[zvec==1,]) + Q.b.y
    m.beta.y = solve(V.beta.y.inv, crossprod(Xy[zvec==1,],y.tilde[zvec==1]-eta.p[zvec==1])+crossprod(Q.b.y,mu.b.y))
    beta.y = m.beta.y + solve(chol(V.beta.y.inv), rnorm(ncol(Xy),0,1))
    if(m>burnin) beta.y.stor[m-burnin,] = beta.y
    
    #browser()
    
    ### RANDOM EFFECT UPDATES ###
    if(is.phi.re){
      ### ALPHA.PHI UPDATE ###
      Xbeta = Xz[idx.z.tilde,]%*%beta.z
      V.alpha.phi.inv = crossprod(K.phi[idx.z.tilde,]) + Q.alpha.phi
      m.alpha.phi = solve(V.alpha.phi.inv, crossprod(K.phi[idx.z.tilde,],z.tilde[idx.z.tilde]-Xbeta))
      alpha.phi = m.alpha.phi + solve(chol(V.alpha.phi.inv), rnorm(ncol(K.phi),0,1))
      eta.phi = K.phi%*%alpha.phi
      ### TAU.PHI UPDATE
      quad.phi = suppressMessages(crossprod(m.phi.re*alpha.phi,Q.phi))%*%(m.phi.re*alpha.phi)/2
      tau.phi = rgamma(length(n.phi.re), rnks.phi/2 + a.phi, as.numeric(quad.phi) + b.phi)
      Tau.phi.mat = Diagonal(x=rep(tau.phi, n.phi.re))
      Q.alpha.phi = Tau.phi.mat%*%Q.phi
      if(m>burnin) {
        alpha.phi.stor[m-burnin,] = as.vector(alpha.phi)
        eta.phi.stor[m-burnin,] = as.vector(eta.phi)
        tau.phi.stor[m-burnin,] = as.vector(tau.phi)
      }
    }
    if(is.p.re){
      ### ALPHA.P UPDATE ###
      Xbeta = Xy[zvec==1,]%*%beta.y
      V.alpha.p.inv = crossprod(K.p[zvec==1,]) + Q.alpha.p
      m.alpha.p = solve(V.alpha.p.inv, crossprod(K.p[zvec==1,],y.tilde[zvec==1]-Xbeta))
      alpha.p = m.alpha.p + solve(chol(V.alpha.p.inv), rnorm(ncol(K.p),0,1))
      eta.p = K.p%*%alpha.p
      ### TAU.P UPDATE
      quad.p = suppressMessages(crossprod(m.p.re*alpha.p,Q.p))%*%(m.p.re*alpha.p)/2
      tau.p = rgamma(length(n.p.re), rnks.p/2 + a.p, as.numeric(quad.p) + b.p)
      Tau.p.mat = Diagonal(x=rep(tau.p, n.p.re))
      Q.alpha.p = Tau.p.mat%*%Q.p
      if(m>burnin) {
        alpha.p.stor[m-burnin,] = as.vector(alpha.p)
        eta.p.stor[m-burnin,] = as.vector(eta.p)
        tau.p.stor[m-burnin,] = as.vector(tau.p)
      }
    }
    
    ### TIMING OF SAMPLER ###
    if(m==30){
      tpi = as.numeric(difftime(Sys.time(), st, units="secs"))/30
      ttc = round((tot.iter-30)*tpi/3600, 2)
      if(ttc>=1) cat("Approximate time till completion: ", ttc, " hours\n")
      else cat("Approximate time till completion: ", ttc*60, " minutes\n")
    }
    if(100*(m/tot.iter) >= 10 & (100*(m/tot.iter))%%10==0) 
    {
      cat(100*(m/tot.iter), "% completed\n")
      flush.console()
    }
  }
  
  message("MCMC complete, processing output...")
  
  ### MAKE SOME OUTPUT ### jll 2 July changed to summarize betas as well 
  ### 21 Sept jll removed real computations; it is now in compute.real called from crm
  phibeta.mcmc = mcmc(beta.z.stor)
  #summ.phi = summary(phibeta.mcmc)
  hpd.phi = HPDinterval(phibeta.mcmc)
  beta.phi = data.frame(mode = apply(beta.z.stor, 2, mcmc_mode), mean=apply(beta.z.stor, 2, mean), 
                        sd=apply(beta.z.stor,2,sd),CI.lower=hpd.phi[,1], CI.upper=hpd.phi[,2])
  pbeta.mcmc = mcmc(beta.y.stor)
  #summ.p = summary(pbeta.mcmc)
  hpd.p = HPDinterval(pbeta.mcmc)
  beta.p = data.frame( mode = apply(beta.y.stor, 2, mcmc_mode), 
                       mean=apply(beta.y.stor, 2, mean), 
                       sd=apply(beta.y.stor,2,sd),
                       CI.lower=hpd.p[,1], CI.upper=hpd.p[,2])  
  if(is.phi.re){
    vc.phi=mcmc(1/tau.phi.stor)
    colnames(vc.phi) = names(n.phi.re)
    vc.phi.df=data.frame(mode=apply(vc.phi, 2, mcmc_mode),mean=apply(vc.phi, 2, mean),
                         sd=apply(vc.phi,2,sd),
                         CI.lower=HPDinterval(vc.phi)[,1], CI.upper=HPDinterval(vc.phi)[,2])
    re.struc.Phi=list(K=K.phi, Q=Q.phi)
    phialpha.mcmc=mcmc(alpha.phi.stor)
    hpd.phi = HPDinterval(phialpha.mcmc)
    alpha.phi = data.frame( mode = apply(alpha.phi.stor, 2, mcmc_mode), 
                            mean=apply(alpha.phi.stor, 2, mean), 
                            sd=apply(alpha.phi.stor,2,sd),
                            CI.lower=hpd.phi[,1], CI.upper=hpd.phi[,2])  
    
  } else {
    vc.phi=NULL
    vc.phi.df=NULL
    re.struc.Phi=NULL
    alpha.phi=NULL
    phialpha.mcmc=NULL
  }
  if(is.p.re){
    vc.p=mcmc(1/tau.p.stor)
    colnames(vc.p) = names(n.p.re)
    vc.p.df=data.frame(mode=apply(vc.p, 2, mcmc_mode),mean=apply(vc.p, 2, mean),
                       sd=apply(vc.p,2,sd),
                       CI.lower=HPDinterval(vc.p)[,1], CI.upper=HPDinterval(vc.p)[,2])
    re.struc.p=list(K=K.p, Q=Q.p)
    palpha.mcmc=mcmc(alpha.p.stor)
    hpd.p = HPDinterval(palpha.mcmc)
    alpha.p = data.frame( mode = apply(alpha.p.stor, 2, mcmc_mode), 
                          mean=apply(alpha.p.stor, 2, mean), 
                          sd=apply(alpha.p.stor,2,sd),
                          CI.lower=hpd.p[,1], CI.upper=hpd.p[,2])  
    
  } else {
    vc.p=NULL
    vc.p.df=NULL
    re.struc.p=NULL
    alpha.p=NULL
    palpha.mcmc=NULL
  }
  res=list(beta=list(Phi=beta.phi,p=beta.p),
           alpha=list(Phi=alpha.phi,p=alpha.p),
           var.comp=list(Phi=vc.phi.df, p=vc.p.df),
           beta.mcmc=list(Phi=phibeta.mcmc, p=pbeta.mcmc), 
           alpha.mcmc=list(Phi=phialpha.mcmc, p=palpha.mcmc),
           var.comp.mcmc=list(Phi=vc.phi, p=vc.p), 
           model_data=list(Phi.dm=dml$Phi$fe,p.dm=dml$p$fe), 
           re.struct=list(Phi=re.struc.Phi, p=re.struc.p)
  )
  #Make Phi.dm equal to cBind(dml$Phi$fe, K.phi) and include alpha.phi with the mcmc output
  class(res)=c("crm","mcmc","probitCJS")
  return(res)
}	### END OF FUNCTION ###
