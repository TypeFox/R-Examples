#' Fit a spatial occupancy model using Gibbs sampling
#' 
#' This function fits a spatial occupancy model where the true occupancy is a
#' function of a spatial process.  An efficient Gibbs sampling algorithm is
#' used by formulating the detection and occupancy process models with a probit
#' model instead of the traditional logit based model. 
#' 
#' A Gibbs sampler is run to draw an MCMC sample of the spatial occupancy
#' parameters \code{beta} (detection parameters), \code{gamma} (the occupancy
#' parameters), \code{psi} (the model occupancy generating process), and the
#' realized occupancy.
#' 
#' @param detection.model A formula object describing the detection portion of
#' the occupancy model. The variables described by the detection model are
#' located in the \code{visit} data frame of the \code{so.data}.
#' @param occupancy.model A formula object describing the fixed effects portion
#' of the spatial occupancy process. The variables described by the occupancy
#' model are located in the \code{site} data frame of an \code{so.data} object.
#' @param spatial.model A named list object describing the spatial component of
#' the occupancy process. Currently the only possible models are ICAR, restricted spatial regression,
#' process convolution models, and no spatial model (i.e., eta = 0). Thus, \code{spatial.model=list(model="icar",
#' threshold= , rho=1)}, \code{spatial.model=list(model="rsr", threshold=, moran.cut=)}, 
#' \code{spatial.model=list(model="proc.conv", knots=)}, and \code{spatial.model=list(model="none")}
#' are the only forms that are accepted at present.  The \code{threshold}
#' component is used the create neighborhoods in the ICAR and RSR models.  All sites
#' within distance \code{threshold} of site i are considered neighbors of site
#' i.  The \code{moran.cut} component is the cut-off for selecting 
#' the spatial harmonics used in the restricted spatial regression model. The value must be between 1 and N and implies that
#' the eigen vectors associated with the largest \code{moan.cut} eigen values are used for the basis functions.
#' The item \code{knots} are xy locations of the discrete process
#' convolution knots. 
#' @param so.data An \code{so.data} object containing the observed occupancies,
#' detection covariates, site covariates, and site coordinates. This is created
#' via the \code{\link{make.so.data}}
#' @param prior A named list that provides the parameter values for the prior
#' distributions. At the current time the elements of the list must contain
#' \code{a.tau} and \code{b.tau} which are the parameters for the gamma prior on the spatial
#' process parameter in the occupancy model. Other elements may include
#' \code{Q.b} and \code{mu.b} which are the tolerance and mean for the beta
#' vector (detection parameters). Also \code{Q.g} and \code{mu.g} which are the
#' prior parameters for the occupancy model. If the Q.b and Q.g are left out,
#' the default is Q.b = 0 and Q.g = 0 (i.e., flat priors). If mu.b and mu.g are left out,
#' the default is zero vectors.
#' @param control A named list with the control parameters for the MCMC. The
#' elements of the list must include: (1) \code{burnin} is the number of
#' iterations of burnin, (2) \code{iter} is the total number of iterations
#' retained for the MCMC sample, and (3) \code{thin} is the thining rate of the
#' chain. The real number of MCMC iterations is equal to \code{iter*thin} of
#' which \code{iter - burnin} are retained for posterior summary.
#' @param initial.values A named list that can include any or all of the following vectors or scalers 
#' (1) \code{beta}, a vector of initial values for the detection parameters, (2) \code{gamma}, a vector or
#' initial values for the occupancy model, and (3) \code{tau}, an initial value for the spatial precision 
#' parameter.
#' @return A list with the following elements: \item{beta}{An object of class
#' \code{mcmc}. The detection model parameters.} \item{gamma}{An object of
#' class \code{mcmc}. The occupancy model parameters.} \item{psi}{An object of
#' class \code{mcmc}. The occupancy generating process} \item{real.occ}{An
#' object of class \code{mcmc}. The realized occupancy at the time of the
#' survey}  \item{tau}{An object of
#' class \code{mcmc}. The variance parameter for the spatial model}
#' \item{occupancy.df}{A data frame with the spatial
#' coordinates, site id, and posterior mean and variance of psi, eta, and real.occ}
#' \item{D.m}{The posterior predictive loss criterion of Gelfand and Ghosh
#' (1998; Biometrika 85:1-11) for model selection. The criterion is a combination of a goodness-of-fit measure, G.m, and 
#' a complexity measure, P.m, similar information criteria such as AIC and BIC. D.m = G.m + P.m. Lower values of D.m imply lower 
#' expected loss in predicting new data with the posterior model parameters.} 
#' \item{G.m}{The
#' goodness-of-fit portion of D.m} \item{P.m}{The model complexity component of
#' D.m} \item{detection.model}{The detection model call.}
#' \item{occupancy.model}{The occupancy model call.} \item{model}{A character
#' version of the joint occupancy and detection model call. This is useful for
#' saving results.}
#' @author Devin S. Johnson <devin.johnson@@noaa.gov>
#' @export
#' @import truncnorm
#' @import coda
#' @import Matrix
#' @import fields
#' @import rARPACK
spatial.occupancy <-
  function(detection.model, occupancy.model, spatial.model, so.data, prior, control, initial.values=NULL){
    #Packages
    
    #Matrix construction, etc...
    site <- so.data$site
    visit <- so.data$visit
    site$site.idx <- factor(site[,attr(site,"site")])
    visit$site.idx <- factor(visit[,attr(visit,"site")], levels=levels(site$site.idx))
    xy <- as.matrix(site[,attr(site,"coords")])
    Xy <- Matrix(model.matrix(detection.model, visit))
    Xz <- Matrix(model.matrix(occupancy.model, site))
    n.site <- nrow(Xz)
    n.knots <- nrow(spatial.model$knots)
    M  <- Matrix(model.matrix(~site.idx-1, visit))	
    
    z.master=as.vector(ifelse(table(visit$site.idx,visit[,attr(visit,"obs")])[,"1"]>0, 1, NA))
    z <- as.numeric(!is.na(z.master))
    z.obs <- site$site.idx %in% visit$site.idx
    n.obs <- sum(z.obs)
    
    y <- visit[,attr(visit,"obs")]
    
    ## spatial model
    if(spatial.model$model=="icar"){
      Q.eta <- Matrix(icar.Q(xy,spatial.model$threshold, 1))
    }
    #	else if(spatial.model$model=="car"){
    #		Q.fun <- icar.Q(xy,spatial.model$threshold, fun=TRUE)
    #		rho<- 0.95
    #		Q.eta <- Matrix(Q.fun(rho))
    #		ln.det.Q.eta <- determinant(Q.eta, log=TRUE)$modulus
    #		alpha <- rnorm(n.site,0,0.01)
    #		#alpha <- alpha-mean(alpha)
    #		sigma <- 1 
    #		eta <- sigma*alpha
    #	}
    else if(spatial.model$model=="proc.conv"){
      Q.alpha <- Diagonal(nrow(spatial.model$knots))
      if(is.null(spatial.model$pc.scale)) pc.scale <- min(dist(knots))
      else pc.scale <- spatial.model$pc.scale
      K <- Matrix(exp(-0.5*(rdist(xy,spatial.model$knots)/pc.scale)^2))
      KtK <- crossprod(K)
      a <- rep(0,n.knots)
    }
    else if(spatial.model$model=="rsr"){
      cat("\nCreating (R)estricted (S)patial (R)egression matrices ...\n")
      if(as.integer(spatial.model$moran.cut)<1  | as.integer(spatial.model$moran.cut)>n.site){
        stop("Invalid value for 'moran.cut' specified. See documentation\n")
      }
      Q <- Matrix(icar.Q(xy,spatial.model$threshold, rho=1))
      A <- Matrix(diag(diag(Q)) - Q)
      P <- diag(n.site) - Xz %*% solve(crossprod(Xz), t(Xz))
      Op <- (nrow(A)/sum(A))*(P %*% (A %*% P))
      e <- rARPACK::eigs(Op, as.integer(spatial.model$moran.cut))
      K <- e$vectors
      KtK <- diag(ncol(K)) 
      Q.alpha <- as.matrix(t(K) %*% Q %*% K)
      a <- rep(0,nrow(Q.alpha))
    }
    
    ## Priors
    a.tau <- prior$a.tau
    b.tau <- prior$b.tau
    
    if(!is.null(prior$mu.b)) mu.b <- prior$mu.b
    else mu.b <- rep(0,ncol(Xy))
    
    Q.b = prior$Q.b
    if(is.null(Q.b)) {
      Q.b = matrix(0, ncol(Xy), ncol(Xy))
    } else if(is.vector(Q.b)){
      Q.b=diag(Q.b, ncol(Xy))
    } 
    
    if(!is.null(prior$mu.g)) mu.g <- prior$mu.g
    else mu.g <- rep(0,ncol(Xz))
    
    Q.g = prior$Q.g
    if(is.null(Q.g)) {
      Q.g = matrix(0, ncol(Xz), ncol(Xz))
    } else if(is.vector(Q.g)){
      Q.g=diag(Q.g, ncol(Xz))
    } 
    
    
    # Storage
    beta <- matrix(nrow=control$iter-control$burnin,ncol=ncol(Xy))
    colnames(beta) <- colnames(Xy)
    gamma <- matrix(nrow=control$iter-control$burnin,ncol=ncol(Xz))
    colnames(gamma) <- colnames(Xz)
    psi <- matrix(nrow=control$iter-control$burnin,ncol=n.site)
    eta.stor <- matrix(nrow=control$iter-control$burnin,ncol=n.site)
    real.occ <- matrix(nrow=control$iter-control$burnin,ncol=n.site)
    tau.stor <- matrix(nrow=control$iter-control$burnin,ncol=1)
    colnames(tau.stor) <- "tau"
    if(spatial.model$model %in% c("proc.conv", "rsr")){
      alpha <- matrix(nrow=control$iter-control$burnin,ncol=nrow(Q.alpha))
    }
    else alpha <- NULL
    #if(spatial.model$model=="car") rho.eta <- matrix(nrow=control$iter-control$burnin,ncol=1)
    #yz.ln.lik <- matrix(nrow=control$iter-control$burnin,ncol=1)
    y.rep.stor <- matrix(nrow=control$iter-control$burnin,ncol=nrow(Xy))
    
    V.g.inv <- crossprod(Xz) + Q.g
    I.n <- Matrix(diag(n.site))
    
    # Starting values
    if(is.null(initial.values$beta)) {
      b <- rep(0,ncol(Xy))
    } else b <- initial.values$beta
    if(length(b)!=ncol(Xy)) {
      stop("Error: Length of initial values for beta does not match the detection model!\n")
    }
    if(is.null(initial.values$gamma)) {
      g <- rep(0,ncol(Xz))
    } else g <- initial.values$gamma
    if(length(g)!=ncol(Xz)) stop("Error: Length of initial values for gamma does not match the occupancy model!\n")
    if(is.null(initial.values$tau)) {
      tau <- 1
    } else tau <- initial.values$tau 
    eta <- rep(0,n.site)
    
    
    iter <- control$iter*control$thin
    burnin <- control$burnin*control$thin
    
    cat("\nBeginning MCMC routine ...\n")
    st <- Sys.time()
    for(i in 1:iter){
      Xge <- as.numeric(Xz%*%g+eta)
      Xb <- as.numeric(Xy%*%b)
      
      #Update missing z
      ln.qz <- pnorm(rep(0,n.obs), Xge, 1, lower.tail=FALSE, log.p=TRUE)
      ln.qy <- crossprod(M,pnorm(0, Xb, 1, log.p=TRUE))
      p.z <- exp(ln.qz + ln.qy)/(exp(ln.qz + ln.qy) + pnorm(rep(0,n.obs), Xge, 1))
      z[is.na(z.master)] <- rbinom(sum(is.na(z.master)), 1, p.z[is.na(z.master)])
      
      #Update y.tilde
      y.tilde <- rtruncnorm(nrow(Xy), a=ifelse(y==0,-Inf,0), b=ifelse(y==0,0,Inf), Xb, 1)
      
      #Update z.tilde
      z.tilde <- rtruncnorm(n.site, a=ifelse(z==0,-Inf,0), b=ifelse(z==0,0,Inf), Xge, 1)
      
      #Update b
      idx <- visit$site.idx %in% (site$site.idx[z==1])
      #Q.b <- crossprod(Xy[idx,])
      V.b.inv <- crossprod(Xy[idx,],Xy[idx,]) + Q.b
      m.b <- solve(V.b.inv, crossprod(Xy[idx,],y.tilde[idx])+crossprod(Q.b,mu.b))
      b <- m.b + solve(chol(V.b.inv), rnorm(ncol(Xy),0,1))
      
      #Update g
      m.g <- solve(V.g.inv, crossprod(Xz,z.tilde-eta)+crossprod(Q.g,mu.g))
      g <- m.g + solve(chol(V.g.inv), rnorm(ncol(Xz),0,1))
      
      #Update eta
      if(spatial.model$model %in% c("car","icar")){
        V.eta.inv <- I.n + tau*Q.eta
        m.eta <- solve(V.eta.inv, z.tilde-Xz%*%g)
        eta <- m.eta + solve(chol(V.eta.inv), rnorm(n.site,0,1), tol=1e-10)
        if(spatial.model$model=="icar") eta <- eta - mean(eta)
      }
      else if(spatial.model$model %in% c("proc.conv","rsr")){
        V.alpha.inv <- KtK + tau*Q.alpha
        m.alpha <- solve(V.alpha.inv, crossprod(K,z.tilde-Xz%*%g))
        #cat("\n", min((diag(chol(V.alpha.inv)))),"           ",tau,"    ",length(tau),"\n")
        a <- m.alpha + backsolve(chol(V.alpha.inv), rnorm(nrow(Q.alpha),0,1))
        eta <- K%*%a
      }
      else NULL
      
      #Update tau
      if(spatial.model$model=="icar") tau <- rgamma(1, (n.site-1)/2 + a.tau, as.numeric(crossprod(eta, Q.eta %*% eta)/2) + b.tau)
      else if(spatial.model$model %in% c("proc.conv","rsr")) tau <- rgamma(1, length(a)/2 + a.tau, as.numeric(crossprod(a, Q.alpha %*% a)/2) + b.tau)
      else NULL 
      
      #		#Update rho for "car" models
      #		if(spatial.model$model=="car"){
      #			rho.p <- rbeta(1,a.rho.tune,1)
      #			Q.eta.p <- Q.fun(rho.p)
      #			ln.det.Q.eta.p <- determinant(Q.eta.p)$modulus
      #			mh <- exp(
      #						(ln.det.Q.eta.p - crossprod(alpha, crossprod(Q.eta.p, alpha)))/2 + dbeta(rho.p, a.rho, b.rho, log=TRUE)
      #						- dbeta(rho.p, a.rho.tune, 1, log=TRUE)
      #						-  (ln.det.Q.eta.p - crossprod(alpha, crossprod(Q.eta.p, alpha)))/2 + dbeta(rho.p, a.rho, b.rho, log=TRUE)
      #						+ dbeta(rho, a.rho.tune, 1, log=TRUE)
      #						)
      #			if(runif(1,0,1)<=as.numeric(mh)){
      #				rho <- rho.p
      #				Q.eta <- Q.eta.p
      #				ln.det.Q.eta <- ln.det.Q.eta.p
      #			}	
      #		}
      
      #Record sample
      if(i>burnin & i%%control$thin==0){
        beta[(i-burnin)/control$thin,] <- as.numeric(b)
        gamma[(i-burnin)/control$thin,] <- as.numeric(g)
        psi[(i-burnin)/control$thin,] <- pnorm(rep(0,n.site), as.numeric(Xz%*%g+eta), 1, lower.tail=FALSE)
        eta.stor[(i-burnin)/control$thin,] <- as.numeric(eta)
        real.occ[(i-burnin)/control$thin,] <- as.numeric(z)
        tau.stor[(i-burnin)/control$thin] <- as.numeric(tau)
        if(spatial.model$model %in% c("proc.conv","rsr")){
          alpha[(i-burnin)/control$thin,] <- as.numeric(a)
        }
        #if(spatial.model$model=="car") rho.eta[(i-burnin)/control$thin] <- rho
        #ln.p.ygz <- crossprod(M, dbinom(y,1, pnorm(0, as.vector(Xy%*%b), 1, lower.tail=FALSE), log=TRUE))[z.obs]
        #ln.p.z <- pnorm(rep(0,sum(z.obs)), as.numeric(Xz%*%g+eta)[z.obs], 1, lower.tail=FALSE, log.p=TRUE)
        #yz.ln.lik[(i-burnin)/control$thin] <- sum(log(exp(ln.p.ygz + ln.p.z) + ifelse(is.na(z.master)[z.obs], 1-exp(ln.p.z), 0)))
        ##Replicated data
        idx.rep <- visit$site.idx %in% site$site.idx[rnorm(n.site, as.numeric(Xz%*%g + eta), sd=1)<=0]
        y.rep <- 1.0*(rnorm(nrow(Xy), as.numeric(Xy%*%b), sd=1)>0)
        y.rep[idx.rep] <- 0
        y.rep.stor[(i-burnin)/control$thin,] <- y.rep
      }
      if(i==15){
        tpi <- as.numeric(difftime(Sys.time(), st, units="secs"))/15
        ttc <- round((iter-15)*tpi/3600, 2)
        cat("\nApproximate time till completion: ", ttc, " hours\n")
      }
      if(100*(i/iter) >= 10 & (100*(i/iter))%%10==0) cat("\n", 100*(i/iter), "% completed\n")
    }	
    
    beta <- mcmc(beta)
    gamma <- mcmc(gamma)
    psi <- mcmc(psi)
    real.occ <- mcmc(real.occ)
    tau <- mcmc(tau.stor)
    #cross.eta <- mcmc(cross.eta)
    #if(spatial.model$model=="car") rho.eta <- mcmc(rho.eta)
    #else rho.eta <- NULL 
    #ln.prior.pred <- mean(yz.ln.lik)-log(mean(exp(-(yz.ln.lik-mean(yz.ln.lik)))))
    #l.max <- mean(yz.ln.lik) + var(yz.ln.lik)
    #d.hat <- 2*var(yz.ln.lik)
    #BICM <- -2*l.max + ncol(Xy)*log(sum(y==0)+1) + ncol(Xz)*log(sum(!is.na(z.master))+1) + (d.hat-ncol(Xy)-ncol(Xz))*log(sum(samp)+1)
    #AICM <- -2*l.max + 2*d.hat
    G.m <- crossprod(apply(y.rep.stor, 2, mean)-y)
    P.m <- crossprod(apply(y.rep.stor, 2, sd))
    D.m <- G.m + P.m
    occupancy.df <- data.frame(
      site[attr(site,"site")], 
      site[,attr(site,"coords")], 
      samp=as.numeric(site$site.idx %in% unique(visit$site.idx)), 
      psi.est=apply(psi,2,mean), 
      psi.se=apply(psi,2,sd), 
      eta.est=apply(eta.stor,2,mean),
      eta.se=apply(eta.stor,2,sd),
      real.occ.est=apply(real.occ,2,mean), 
      real.occ.se=apply(real.occ, 2, sd)
    )
    out <- list(
      beta=beta,
      gamma=gamma,
      psi=psi, 
      real.occ=real.occ, 
      tau=tau, 
      alpha=alpha, 
      occupancy.df=occupancy.df, 
      #ln.prior.pred=ln.prior.pred, yz.ln.lik=yz.ln.lik, BICM=BICM, AICM=AICM, 
      D.m=D.m, 
      G.m=G.m, 
      P.m=P.m,
      detection.model=detection.model, 
      occupancy.model=occupancy.model, 
      model=paste(
        c(
          paste(c("DET", as.character(detection.model)[-1]),collapse="~"),
          paste(c("OCC", as.character(occupancy.model)[-1]),collapse="~")
        ),
        collapse='-'
      ),
      so.data=so.data
    )
    
    class(out) <- "spat.occ"
    return(out)
  }

