crmh = function(a,x,y,w,s) {  ## posterior
  v = exp(-a^2/2/s^2) 
  for (i in 1:length(x)) {
	v = v * ((x[i]^exp(a))^y[i])*(((1-w[i]*x[i]^exp(a))^(1-y[i])))
  }
  return(v)
}
crmht = function(a,x,y,w,s) { ## posterior times x
  v = a * exp(-a^2/2/s^2) 
  for (i in 1:length(x)) {
	v = v * ((x[i]^exp(a))^y[i])*(((1-w[i]*x[i]^exp(a))^(1-y[i])))
  }
  return(v)
}
crmht2 = function(a,x,y,w,s) { ## posterior times x^2
  v = a^2 * exp(-a^2/2/s^2) 
  for (i in 1:length(x)) {
	v = v * ((x[i]^exp(a))^y[i])*(((1-w[i]*x[i]^exp(a))^(1-y[i])))
  }
  return(v)
}

  
crmhlgt <- function(a,x,y,w,s,alp0)  { ## posterior logit model
  v = exp(-a^2/2/s^2) 
  for (i in 1:length(x)) {
 	PSI <- (1 + exp(-alp0-exp(a)*x[i]))^{-1}
    	v <- v * (PSI^y[i]) * (1-w[i]*PSI)^(1-y[i])
  }
  return(v)
}
crmhtlgt <- function(a,x,y,w,s,alp0)  { ## posterior times x
  v = a * exp(-a^2/2/s^2) 
  for (i in 1:length(x)) {
 	PSI <- (1 + exp(-alp0-exp(a)*x[i]))^{-1}
    	v <- v * (PSI^y[i]) * (1-w[i]*PSI)^(1-y[i])
  }
  return(v)
}
crmht2lgt <- function(a,x,y,w,s,alp0)  { ## posterior times x^2
  v = a^2 * exp(-a^2/2/s^2) 
  for (i in 1:length(x)) {
 	PSI <- (1 + exp(-alp0-exp(a)*x[i]))^{-1}
    	v <- v * (PSI^y[i]) * (1-w[i]*PSI)^(1-y[i])
  }
  return(v)
}
 


lcrm <- function(a,x,y,w) { #loglikelihood of empiric function
  v <- 0
  for (i in 1:length(x))
    v <- v + y[i]*log(x[i])*exp(a) + (1-y[i])*log(1 - w[i]*x[i]^exp(a))
  return(v)
}
lcrmlgt <- function(a,x,y,w,alp0) { #loglikelihood of logit function
  v <- 0
  for (i in 1:length(x)) {
    PSI <- (1 + exp(-alp0-exp(a)*x[i]))^{-1}
    v <- v + y[i]*log(PSI) + (1-y[i])*log(1-w[i]*PSI)
  }
  return(v)
}

crm <- function(prior, target, tox, level, n=length(level),
                dosename=NULL, include=1:n, pid=1:n, conf.level=0.90,
                method="bayes", model="empiric", intcpt=3,
                scale=sqrt(1.34), model.detail=TRUE, patient.detail=TRUE, var.est=TRUE) {

  y1p <- tox[include]
  w1p <- rep(1,length(include))
  if (model=="empiric") {
    dosescaled <- prior
    
    x1p <- prior[level[include]]
    if (method=="mle") {
      if (sum(y1p)==0 | sum(y1p)==length(y1p)) stop(" mle does not exist!")
      est <- optimize(lcrm,c(-10,10),x1p,y1p,w1p,tol=0.0001,maximum=TRUE)$max
      if (var.est) { e2 <- integrate(crmht2,-100,100,x1p,y1p,w1p,500,abs.tol=0)[[1]] / integrate(crmh,-10,10,x1p,y1p,w1p,500,abs.tol=0)[[1]]; }
    }
    else if (method=="bayes") {
      den <- integrate(crmh,-Inf,Inf,x1p,y1p,w1p,scale,abs.tol=0)[[1]]
      est <- integrate(crmht,-10,10,x1p,y1p,w1p,scale,abs.tol=0)[[1]] / den
      if (var.est) { e2 <- integrate(crmht2,-10,10,x1p,y1p,w1p,scale,abs.tol=0)[[1]] / den; }
    }
    else { stop(" unknown estimation method"); }
    ptox <- prior^exp(est)
    if (var.est) {
      post.var <- e2-est^2
      crit <- qnorm(0.5+conf.level/2)
      lb <- est - crit*sqrt(post.var)
      ub <- est + crit*sqrt(post.var)
      ptoxL <- prior^exp(ub)
      ptoxU <- prior^exp(lb)
    }
  }
  else if (model=="logistic") {
    dosescaled <- log(prior/(1-prior)) - intcpt
    if (!all(dosescaled<0)) {
      stop( "Intercept parameter in logit model is too small: scaled doses > 0!")
    }
   
    x1p <- dosescaled[level[include]]
 
    if (method=="mle") {
      if (sum(y1p)==0 | sum(y1p)==length(y1p)) stop(" mle does not exist!")
      est <- optimize(lcrmlgt,c(-10,10),x1p,y1p,w1p,intcpt,tol=0.0001,maximum=TRUE)$max
      if (var.est) { e2 <- integrate(crmht2lgt,-100,100,x1p,y1p,w1p,500,intcpt,abs.tol=0)[[1]] / integrate(crmhlgt,-10,10,x1p,y1p,w1p,500,intcpt,abs.tol=0)[[1]]; }
    }
    else if (method=="bayes") {
      den <- integrate(crmhlgt,-Inf,Inf,x1p,y1p,w1p,scale,intcpt,abs.tol=0)[[1]]
      est <- integrate(crmhtlgt,-10,10,x1p,y1p,w1p,scale,intcpt,abs.tol=0)[[1]] / den
      if (var.est) { e2 <- integrate(crmht2lgt,-10,10,x1p,y1p,w1p,scale,intcpt,abs.tol=0)[[1]] / den; }
    }
    else { stop(" unknown estimation method"); }
    ptox <- (1 + exp(-intcpt-exp(est)*dosescaled))^{-1}
    if (var.est) {
      post.var <- e2-est^2
      crit <- qnorm(0.5+conf.level/2)
      lb <- est - crit*sqrt(post.var)
      ub <- est + crit*sqrt(post.var)
      ptoxL <- (1 + exp(-intcpt-exp(ub)*dosescaled))^{-1}
      ptoxU <- (1 + exp(-intcpt-exp(lb)*dosescaled))^{-1}
    }
  }
  else { stop(" model specified not available."); }

  if (all(ptox<=target)) { rec <- length(prior); }
  else if (all(ptox>=target)) { rec <- 1; }
  else { rec <- order(abs(ptox-target))[1]; }
  if (!var.est) { post.var <- ptoxU <- ptoxL <- NA; }
  foo <- list(prior=prior, target=target, tox=tox, level=level,
              dosename=dosename, subset=pid[include], estimate=est,
              model=model, prior.var=scale^2, post.var=post.var,method=method,
              mtd=rec, include=include, pid=pid, model.detail=model.detail,intcpt=intcpt,
              ptox=ptox, ptoxL=ptoxL, ptoxU=ptoxU, conf.level=conf.level,
              patient.detail=patient.detail,tite=FALSE,dosescaled=dosescaled,var.est=var.est)
  class(foo) <- "mtd"
  foo
}

crmsim <- function(PI, prior, target, n, x0, nsim=1, mcohort=1, restrict=TRUE, count=TRUE,
                    method="bayes", model="empiric", intcpt=3, scale=sqrt(1.34),seed=1009) {
  if (nsim==1) {
    foo <- onetrial(PI,prior,target,n,x0,mcohort=mcohort,restrict=restrict,method=method,
                    model=model,intcpt=intcpt,scale=scale,seed=seed)
  }
  else {
    foo <- mtrials(nsim,PI,prior,target,n,x0,mcohort=mcohort,restrict=restrict,method=method,
                   model=model,intcpt=intcpt,scale=scale,seed=seed,count=count)
  }
  foo
}

mtrials <- function(nsim, PI, prior, target, n, x0, mcohort=1, restrict=TRUE, count=TRUE,
                    method="bayes", model="empiric", intcpt=3, scale=sqrt(1.34),seed=1009) {
  set.seed(seed)
  if (nsim<=1) { stop(" nsim less than 2!  For single simulation replicate, use function `onetrial'."); }
  if ((n/mcohort != round(n/mcohort))) { stop(" Number of patients is not multiple of cohort size!"); }
  M <- n/mcohort

  nexpt <- ntox <- sel <- rep(0,length(prior))
  BETAHAT <- matrix(rep(NA,nsim*n),nrow=nsim)
  final.est <- rep(NA,nsim)
  
  for (r in 1:nsim) {
    if (count) { cat("simulation number:",r,"\n"); }
                     
    bethat <- y <- level <- rep(NA,n)
    if (length(x0)>1) {  # length(x0)>1 implies it's a user-input two-stage design
      if (length(x0)!=n) { stop(" Initial design has a different sample size than that is specified!"); }
      
      stage1 <- TRUE
      for (m in 1:M) {
        pos <- ((m-1)*mcohort+1):(m*mcohort)
        
        if (stage1) {
          bethat[pos] <- 0
          level[pos] <- x0[pos]
          ynew <- rbinom(mcohort, 1, PI[level[pos]])
          y[pos] <- ynew
          
          xcur <- level[1:pos[mcohort]]
          ycur <- y[1:pos[mcohort]]
          
          if (any(ynew==1)) { stage1 <- FALSE; }
        }
        else {
          if (mean(ycur)==1) {
            obj <- crm(prior, target,ycur,xcur,model=model,intcpt=intcpt,scale=500,var.est=FALSE)
          }
          else {
            obj <- crm(prior, target, ycur, xcur, method=method, model=model, intcpt=intcpt, scale=scale, var.est=FALSE)
          }
          cur <- obj$mtd
          bethat[pos] <- obj$est
          if (restrict) {
            if (sum(ynew)/mcohort >= target) { cur <- min(cur, xcur[length(xcur)]); } # coherent restriction
            else { cur <- min(cur, xcur[length(xcur)]+1); } # no skipping dose restriction
          }
          
          level[pos] <- cur
          ynew <- rbinom(mcohort, 1, PI[level[pos]])
          y[pos] <- ynew
          
          xcur <- level[1:pos[mcohort]]
          ycur <- y[1:pos[mcohort]]
        }
      }
      BETAHAT[r,] <- bethat
    }
    else {
      if (method=="mle") { stop(" Require an initial design for mle-CRM!"); }
      pos <- 1:mcohort
      bethat[pos] <- 0
      level[pos] <- x0    
      ynew <- rbinom(mcohort, 1, PI[level[pos]])
      y[pos] <- ynew
      
      xcur <- level[1:pos[mcohort]]
      ycur <- y[1:pos[mcohort]]
      for (m in 2:M) {
        pos <- ((m-1)*mcohort+1):(m*mcohort)
        obj <- crm(prior, target, ycur, xcur, method=method, model=model, intcpt=intcpt, scale=scale, var.est=FALSE)
        cur <- obj$mtd
        bethat[pos] <- obj$est
        if (restrict) {
          if (sum(ynew)/mcohort >= target) { cur <- min(cur, xcur[length(xcur)]); }
          else { cur <- min(cur, xcur[length(xcur)]+1); }
        }
        
        level[pos] <- cur
        ynew <- rbinom(mcohort, 1, PI[level[pos]])
        y[pos] <- ynew
        
        xcur <- level[1:pos[mcohort]]
        ycur <- y[1:pos[mcohort]]
      }
      BETAHAT[r,] <- bethat
    }

    if (method=="mle" & (sum(ycur)==0|mean(ycur)==1)) {
      finalobj <- crm(prior, target, ycur, xcur, model=model, intcpt=intcpt,scale=500, var.est=FALSE)
      msg <- "Warning: mle is approximated"
    }
    else {
      finalobj <- crm(prior, target, ycur, xcur, method=method,model=model,intcpt=intcpt,scale=scale, var.est=FALSE)
      msg <- "Okay"
    }
    rec <- finalobj$mtd
    sel[rec] <- sel[rec] + 1
    for (k in 1:length(prior)) {
      nexpt[k] <- nexpt[k] + length(which(xcur==k))
      ntox[k] <- ntox[k] + length(which(xcur==k & ycur==1))
    }
    final.est[r] <- finalobj$est
  }
  sel <- sel/nsim
  nexpt <- nexpt/nsim
  ntox <- ntox/nsim

  if (length(x0)==1) {
    if (mcohort==1) { design <- paste("CRM starting at dose",x0); }
    else { design <- paste("Group CRM starting at dose",x0,"with group size",mcohort); }
  }
  else {
    if (mcohort==1) { design <- "Two-stage CRM"; }
    else { design <- paste("Two-stage group CRM with group size",mcohort); }
  }
      
  foo <- list(PI=PI, prior=prior, target=target, n=n, x0=x0, nsim=nsim,
              MTD=sel, level=nexpt, tox=ntox, beta.hat=BETAHAT, final.est=final.est, 
              design=design, method=method, prior.var=scale^2, model=model, intcpt=intcpt,
              mcohort=mcohort, restriction=restrict,seed=seed,tite=FALSE,dosescaled=finalobj$dosescaled,msg=msg)
  class(foo) <- "sim"
  foo
}

onetrial <- function(PI, prior, target, n, x0, mcohort=1, restrict=TRUE,
                     method="bayes", model="empiric", intcpt=3, scale=sqrt(1.34),seed=1009) {
  set.seed(seed)
  if ((n/mcohort != round(n/mcohort))) { stop(" Number of patients is not multiple of cohort size!"); }
  M <- n/mcohort
  bethat <- y <- level <- rep(NA,n)

  if (length(x0)>1) {  # length(x0)>1 implies it's a user-input two-stage design
    if (length(x0)!=n) { stop(" Initial design has a different sample size than that is specified!"); }

    stage1 <- TRUE
    for (m in 1:M) {
      pos <- ((m-1)*mcohort+1):(m*mcohort)

      if (stage1) {
        bethat[pos] <- 0
        level[pos] <- x0[pos]
        ynew <- rbinom(mcohort, 1, PI[level[pos]])
        y[pos] <- ynew
        
        xcur <- level[1:pos[mcohort]]
        ycur <- y[1:pos[mcohort]]
        
        if (any(ynew==1)) { stage1 <- FALSE; }
      }
      else {
        if (mean(ycur)==1) {
          obj <- crm(prior,target,ycur,xcur,model=model,intcpt=intcpt,scale=500,var.est=FALSE)
        }
        else {
          obj <- crm(prior, target, ycur, xcur, method=method, model=model, intcpt=intcpt, scale=scale, var.est=FALSE)
        }
        
        cur <- obj$mtd
        bethat[pos] <- obj$est
        if (restrict) {
          if (sum(ynew)/mcohort >= target) { cur <- min(cur, xcur[length(xcur)]); } # coherent restriction
          else { cur <- min(cur, xcur[length(xcur)]+1); } # no skipping dose restriction
        }

        level[pos] <- cur
        ynew <- rbinom(mcohort, 1, PI[level[pos]])
        y[pos] <- ynew

        xcur <- level[1:pos[mcohort]]
        ycur <- y[1:pos[mcohort]]
      }
    }
  }
  else {
    if (method=="mle") { stop(" Require an initial design for mle-CRM!"); }
    pos <- 1:mcohort
    bethat[pos] <- 0
    level[pos] <- x0    
    ynew <- rbinom(mcohort, 1, PI[level[pos]])
    y[pos] <- ynew
    
    xcur <- level[1:pos[mcohort]]
    ycur <- y[1:pos[mcohort]]
    for (m in 2:M) {
      pos <- ((m-1)*mcohort+1):(m*mcohort)
      obj <- crm(prior, target, ycur, xcur, method=method, model=model, intcpt=intcpt, scale=scale, var.est=FALSE)
      cur <- obj$mtd
      bethat[pos] <- obj$est
      if (restrict) {
        if (sum(ynew)/mcohort >= target) { cur <- min(cur, xcur[length(xcur)]); }
        else { cur <- min(cur, xcur[length(xcur)]+1); }
      }

      level[pos] <- cur
      ynew <- rbinom(mcohort, 1, PI[level[pos]])
      y[pos] <- ynew

      xcur <- level[1:pos[mcohort]]
      ycur <- y[1:pos[mcohort]]
    }
  }

  if (method=="mle" & (sum(ycur)==0|mean(ycur)==1)) {
    finalobj <- crm(prior,target,ycur,xcur,model=model,intcpt=intcpt,scale=500)
    msg <- "Warning: mle is approximated"
  }
  else {
    finalobj <- crm(prior, target, ycur, xcur, method=method,model=model,intcpt=intcpt,scale=scale)
    msg <- "Okay"
  }
  
  cur <- finalobj$mtd; est <- finalobj$est;
  if (length(x0)==1) {
    if (mcohort==1) { design <- paste("CRM starting at dose",x0); }
    else { design <- paste("Group CRM starting at dose",x0,"with group size",mcohort); }
  }
  else {
    if (mcohort==1) { design <- "Two-stage CRM"; }
    else { design <- paste("Two-stage group CRM with group size",mcohort); }
  }
  
  foo <- list(PI=PI, prior=prior, target=target, n=n, x0=x0, nsim=1, 
              MTD=cur, level=level, tox=y, beta.hat=bethat, final.est=est,
              design=design, method=method, prior.var=scale^2, model=model, intcpt=intcpt,
              mcohort=mcohort, restriction=restrict, seed=seed, tite=FALSE, dosescaled=finalobj$dosescaled, msg=msg,
              post.var=finalobj$post.var, ptox=finalobj$ptox, ptoxL=finalobj$ptoxL, ptoxU=finalobj$ptoxU,
              conf.level=finalobj$conf.level)
  class(foo) <- "sim"
  foo
}

titecrm <- function(prior, target, tox, level, n=length(level),
                    weights=NULL, followup=NULL, entry=NULL, exit=NULL,
                    obswin=NULL, scheme="linear", conf.level=0.90,
                    dosename=NULL, include=1:n, pid=1:n, method="bayes",model="empiric",var.est=TRUE,
                    scale=sqrt(1.34), intcpt=3, model.detail=TRUE, patient.detail=TRUE, tite=TRUE) {
  if (is.null(weights)) {
    if (is.null(followup)) { followup <- exit-entry; }
    if (scheme=="linear") { weights <- followup/obswin; }
    else if (scheme=="adaptive") {
      support <- sort(followup[tox==1])
      z <- length(support)
      if (z) {
        for (i in 1:n) {
          m <- length(support[support<=followup[i]])
          if (!m) weights[i] <- followup[i] / support[1] / (z+1)
          else if (m==z) weights[i] <- (z + (followup[i]-support[z])/(obswin-support[z]))/(z+1)
          else weights[i] <- (m + (followup[i]-support[m])/(support[m+1]-support[m]))/(z+1)
        }
      }
      else { weights <- followup/obswin; }
    }
    else { stop(" Weighting scheme undefined!"); }
    weights <- pmin(weights, 1)
  }
  
  if (any(weights>1) | any(weights<0)) stop(" Weights have to be between 0 and 1!")
  
  if (is.null(pid)) {
    if (! (length(tox)==length(level) & length(tox)==length(weights)))
      stop(" tox, level, and weights are of different lengths!")
  }
  else {
    if (! (length(tox)==length(level) & length(tox)==length(weights) & length(tox)==length(pid)) )
      stop(" pid, tox, level, and weights are of different lengths!")
  }
  weights[tox==1] <- 1
  
  y1p <- tox[include]
  w1p <- weights[include]
  if (model=="empiric") {
    dosescaled <- prior
    
    x1p <- prior[level[include]]
    if (method=="mle") {
      if (sum(y1p)==0 | sum(y1p)==length(y1p)) stop(" mle does not exist!")
      est <- optimize(lcrm,c(-10,10),x1p,y1p,w1p,tol=0.0001,maximum=TRUE)$max
      if (var.est) { e2 <- integrate(crmht2,-10,10,x1p,y1p,w1p,500,abs.tol=0)[[1]] / integrate(crmh,-10,10,x1p,y1p,w1p,500,abs.tol=0)[[1]]; }
    }
    else if (method=="bayes") {
      den <- integrate(crmh,-Inf,Inf,x1p,y1p,w1p,scale,abs.tol=0)[[1]]
      est <- integrate(crmht,-10,10,x1p,y1p,w1p,scale,abs.tol=0)[[1]] / den
      if (var.est) { e2 <- integrate(crmht2,-10,10,x1p,y1p,w1p,scale,abs.tol=0)[[1]] / den; }
    }
    else { stop(" unknown estimation method"); }
    ptox <- prior^exp(est)
    if (var.est) {
      post.var <- e2-est^2
      crit <- qnorm(0.5+conf.level/2)
      lb <- est - crit*sqrt(post.var)
      ub <- est + crit*sqrt(post.var)
      ptoxL <- prior^exp(ub)
      ptoxU <- prior^exp(lb)
    }
  }
  else if (model=="logistic") {
    dosescaled <- log(prior/(1-prior)) - intcpt
    if (!all(dosescaled<0)) {
      stop( "intercept parameter in logit model is too small: scaled doses > 0!")
    }
#    LB <- log(  (log((1+target)/(1-target)) - intcpt)/dosescaled[1] ) - 3
#    UB <- log(  (log((target/2)/(1-target/2)) - intcpt)/dosescaled[length(prior)] ) + 3
    
    x1p <- dosescaled[level[include]]

    if (method=="mle") {
      if (sum(y1p)==0 | sum(y1p)==length(y1p)) stop(" mle does not exist!")
      est <- optimize(lcrmlgt,c(-10,10),x1p,y1p,w1p,intcpt,tol=0.0001,maximum=TRUE)$max
      if (var.est) { e2 <- integrate(crmht2lgt,-10,10,x1p,y1p,w1p,500,intcpt,abs.tol=0)[[1]] / integrate(crmhlgt,-10,10,x1p,y1p,w1p,500,abs.tol=0)[[1]]; }
    }
    else if (method=="bayes") {
      den <- integrate(crmhlgt,-Inf,Inf,x1p,y1p,w1p,scale,intcpt,abs.tol=0)[[1]]
      est <- integrate(crmhtlgt,-10,10,x1p,y1p,w1p,scale,intcpt,abs.tol=0)[[1]] / den
      if (var.est) { e2 <- integrate(crmht2lgt,-10,10,x1p,y1p,w1p,scale,intcpt,abs.tol=0)[[1]] / den; }
    }
    else { stop(" unknown estimation method"); }
#    est <- min(UB, max(LB, est))
    ptox <- (1 + exp(-intcpt-exp(est)*dosescaled))^{-1}
    if (var.est) {
      post.var <- e2-est^2
      crit <- qnorm(0.5+conf.level/2)
      lb <- est - crit*sqrt(post.var)
      ub <- est + crit*sqrt(post.var)
      ptoxL <- (1 + exp(-intcpt-exp(ub)*dosescaled))^{-1}
      ptoxU <- (1 + exp(-intcpt-exp(lb)*dosescaled))^{-1}
    }
  }
  else { stop(" model specified not available."); }

  if (all(ptox<=target)) { rec <- length(prior); }
  else if (all(ptox>=target)) { rec <- 1; }
  else { rec <- order(abs(ptox-target))[1]; }
  if (!var.est) { post.var <- ptoxL <- ptoxU <- NA; }
  foo <- list(prior=prior, target=target, tox=tox, level=level,
              dosename=dosename, subset=pid[include], estimate=est,
              weights=weights, followup=followup, entry=entry, exit=exit,
              obswin=obswin, scheme=scheme,
              model=model, prior.var=scale^2, post.var=post.var,method=method,
              mtd=rec, include=include, pid=pid, model.detail=model.detail,intcpt=intcpt,
              ptox=ptox, ptoxL=ptoxL, ptoxU=ptoxU, conf.level=conf.level,
              patient.detail=patient.detail,tite=tite,dosescaled=dosescaled)
  class(foo) <- "mtd"
  foo
}

print.mtd <- function(x, dgt=3, model.detail=x$model.detail, patient.detail=x$patient.detail, ...) {
  cat("Today: ", date(), "\n")
  n <- length(x$pid)
  used <- rep(0,n)
  used[x$include] <- 1
  ptox <- round(x$ptox, digits=dgt)
  ptoxL <- round(x$ptoxL, digits=dgt)
  ptoxU <- round(x$ptoxU, digits=dgt)
  if (x$tite) {
    wts <- round(x$weights, digits=dgt)
    if (is.null(x$followup)) { followup <- rep("N/A",n); }
    else { followup <- round(x$followup, digits=dgt); }
    if (patient.detail) {
      cat("DATA SUMMARY (TITE-CRM) \n")
      cat("PID", "\t", "Level", "\t","Toxicity","\t","f/u","\t","Weight","\t","Included","\n")
      for (i in 1:n)
        cat(x$pid[i],"\t",x$level[i],"\t",x$tox[i],"\t\t",followup[i],"\t",wts[i],"\t\t",used[i],"\n")
    }
  }
  else {
    if (patient.detail) {
      cat("DATA SUMMARY (CRM) \n")
      cat("PID", "\t", "Level", "\t", "Toxicity", "\t", "Included", "\n")
      for (i in 1:n) {
        cat(x$pid[i],"\t",x$level[i],"\t",x$tox[i],"\t\t",used[i],"\n")
      }
    }
  }
  cat("\nToxicity probability update (with",x$conf.level*100,"percent probability interval):","\n")
  if (is.null(x$dosename)) {
    cat("Level", "\t", "Prior","\t","n","\t","total.wts","\t","total.tox","\t","Ptox", "\t",
        "LoLmt","\t", "UpLmt", "\n")
    K <- length(x$prior)
    for (k in 1:K) {
      expt <- which(x$level==k & used==1)
      if (!x$tite) { totwts <- length(expt); }
      else { totwts <- round(sum(x$weights[expt]),digits=dgt); }
      cat(k,"\t", x$prior[k], "\t", length(expt), "\t", totwts, "\t\t", sum(x$tox[expt]),"\t\t",ptox[k],"\t",ptoxL[k],"\t",ptoxU[k],"\n")
    }
    cat("Next recommended dose level:", x$mtd,"\n")
  }
  else {
    cat("Dose","\t\t","Level", "\t", "Prior","\t","n","\t","total.wts","\t","total.tox","\t","Ptox", "\t",
        "LoLmt","\t", "UpLmt", "\n")
    K <- length(x$prior)
    for (k in 1:K) {
      expt <- which(x$level==k & used==1)
      if (!x$tite) { totwts <- length(expt); }
      else { totwts <- round(sum(x$weights[expt]),digits=dgt); }
      cat(x$dosename[k],"\t", k,"\t", x$prior[k], "\t", length(expt), "\t", totwts, "\t\t", sum(x$tox[expt]),"\t\t",ptox[k],"\t",ptoxL[k],"\t",ptoxU[k],"\n")
    }
    cat("Next recommended dose level:", x$mtd,"\n")
  }
  cat("Recommendation is based on a target toxicity probability of",x$target,"\n")
  
  if (model.detail) {
    cat("\nEstimation details:\n")
    if (x$model=="empiric") { cat("Empiric dose-toxicity model: p = dose^{exp(beta)}\n"); }
    else { cat("Logistic dose-toxicity model: p = {1 + exp(-a-exp(beta)*dose)}^{-1} with a =",x$intcpt,"\n"); }
    cat("dose =", round(x$dosescaled,digits=dgt),"\n")
    if (x$method=="bayes") {
      cat("Normal prior on beta with mean 0 and variance",x$prior.var,"\n")
      cat("Posterior mean of beta:", round(x$estimate,digits=dgt),"\n")
      cat("Posterior variance of beta:",round(x$post.var,digits=dgt),"\n")
    }
    else if (x$method=="mle") {
      cat("MLE of beta:", round(x$estimate,digits=dgt),"\n")
      cat("Approximate variance of beta:", round(x$post.var,digits=dgt),"\n")
    }
  }
}

plot.mtd <- function(x, ask=FALSE, den=40, ci=TRUE, drlegend=FALSE, ...) {
  ptox <- x$ptox
  ptoxL <- x$ptoxL
  ptoxU <- x$ptoxU
  conf.level <- x$conf.level
  ptox0 <- x$prior
  target <- x$target
  K <- length(ptox)
  include <- x$include
  status <- x$tox[include]
  level <- x$level[include]
  if (x$tite) {
    followup <- x$followup[include]
    exit <- x$exit[include]
    entry <- x$entry[include]
    obswin <- x$obswin
    scheme <- x$scheme
    if (!is.null(exit) & !is.null(entry)) { exit <- x$exit[include]; entry <- x$entry[include]; }
    if (!(is.null(obswin) | is.null(followup)))  {
      hor <- seq(0,obswin,len=den)
      ver <- hor/obswin
      if (scheme=="adaptive") {
        support <- sort(followup[status==1])
        z <- length(support)
        if (z) {
          for (j in 1:den) {
            m <- length(support[support<=hor[j]])
            if (!m)  ver[j] <- hor[j] / support[1] / (z+1)
            else if (m==z)
            ver[j] <- (z + (hor[j]-support[z])/(obswin-support[z])) / (z+1)
            else
              ver[j] <- (m + (hor[j]-support[m])/(support[(m+1)]-support[m])) / (z+1)
          }
        }
      }
    }
  }

  if (x$tite) { choices <- c("all","dose-response curves","trial summary","trial summary (study time)","weight function"); }
  else { choices <- c("all", "dose-response curves","trial summary"); }
  choices <- substring(choices,1,40)
  tmenu <- paste("Plot",choices)
  maxpick <- pick <- 2
  ask.now <- ask
  while (pick <= length(tmenu) + 1) {
    if (ask.now) pick <- menu(tmenu,title="\nMake a plot selection (or 0 to exit):") + 1
    switch(pick,
           return(invisible(x)),   # exit plotting
           ask.now <- FALSE,         # plot next
           {
             if (drlegend) { cat("Please point-click location of legend.\n"); }
             par(las=1,mgp=c(2.8,1,0))
             plot(1:K,ptox, type="l",xlab="Dose level",ylab="Probability of toxicity",ylim=c(0,1))
             points(1:K,ptox, pch="X")
             lines(1:K,ptox0,lty=2)
             points(1:K,ptox0,pch="O")
             abline(h=target,lwd=2)
             if (ci) { lines(1:K,ptoxL,lty=3); lines(1:K,ptoxU,lty=3); }
             if (drlegend) {
               legend(locator(1), c("Updated curve","Prior curve"),lty=c(1,2),pch=c("X","O"),bty="n");
               mtext("Dose-toxicity curves based on",length(status),"subjects",line=0.5);
             }
             else {  mtext("Prior (dashed) and updated (solid) dose-toxicity curves",line=0.5); }
           }
           ,
           {
             par(las=1)
             nontox <- which(status==0)
             plot(nontox,level[nontox],pch="O",ylim=c(1,K),xlab="Patient number",ylab="Dose level",xlim=c(1,length(status)))
             points((1:length(status))[-nontox],level[-nontox],pch="X")
             mtext("Each point represents a patient",line=2)
             mtext("A circle indicates no toxicity, a cross toxicity",line=0.5)
           }
           ,
           {
             if (is.null(entry) | is.null(exit)) {
               cat("\nNo plot output of the study-time trial summary:\n")
               cat("Patient entry and exit information is required.\n");
             }
             else {
               n <- length(entry)
               entrycal <- entry - entry[1]
               exitcal <- exit - entry[1]
               xmax <- max( max(exitcal[status==1]), max(entrycal) )
               plot(entrycal, level, type="n", xlab="Study time (days)", ylab="Dose level",xlim=c(0,xmax),ylim=c(1,K))
               mtext("Each number represents a patient",line=2)
               mtext("An uncircled number indicates entry time, circled time of toxicity",line=0.5)
               text(entrycal,level,as.character(1:n),cex=0.7)
               leveljit <- level + 0.1
               text(exitcal[status==1],leveljit[status==1], as.character((1:n)[status==1]), cex=0.7)
               points(exitcal[status==1],leveljit[status==1], pch="O", cex =1.7 )
             }
           }
           ,
           {
             if (! (is.null(followup) | is.null(obswin)) & scheme=="adaptive") {
               plot(hor,ver,type="n",xlab="Follow-up time of patient (days)",ylab="Weight",ylim=c(0,1))
               lines(hor,ver,lty=5)
               lines(c(0,obswin),c(0,1))
               mtext("An adaptive weight function (dashed line)",line=0.5)
               if (z) points(myjitter(support),rep(0,z),pch=16)
             }
             else {
               cat("\nNo plot output of weight function:\n")
               cat("Times-to-toxicity and censoring times are not available.\n")
             }
           }
           )
    if (!ask.now)  pick <- pick + 1
  }
  return(invisible(x))
}

myjitter <- function(x,factor=1) {
  z <- diff(range(x[!is.na(x)]))
  z <- factor * (z/50)
  if (!z)  z <- abs(mean(x)/50)
  x + runif(length(x),-z,z)
}

cohere <- function(prior,target,x0,method="bayes",model="empiric",intcpt=3,
                   scale=sqrt(1.34), detail=TRUE) {
  n <- length(x0)
  K <- length(prior)
  if (any(x0>K)) { stop(" Initial design escalates beyond the highest test dose"); }
  if (all(x0<K)) { stop(" Initial design never reaches the highest test dose"); }
  
  coh <- rep(TRUE,(n-1)); vlevel <- rep(NA,(n-1))
  for (i in 1:(n-1)) {
    if (method=="bayes" | i>1) {
      level <- x0[1:i]
      y <- c(rep(0,(i-1)),1)
      obj <- crm(prior,target,y,level,method=method,model=model,intcpt=intcpt,scale=scale,var.est=FALSE)
      est <- obj$estimate
      if (model=="logistic") {
	  pest = exp( intcpt + exp(est)*( log(prior/(1-prior)) - intcpt ) )
	  pest = pest/(1+pest)
 	}
	else {
	  pest <- prior^exp(est)
	}
      cur <- order(abs(pest-target))[1]
      if (cur > level[i]) { 
        coh[i] <- FALSE; vlevel[i] <- cur;
        if (detail) {
          cat("\nIncoherent escalation could occur after n =",i,"patients; e.g.\n")
          cat("Level","\t",level,"\n"); cat("tox","\t",y,"\n");
          cat("Recommended level:",cur,"\n\n")
        }
      }
    }
  }
  if (all(coh)) { msg <- "Coherent" }
  else {
    ind <- which(!coh)
    m <- length(ind)
    msg <- "Incoherent: Take a less conservative initial design.  Try the function `getinit'."
  }
  
  foo <- list(message=msg,coherent=coh,x0=x0,
              prior=prior,target=target,method=method,model=model,intcpt=intcpt,prior.var=scale^2,property="cohere")
  class(foo) <- "dxcrm"
  foo
}

getinit <- function(prior, target, n, nK=round(n/3), method="bayes", model="empiric",intcpt=3,
                    scale=sqrt(1.34), detail=FALSE) {
  K <- length(prior)
  if (K<=2) { stop(" Number of dose levels less than 3!"); }
  xtab <- c(rep(0,(K-2)),1,(n-1))
  x0 <- rep(1:K, xtab) 
  x1 <- x0
  counter <- K-2
  obj <- cohere(prior,target,x1,method=method,model=model,intcpt=intcpt,scale=scale,detail=detail)
  while (obj$message=="Coherent") {
    xtab[counter] <- xtab[counter]+1
    xtab[K] <- xtab[K]-1
    if (xtab[K] < nK) { x0 <- x1; break; }
    if (counter==1) { counter <- K-1; }
    else { counter <- counter - 1; }
    if (detail) { cat("An intermediate design:", xtab,"\n"); }
    x0 <- x1
    x1 <- rep(1:K, xtab)
    obj <- cohere(prior,target,x1,method=method,model=model,intcpt=intcpt,scale=scale,detail=detail)
  }
  return(x0)
}
  
print.dxcrm <- function(x, dgt=3, ...) {
  if (x$property=="cohere") {
    x0 <- x$x0; 
    cat("\nMessage:", x$message,"\n")
    cat("\nInitial design:", x0, "\n")
    cat("\nCRM setup:\n")
    if (x$model=="empiric" & x$method=="bayes") {
      obj <- list(prior=x$prior, target=x$target, model=x$model, method=x$method, prior.var=x$prior.var);
    }
    else if (x$model=="empiric" & x$method=="mle") {
      obj <- list(prior=x$prior, target=x$target, model=x$model, method=x$method);
    }
    else if (x$model=="logistic" & x$method=="bayes") {
      obj <- list(prior=x$prior, target=x$target, model=x$model, intercept=x$intcpt, method=x$method, prior.var=x$prior.var);
    }
    else if (x$model=="logistic" & x$method=="mle") {
      obj <- list(prior=x$prior, target=x$target, model=x$model, intercept=x$intcpt,method=x$method);
    }
    print(obj)
  }
  else if (x$property=="sens") {
    K <- length(x$prior)
    if (x$detail) {
      cat("\nThe H sets for the model setup:\n")
      rownames(x$Hset) <- paste("H -",1:K)
      colnames(x$Hset) <- c("LoLmt","UpLmt")
      print(x$Hset)
    }

    cat("\nTrue","\t\t","Indifference interval","\n")
    cat("MTD","\t\t","LoLmt","\t\t","UpLmt","\n")
    for (k in 1:K) cat(k,"\t\t",round(x$iint[k,1],digits=dgt),"\t\t",round(x$iint[k,2],digits=dgt),"\n")
    cat("for the",x$model,"model with prior guesses",x$prior,"\n")
    if (x$model=="logistic") { cat("and intercept", x$intcpt,"\n"); }
    cat("With this model, the CRM will eventually choose a dose with \n")
    cat("  toxicity probability between", signif(min(x$iint[(2:K),1]),digits=dgt),"and",signif(max(x$iint[(1:(K-1)),2]),digits=dgt),
        "while targeting at",x$target,"\n\n")
    cat("Consistency will hold if\n")
    cat("(i) the dose below the MTD does not exceed the Lower limit\n")
    cat("(ii) the dose above the MTD is more toxic than the Upper limit\n\n")
  }
}

titesim <- function(PI, prior, target, n, x0, nsim=1, restrict=TRUE, obswin=1, tgrp=obswin, rate=1, accrual="fixed",
                    surv="uniform", scheme="linear", count=TRUE, method="bayes",model="empiric",intcpt=3,scale=sqrt(1.34),seed=1009) {
  if (nsim==1) {
    foo <- onetite(PI,prior,target,n,x0,restrict=restrict,obswin=obswin,tgrp=tgrp,rate=rate,accrual=accrual,
                   surv=surv,scheme=scheme,method=method,model=model,intcpt=intcpt,scale=scale,seed=seed)
  }
  else {
    foo <- mtite(nsim,PI,prior,target,n,x0,restrict=restrict,obswin=obswin,tgrp=tgrp,rate=rate,accrual=accrual,count=count,
                 surv=surv,scheme=scheme,method=method,model=model,intcpt=intcpt,scale=scale,seed=seed)
  }
  foo
}

onetite <- function(PI, prior, target, n,x0, restrict=TRUE,
                    obswin=1, tgrp=obswin, rate=1, accrual="fixed", surv="uniform", scheme="linear",
                    method="bayes", model="empiric", intcpt=3, scale=sqrt(1.34), seed=1099) {
  set.seed(seed)
  if (accrual=="fixed") { next.arrival <- obswin/rate; }
  else if (accrual=="poisson") { next.arrival <- rexp(1, rate/obswin); }
  
  if (length(x0)>1) {
    if (length(x0)!=n) { stop(" Initial design has a different sample size than that is specified!"); }
    bethat <- rep(0,n)
    u <- y <- level <- arrival <- rep(NA,n)
    m <- 1
    while (TRUE) {
      level[m] <- cur <- x0[m]
      if (is.na(arrival[m])) { arrival[m] <- next.arrival; }
      if (is.na(y[m])) {
        if (surv=="uniform") {
          y[m] <- ynew <- rbinom(1,1,PI[cur])
          if (ynew) { unew <- runif(1,0,obswin); }
          else { unew <- Inf; }
          u[m] <- unew
          utox <- u + arrival
        }
      }
      if (accrual=="fixed") { next.arrival <- next.arrival + obswin/rate; }
      else if (accrual=="poisson") { next.arrival <- next.arrival + rexp(1, rate/obswin); }
      B <- rep(0,n); B[utox<=next.arrival] <- 1;
      if (sum(B)>0 | m==(n-1)) { break; }
      if (x0[m+1]==cur | (next.arrival-arrival[m])>=tgrp) { m <- m+1; }
    }   
    if (m==(n-1)) {
      if (sum(B)==0) { cur <- x0[n]; }
      else {
        censor <- pmin(next.arrival, utox) - arrival;
        followup <- pmin(censor, obswin)
        if (mean(B[1:m])==1) {
          obj <- titecrm(prior, target, B[1:m], level[1:m], followup=followup[1:m],obswin=obswin,scheme=scheme,model=model,
                         intcpt=intcpt,scale=500,var.est=FALSE)
        }
        else {
          obj <- titecrm(prior, target, B[1:m], level[1:m], followup=followup[1:m],obswin=obswin, scheme=scheme,method=method,
                         model=model,intcpt=intcpt,scale=scale,var.est=FALSE)
        }
        if (restrict) { cur <- min(obj$mtd, (cur+1)); }
        else { cur <- obj$mtd; }
        bethat[m+1] <- obj$est
      }
      arrival[n] <- next.arrival
      level[n] <- cur
      if (surv=="uniform") {
        y[n] <- ynew <- rbinom(1,1, PI[cur])
        if (ynew) { unew <- runif(1,0,obswin); }
        else { unew <- Inf; }
        u[n] <- unew; utox <- u + arrival
      }
    }
    else {
      censor <- pmin(next.arrival, utox) - arrival
      followup <- pmin(censor, obswin)
      if (mean(B[1:m])==1) {
        obj <- titecrm(prior, target, B[1:m], level[1:m], followup=followup[1:m],obswin=obswin,scheme=scheme,model=model,
                       intcpt=intcpt,scale=500,var.est=FALSE)
      }
      else {
        obj <- titecrm(prior, target, B[1:m], level[1:m], followup=followup[1:m],obswin=obswin, scheme=scheme,method=method,model=model,
                     intcpt=intcpt,scale=scale,var.est=FALSE)
#        if (accrual=="fixed") { next.arrival <- next.arrival + obswin/rate; }
#        else if (accrual=="poisson") { next.arrival <- next.arrival + rexp(1, rate/obswin); }
      }
      if (restrict) { cur <- min(obj$mtd, (cur+1)); }
      else { cur <- obj$mtd; }
      bethat[m+1] <- obj$est
      
      for (i in (m+1):(n-1)) {
        arrival[i] <- next.arrival; level[i] <- cur;
        if (surv=="uniform") {
          y[i] <- ynew <- rbinom(1,1,PI[cur])
          if (ynew) unew <- runif(1,0,obswin)
          else unew <- Inf
          u[i] <- unew; utox <- u + arrival;
        }
        if (accrual=="fixed") next.arrival <- next.arrival + obswin/rate
        else if (accrual=="poisson") next.arrival <- next.arrival + rexp(1,rate/obswin)
        B <- rep(0,n);  B[utox<=next.arrival] <- 1;
        censor <- pmin(next.arrival,utox) - arrival; followup <- pmin(censor,obswin);
        if (mean(B[1:i])==1) {
          obj <- titecrm(prior, target, B[1:i], level[1:i], followup=followup[1:i],obswin=obswin,scheme=scheme,model=model,
                         intcpt=intcpt,scale=500,var.est=FALSE)
        }
        else {
          obj <- titecrm(prior, target, B[1:i], level[1:i], followup=followup[1:i],obswin=obswin, scheme=scheme,method=method,
                         model=model,
                         intcpt=intcpt,scale=scale,var.est=FALSE)
        }
        if (restrict) { cur <- min(obj$mtd, (cur+1)); }
        else { cur <- obj$mtd; }
        bethat[i+1] <- obj$est
      }

      arrival[n] <- next.arrival; level[n] <- cur;
      if (surv=="uniform") {
        y[n] <- ynew <- rbinom(1,1,PI[cur])
        if (ynew) unew <- runif(1,0,obswin)
        else unew <- Inf
        u[n] <- unew; utox <- u + arrival
      }
    }
    if (method=="mle" & (sum(y)==0|mean(y)==1)) {
      finalobj <- titecrm(prior,target,y,level,weights=rep(1,n),model=model,intcpt=intcpt,scale=500)
      msg <- "Warning: mle is approximated"
    }
    else {
      finalobj <- titecrm(prior, target, y, level, weights=rep(1,n),method=method,model=model,intcpt=intcpt,scale=scale)
      msg <- "Okay"
    }
    cur <- finalobj$mtd
    est <- finalobj$est
  }
  else {
    bethat <- 0
    u <- y <- level <- arrival <- NULL
    cur <- x0
    if (method=="mle") { stop(" Require an initial design for mle-CRM!"); }
    for (i in 1:(n-1)) {
      arrival <- c(arrival, next.arrival)
      level <- c(level, cur)
      if (surv=="uniform") {
        ynew <- rbinom(1,1,PI[cur])
        if (ynew) unew <- runif(1,0,obswin)
        else unew <- Inf
        y <- c(y,ynew); u <- c(u,unew); utox <- u + arrival;
      }
      if (accrual=="fixed") { next.arrival <- next.arrival + obswin/rate; }
      else if (accrual=="poisson") { next.arrival <- next.arrival + rexp(1, rate/obswin); }
      B <- rep(0,length(y));  B[utox<=next.arrival] <- 1;
      censor <- pmin(next.arrival,utox) - arrival
      followup <- pmin(censor,obswin)
      obj <- titecrm(prior,target,B,level,followup=followup,obswin=obswin, scheme=scheme, method=method,model=model,intcpt=intcpt,scale=scale,var.est=FALSE)
      if (restrict) { cur <- min(obj$mtd, (cur+1)); }
      else { cur <- obj$mtd; }
      bethat <- c(bethat, obj$est)
    }
    arrival <- c(arrival, next.arrival)
    level <- c(level,cur)
    if (surv=="uniform") {
      ynew <- rbinom(1,1,PI[cur])
      if (ynew) unew <- runif(1,0,obswin)
      else unew <- Inf
      y <- c(y,ynew); u <- c(u,unew); utox <- u + arrival;
    }
    finalobj <- titecrm(prior,target,y,level,weights=rep(1,n), method=method,model=model,intcpt=intcpt,scale=scale)
    cur <- finalobj$mtd
    est <- finalobj$est
    msg <- "Okay"
  }
  if (length(x0)==1) { design <- paste("TITE-CRM starting at dose",x0); }
  else { design <- "Two-stage TITE-CRM"; }
      
  foo <- list(PI=PI, prior=prior, target=target, n=n, x0=x0, nsim=1,
              MTD=cur, level=level, tox=y, beta.hat=bethat, final.est=est, 
              arrival=arrival,  toxicity.time=u, toxicity.study.time=utox,
              design=design, method=method, prior.var=scale^2, model=model, intcpt=intcpt,
              restriction=restrict,  seed=seed, tite=TRUE, dosescaled=finalobj$dosescaled, msg=msg,
              obswin=obswin, tgrp=tgrp, rate=rate,accrual=accrual, scheme=scheme,
              post.var=finalobj$post.var, ptox=finalobj$ptox, ptoxL=finalobj$ptoxL, ptoxU=finalobj$ptoxU,
              conf.level=finalobj$conf.level)

  class(foo) <- "sim"
  foo
}

mtite <- function(nsim, PI, prior, target, n, x0, restrict=TRUE, obswin=1, tgrp=obswin, rate=1, accrual="fixed",
                  surv="uniform", scheme="linear", count=TRUE, method="bayes",model="empiric",intcpt=3,scale=sqrt(1.34),seed=1009) {
  set.seed(seed)
  if (nsim<=1) { stop(" nsim less than 2!  For single simulation replicate, use function `onetite'."); }
  nexpt <- ntox <- sel <- rep(0,length(prior))
  BETAHAT <- matrix(rep(NA,nsim*n), nrow=nsim)
  final.est <- DURATION <- rep(NA,nsim)

  for (r in 1:nsim) {
    if (count) { cat("simulation number:",r,"\n"); }
    
    if (accrual=="fixed") { next.arrival <- obswin/rate; }
    else if (accrual=="poisson") { next.arrival <- rexp(1, rate/obswin); }
    
    if (length(x0)>1) {
      if (length(x0)!=n) { stop(" Initial design has a different sample size than that is specified!"); }
      bethat <- rep(0,n)
      u <- y <- level <- arrival <- rep(NA,n)
      m <- 1
      while (TRUE) {
        level[m] <- cur <- x0[m]
        if (is.na(arrival[m])) { arrival[m] <- next.arrival; }
        if (is.na(y[m])) {
          if (surv=="uniform") {
            y[m] <- ynew <- rbinom(1,1,PI[cur])
            if (ynew) { unew <- runif(1,0,obswin); }
            else { unew <- Inf; }
            u[m] <- unew
            utox <- u + arrival
          }
        }
        if (accrual=="fixed") { next.arrival <- next.arrival + obswin/rate; }
        else if (accrual=="poisson") { next.arrival <- next.arrival + rexp(1, rate/obswin); }
        B <- rep(0,n); B[utox<=next.arrival] <- 1;
        if (sum(B)>0 | m==(n-1)) { break; }
        if (x0[m+1]==cur | (next.arrival-arrival[m])>=tgrp) { m <- m+1; }
      }   
      if (m==(n-1)) {
        if (sum(B)==0) { cur <- x0[n]; }
        else {
          censor <- pmin(next.arrival, utox) - arrival;
          followup <- pmin(censor, obswin)
          if (mean(B[1:m])==1) {
            obj <- titecrm(prior, target, B[1:m], level[1:m], followup=followup[1:m],obswin=obswin,scheme=scheme,model=model,
                           intcpt=intcpt,scale=500,var.est=FALSE)
          }
          else {
            obj <- titecrm(prior, target, B[1:m], level[1:m], followup=followup[1:m],obswin=obswin, scheme=scheme,method=method,model=model,
                           intcpt=intcpt,scale=scale,var.est=FALSE)
          }
          if (restrict) { cur <- min(obj$mtd, (cur+1)); }
          else { cur <- obj$mtd; }
          bethat[m+1] <- obj$est
        }
        arrival[n] <- next.arrival
        level[n] <- cur
        if (surv=="uniform") {
          y[n] <- ynew <- rbinom(1,1, PI[cur])
          if (ynew) { unew <- runif(1,0,obswin); }
          else { unew <- Inf; }
          u[n] <- unew; utox <- u + arrival
        }
      }
      else {
        censor <- pmin(next.arrival, utox) - arrival
        followup <- pmin(censor, obswin)
        if (mean(B[1:m])==1) {
          obj <- titecrm(prior, target, B[1:m], level[1:m], followup=followup[1:m],obswin=obswin,scheme=scheme,model=model,
                         intcpt=intcpt,scale=500,var.est=FALSE)
        }
        else {
          obj <- titecrm(prior, target, B[1:m], level[1:m], followup=followup[1:m],obswin=obswin, scheme=scheme,method=method,model=model,
                         intcpt=intcpt,scale=scale,var.est=FALSE)
        }
        if (restrict) { cur <- min(obj$mtd, (cur+1)); }
        else { cur <- obj$mtd; }
        bethat[m+1] <- obj$est
#        if (accrual=="fixed") { next.arrival <- next.arrival + obswin/rate; }
#        else if (accrual=="poisson") { next.arrival <- next.arrival + rexp(1, rate/obswin); }
        
        for (i in (m+1):(n-1)) {
          arrival[i] <- next.arrival; level[i] <- cur;
          if (surv=="uniform") {
            y[i] <- ynew <- rbinom(1,1,PI[cur])
            if (ynew) unew <- runif(1,0,obswin)
            else unew <- Inf
            u[i] <- unew; utox <- u + arrival;
          }
          if (accrual=="fixed") next.arrival <- next.arrival + obswin/rate
          else if (accrual=="poisson") next.arrival <- next.arrival + rexp(1,rate/obswin)
          B <- rep(0,n);  B[utox<=next.arrival] <- 1;
          censor <- pmin(next.arrival,utox) - arrival; followup <- pmin(censor,obswin);
          if (mean(B[1:i])==1) {
            obj <- titecrm(prior, target, B[1:i], level[1:i], followup=followup[1:i],obswin=obswin,scheme=scheme,model=model,
                           intcpt=intcpt,scale=500,var.est=FALSE)
          }
          else {
            obj <- titecrm(prior, target, B[1:i], level[1:i], followup=followup[1:i],obswin=obswin, scheme=scheme,method=method,model=model,
                           intcpt=intcpt,scale=scale,var.est=FALSE)
          }
          if (restrict) { cur <- min(obj$mtd, (cur+1)); }
          else { cur <- obj$mtd; }
          bethat[i+1] <- obj$est
        }
        BETAHAT[r,] <- bethat
        
        arrival[n] <- next.arrival; level[n] <- cur;
        if (surv=="uniform") {
          y[n] <- ynew <- rbinom(1,1,PI[cur])
          if (ynew) unew <- runif(1,0,obswin)
          else unew <- Inf
          u[n] <- unew; utox <- u + arrival
        }
      }
      if (method=="mle" & (sum(y)==0|mean(y)==1)) {
        finalobj <- titecrm(prior,target,y,level,weights=rep(1,n),model=model,intcpt=intcpt,scale=500,var.est=FALSE)
        msg <- "Warning: mle is approximated"
      }
      else {
        finalobj <- titecrm(prior, target, y, level, weights=rep(1,n),method=method,model=model,intcpt=intcpt,scale=scale,var.est=FALSE)
        msg <- "Okay"
      }
      cur <- finalobj$mtd
      est <- finalobj$est
    }
    else {
      if (method=="mle") { stop(" Require an initial design for mle-CRM!"); }
      bethat <- 0
      u <- y <- level <- arrival <- NULL
      cur <- x0
      
      for (i in 1:(n-1)) {
        arrival <- c(arrival, next.arrival)
        level <- c(level, cur)
        if (surv=="uniform") {
          ynew <- rbinom(1,1,PI[cur])
          if (ynew) unew <- runif(1,0,obswin)
          else unew <- Inf
          y <- c(y,ynew); u <- c(u,unew); utox <- u + arrival;
        }
        if (accrual=="fixed") { next.arrival <- next.arrival + obswin/rate; }
        else if (accrual=="poisson") { next.arrival <- next.arrival + rexp(1, rate/obswin); }
        B <- rep(0,length(y));  B[utox<=next.arrival] <- 1;
        censor <- pmin(next.arrival,utox) - arrival
        followup <- pmin(censor,obswin)
        obj <- titecrm(prior,target,B,level,followup=followup,obswin=obswin, scheme=scheme, method=method,model=model,intcpt=intcpt,scale=scale,var.est=FALSE)
        if (restrict) {
          cur <- min(obj$mtd, (cur+1))
        }
        else {
          cur <- obj$mtd
        }
        bethat <- c(bethat, obj$est)
      }
      BETAHAT[r,] <- bethat
      arrival <- c(arrival, next.arrival)
      level <- c(level,cur)
      if (surv=="uniform") {
        ynew <- rbinom(1,1,PI[cur])
        if (ynew) unew <- runif(1,0,obswin)
        else unew <- Inf
        y <- c(y,ynew); u <- c(u,unew); utox <- u + arrival;
      }
      finalobj <- titecrm(prior,target,y,level,weights=rep(1,n), method=method,model=model,intcpt=intcpt,scale=scale,var.est=FALSE)
      cur <- finalobj$mtd
      est <- finalobj$est
      msg <- "Okay"
    }

    sel[cur] <- sel[cur] + 1
    final.est[r] <- est
    DURATION[r] <- max(arrival) + obswin
    for (k in 1:length(prior)) {
      nexpt[k] <- nexpt[k] + length(which(level==k))
      ntox[k] <- ntox[k] + length(which(y==1 & level==k))
    }
  }
  sel <- sel/nsim
  nexpt <- nexpt/nsim
  ntox <- ntox/nsim
  if (length(x0)==1) { design <- paste("TITE-CRM starting at dose",x0); }
  else { design <- "Two-stage TITE-CRM"; }
  
  foo <- list(PI=PI, prior=prior, target=target, n=n, x0=x0, nsim=nsim,
              MTD=sel, level=nexpt, tox=ntox, beta.hat=BETAHAT, final.est=final.est,
              Duration=DURATION,
              design=design, method=method, prior.var=scale^2, model=model, intcpt=intcpt,
              restriction=restrict, seed=seed, tite=TRUE, dosescaled=finalobj$dosescaled, msg=msg,
              obswin=obswin, tgrp=tgrp, rate=rate, accrual=accrual, scheme=scheme)
  
  class(foo) <- "sim"
  foo
}


print.sim <- function(x, dgt=3, patient.detail=TRUE, ...) {
  n <- x$n
  PI <- x$PI
  prior <- round(x$prior, digits=dgt)
  target <- x$target
  K <- length(prior)
  
  if (x$nsim==1) {
    y <- x$tox
    level <- x$level
    
    bethat <- signif(x$beta.hat, digits=1)
    est <- round(x$final.est, digits=dgt)
    ptox <- round(x$ptox, digits=dgt)
    ptoxL <- round(x$ptoxL, digits=dgt)
    ptoxU <- round(x$ptoxU, digits=dgt)

    if (patient.detail) {
      if (x$tite) {
        arrival <- signif(x$arrival, digits=dgt)
        utox <- round(x$toxicity.study.time, digits=dgt)
        u <- round(x$toxicity.time, digits=dgt)
        tevent <- round( c(arrival,utox), digits=dgt)
        pid <- rep(1:n, 2)
        event <- c(rep("enrol",n), rep("TOX",n))
        level2 <- c(level, level)
        est2 <- round( c(bethat, rep(NA,n)), digits=dgt)
        o <- order(tevent)
        tevent <- tevent[o]; pid <- pid[o]; event <- event[o]; level2 <- level2[o]; est2 <- est2[o];
        ind <- which(tevent<Inf)
        tevent <- tevent[ind]; pid <- pid[ind]; event <- event[ind]; level2 <- level2[ind]; est2 <- est2[ind];
        m <- length(ind)
        
        cat("Trial summary on study time\n")
        cat("Time \t PID \t Event \t Level \t Beta \n")
        for (j in 1:m) { cat(tevent[j],"\t", pid[j], "\t", event[j], "\t", level2[j], "\t", est2[j],"\n"); }
        
        cat("\nPatient summary (TITE-CRM) \n")
        cat("PID \t Arrive \t Beta \t Level \t Tox \t Tox.time \n")
        for (i in 1:n) {
          cat(i,"\t", arrival[i], "\t\t", bethat[i],"\t", level[i], "\t", y[i], "\t", u[i],"\n")
        }
      }
      else {
        cat("Patient summary (CRM) \n")
        cat("PID", "\t", "Beta", "\t", "Level", "\t", "Toxicity", "\n")
        for (i in 1:n) {
          cat(i,"\t",bethat[i],"\t", level[i],"\t",y[i],"\n")
        }
      }
    }

    cat("\nToxicity probability summary (with",x$conf.level*100,"percent probability interval):","\n")
    cat("Level", "\t", "Ptrue", "Prior","\t","n","\t","ntox","\t","Posterior", "\t",
        "LoLmt","\t", "UpLmt", "\n")

    ntox <- nexpt <- rep(0,K)
    for (k in 1:K) {
      nexpt[k] <- length(which(level==k))
      ntox[k] <- length(which(level==k & y==1))
      cat(k,"\t", PI[k], "\t", prior[k], "\t", nexpt[k], "\t", ntox[k],"\t",ptox[k],"\t\t",ptoxL[k],"\t",ptoxU[k],"\n")
    }
    PIjit <- PI; PIjit[1] <- PI[1]-(1e-5); PIjit[K] <- PI[K]+(1e-5);
    cat("True MTD:", order(abs(PIjit-target))[1],"\tEstimated MTD:", x$MTD,"\tTarget DLT rate:",target,"\n")

    cat("\nThis trial is generated by a",x$design,"\n")
    if (length(x$x0)>1) {
      cat("Dose escalation proceeds as follows before any toxicity is seen:")
      xtab <- cbind(1:K,rep(NA,K))
 	for (k in 1:K) { xtab[k,2] = length(which(x$x0==k)); }
      colnames(xtab) <- c("dose.level","cohort.size")
      rownames(xtab) <- rep("",K)
      print(t(xtab));
    }

    if (x$restrict) {
      cat("\nRestrictions apply to avoid\n")
      cat("\t (1) Skipping doses in escalation;\n")
      cat("\t (2) Escalation immediately after a toxic outcome.\n")
    }

    cat("\nThe working model is",x$model,"\n")
    if (x$model=="empiric") { cat("\tptox = dose^{exp(beta)} with doses =",round(x$dosescaled,digits=dgt),"\n"); }
    else { cat("\tlogit(ptox) = a + exp(beta)*dose, with a =",x$intcpt,"\n\tand doses =",signif(x$dosescaled,digits=dgt),"\n"); }

    if (x$method=="bayes") {
      cat("\tand beta is estimated by its posterior mean \n\tassuming a normal prior with mean 0 and variance",x$prior.var,"\n")
    }
    else if (x$method=="mle") {
      cat("\tand beta is estimated by its mle\n")
    }

    cat("\nThe final estimate of beta",round(x$final.est, digits=dgt))        
    if (x$method=="bayes") {
      cat(" with posterior variance",round(x$post.var,digits=dgt),"\n")
    }
    else if (x$method=="mle") {
      cat(" with variance", round(x$post.var,digits=dgt),"\n")
      if (x$msg!="Okay") { print(x$msg); }
    }

    if (x$tite) {
      cat("\nThe",x$scheme,"function is used to assign weights to patients.\n\n")
      cat("Patient arrival is modeled as a", x$accrual,"process\n")
      cat("\twith rate",x$rate,"patients per",x$obswin,"time units (= observation window).\n")
      if (length(x$x0)>1) {
        cat("\tA minimum waiting time of",x$tgrp,"time units is imposed\n")
        cat("\tbetween two dose cohorts in the initial stage.\n")
      }
    }
  }
  else if (x$nsim > 1) {
    oc <- t(cbind(PI, prior, x$MTD, x$level, x$tox))
    colnames(oc) <- as.character(1:K)
    rownames(oc) <- c("Truth","Prior","Selected","Nexpt","Ntox")

    cat("\nNumber of simulations:\t",x$nsim,"\n")
    cat("Patient accrued:\t",n,"\n")
    cat("Target DLT rate:\t", target, "\n")
    print(round(oc,digits=dgt))

    if (x$tite) {
      cat("\nThe distribution of trial duration:\n")
      print(summary(x$Dur))
    }
    
    cat("\nThe trials are generated by a",x$design,"\n")
    if (length(x$x0)>1) {
      cat("Dose escalation proceeds as follows before any toxicity is seen:")
	xtab <- cbind(1:K,rep(NA,K))
 	for (k in 1:K) { xtab[k,2] = length(which(x$x0==k)); }
      colnames(xtab) <- c("dose.level","cohort.size")
      rownames(xtab) <- rep("",K)
      print(t(xtab))
    }

    if (x$restrict) {
      cat("\nRestriction apply to avoid\n")
      cat("\t (1) Skipping doses in escalation;\n")
      cat("\t (2) Escalation immediately after a toxic outcome.\n")
    }

    cat("\nThe working model is",x$model,"\n")
    if (x$model=="empiric") { cat("\tptox = dose^{exp(beta)} with doses =",round(x$dosescaled,digits=dgt),"\n"); }
    else { cat("\tlogit(ptox) = a + exp(beta)*dose, with a =",x$intcpt,"\n\tand doses =",signif(x$dosescaled,digits=dgt),"\n"); }
    x0 <- x$x0

        
    if (x$method=="bayes") {
      cat("\tand beta is estimated by its posterior mean \n\tassuming a normal prior with mean 0 and variance",x$prior.var,"\n")
    }
    else if (x$method=="mle") {
      cat("\tand beta is estimated by its mle\n")
    }

    if (x$tite) {
      cat("\nThe",x$scheme,"function is used to assign weights to patients.\n\n")
      cat("Patient arrival is modeled as a", x$accrual,"process\n")
      cat("\twith rate",x$rate,"patients per",x$obswin,"time units (= observation window).\n\n")
      if (length(x$x0)>1) {
        cat("\tA minimum waiting time of",x$tgrp,"time units is imposed\n")
        cat("\tbetween two dose cohorts in the initial stage.\n")
      }
    }

  }
}

plot.sim <- function(x, ask=FALSE, den=40, ci=TRUE, drlegend=FALSE, n0=1, ...) {
  PI <- x$PI
  prior <- x$prior
  target <- x$target
  n <- x$n
  x0 <- x$x0
  nsim <- x$nsim
  K <- length(prior)
  
  if (nsim==1) {
    ptox <- x$ptox
    ptoxL <- x$ptoxL
    ptoxU <- x$ptoxU
    conf.level <- x$conf.level
    ptox0 <- prior
    status <- x$tox
    level <- x$level
    if (x$tite) {
      followup <- x$toxicity.time
      obswin <- x$obswin
      scheme <- x$scheme
      if (!(is.null(obswin) | is.null(followup)))  {
        hor <- seq(0,obswin,len=den)
        ver <- hor/obswin
        if (scheme=="adaptive") {
          support <- sort(followup[status==1])
          z <- length(support)
          if (z) {
            for (j in 1:den) {
              m <- length(support[support<=hor[j]])
              if (!m)  ver[j] <- hor[j] / support[1] / (z+1)
              else if (m==z)
                ver[j] <- (z + (hor[j]-support[z])/(obswin-support[z])) / (z+1)
              else
                ver[j] <- (m + (hor[j]-support[m])/(support[(m+1)]-support[m])) / (z+1)
            }
          }
        }
      }
    }

    if (x$tite) {
      choices <- c("all","dose-response curves","trial summary","trial summary (study time)","weight function");
    }
    else { choices <- c("all", "dose-response curves","trial summary"); }
    choices <- substring(choices,1,40)
    tmenu <- paste("Plot",choices)
    maxpick <- pick <- 2
    ask.now <- ask
    while (pick <= length(tmenu) + 1) {
      if (ask.now) pick <- menu(tmenu,title="\nMake a plot selection (or 0 to exit):") + 1
      switch(pick,
             return(invisible(x)),   # exit plotting
             ask.now <- FALSE,         # plot next
             {
               if (drlegend) { cat("Please point-click location of legend.\n"); }
               par(las=1,mgp=c(2.8,1,0))
               plot(1:K,ptox, type="l",xlab="Dose level",ylab="Probability of toxicity",ylim=c(0,1))
               points(1:K,ptox, pch="X")
               lines(1:K,ptox0,lty=2)
               points(1:K,ptox0,pch="O")
               abline(h=target,lwd=2)
               if (ci) { lines(1:K,ptoxL,lty=3); lines(1:K,ptoxU,lty=3); }
               if (drlegend) {
                 legend(locator(1), c("Updated curve","Prior curve"),lty=c(1,2),pch=c("X","O"),bty="n");
                 mtext("Dose-toxicity curves based on",length(status),"subjects",line=0.5);
               }
               else {  mtext("Prior (dashed) and updated (solid) dose-toxicity curves",line=0.5); }
             }
             ,
             {
               par(las=1)
               nontox <- which(status==0)
               plot(nontox,level[nontox],pch="O",ylim=c(1,K),xlab="Patient number",ylab="Dose level",xlim=c(1,n))
               points((1:length(status))[-nontox],level[-nontox],pch="X")
               mtext("Each point represents a patient",line=2)
               mtext("A circle indicates no toxicity, a cross toxicity",line=0.5)
             }
             ,
             {
               entrycal <- x$arrival
               exitcal <- x$toxicity.study.time
               xmax <- max( max(exitcal[status==1]), max(entrycal) )
               plot(entrycal, level, type="n", xlab="Study time", ylab="Dose level",xlim=c(0,xmax),ylim=c(1,K))
               mtext("Each number represents a patient",line=2)
               mtext("An uncircled number indicates entry time, circled time of toxicity",line=0.5)
               text(entrycal,level,as.character(1:n),cex=0.7)
               leveljit <- level + 0.1
               text(exitcal[status==1],leveljit[status==1], as.character((1:n)[status==1]), cex=0.7)
               points(exitcal[status==1],leveljit[status==1], pch="O", cex =1.7 )
             }
             ,
             {
               if (! (is.null(followup) | is.null(obswin)) & scheme=="adaptive") {
                 plot(hor,ver,type="n",xlab="Follow-up time of patient",ylab="Weight",ylim=c(0,1))
                 lines(hor,ver,lty=5)
                 lines(c(0,obswin),c(0,1))
                 mtext("An adaptive weight function (dashed line)",line=0.5)
                 if (z) points(myjitter(support),rep(0,z),pch=16)
               }
               else {
                 cat("\nNo plot output of weight function:\n")
                 cat("Times-to-toxicity and censoring times are not available.\n")
               }
             }
             )
      if (!ask.now)  pick <- pick + 1
    }
  }
  else {
    AllBeta <- cbind(x$beta.hat, x$final.est)
    if (x$tite) {
      dur <- x$Duration;
      choices <- c("all","beta","beta trajectory","trial duration")
    }
    else {
      choices <- c("all","beta","beta trajectory")
    }
    choices <- substring(choices, 1, 40)
    tmenu <- paste("Plot", choices)
    maxpic <- pick <- 2
    ask.now <- ask
    while (pick <= length(tmenu) + 1) {
      if (ask.now) { pick <- menu(tmenu,title="\nMake a plot selection (or 0 to exit):") + 1; }
      switch(pick,
             return(invisible(x)),   # exit plotting
             ask.now <- FALSE,         # plot next
             {
               if (x$method=="bayes") { maintext <- "Posterior mean of beta";}
               else if (x$method=="mle") { maintext <- "MLE of beta"; }
               maintext <- paste(maintext,"with n =",n)
               par(las=1)
               hist(x$final.est,prob=TRUE,main=maintext,xlab="beta")
             }
             ,
             {
               hor <- n0:(n+1)
               med <- lcb <- ucb <- 1:(n+1)
               objsens <- crmsens(prior,target,model=x$model,intcpt=x$intcpt)
               LB <- objsens$Hset[1,1]; UB <- objsens$Hset[K,2];
               LB <- max(LB, min(AllBeta[,hor]))
               UB <- min(UB, max(AllBeta[,hor]))
               plot(hor, AllBeta[1,hor],type="l",ylim=c(LB,UB),xlab="Patient number",ylab="beta")
               for (i in hor) {
                 lcb[i] <- quantile(AllBeta[,i],0.05)
                 ucb[i] <- quantile(AllBeta[,i],0.95)
                 med[i] <- median(AllBeta[,i])
               }
               for (r in 2:nsim) { lines(hor,AllBeta[r,hor]); }
               lines(hor,lcb[hor],lwd=3); lines(hor,ucb[hor],lwd=3);
               lines(hor,med[hor],lwd=3,lty=2)
             }
             ,
             {
               maintext <- paste("Trial duration for",x$design,"with n =",n)
               hist(dur,prob=TRUE,main=maintext,xlab="Duration")
             }
           )
      if (!ask.now) pick <- pick + 1
    }
  }  
  return(invisible(x))
}

crmsens <- function(prior,target,model="empiric",intcpt=3,eps=1e-6,maxit=100,detail=FALSE) {
  K <- length(prior)
  b <- rep(NA,(K+1))
  Hset <- val <- matrix(rep(NA,K*2),nrow=K)

  if (model=="empiric") {
    dosescaled <- prior
    lb <- log(log((target+1)/2)/log(dosescaled[1]))
    b[1] <- min(-5,lb)
    ub <- log(log(target/2)/log(dosescaled[K]))
    b[K+1] <- max(5,ub)
    LB <- b[1]; UB <- b[K+1];
    for (k in 2:K) {
      lb <- LB; ub <- UB;
      while (TRUE) {
        mid <- (lb+ub)/2
        lhs <- dosescaled[k-1]^exp(mid) + dosescaled[k]^exp(mid)
        if (lhs>(2*target)) { lb <- mid; }
        if (lhs<(2*target)) { ub <- mid; }
        if (abs(lhs-2*target)<eps) { break; }
      }
      LB <- b[k] <- mid
    }
  
    val[1,1] <- NA; val[1,2] <- dosescaled[2]^exp(b[2])
    val[K,1] <- dosescaled[(K-1)]^exp(b[K]); val[K,2] <- NA;
    for (k in 2:(K-1)) {
      val[k,1] <- dosescaled[k-1]^exp(b[k])
      val[k,2] <- dosescaled[k+1]^exp(b[k+1])
    }
  }
  else if (model=="logistic") {
    dosescaled <- log(prior/(1-prior)) - intcpt
    lb <- log(  (log((1+target)/(1-target)) - intcpt)/dosescaled[1] )
    b[1] <- min(-5,lb)
    ub <- log(  (log((target/2)/(1-target/2)) - intcpt)/dosescaled[K] )
    b[K+1] <- max(5, ub)

    LB <- b[1]; UB <- b[K+1];
    for (k in 2:K) {
      lb <- LB; ub <- UB;
      while(TRUE) {
        mid <- (lb+ub)/2
        lhs <- (1 + exp(-intcpt-exp(mid)*dosescaled[k-1]))^{-1} + (1+exp(-intcpt-exp(mid)*dosescaled[k]))^{-1}
        if (lhs>(2*target)) { lb <- mid; }
        if (lhs<(2*target)) { ub <- mid; }
        if (abs(lhs-2*target)<eps) { break; }
      }
      LB <- b[k] <- mid
    }

    val[1,1] <- NA; val[1,2] <- (1 + exp(-intcpt-exp(b[2])*dosescaled[2]))^{-1};
    val[K,1] <- (1 + exp(-intcpt-exp(b[K])*dosescaled[K-1]))^{-1}
    for (k in 2:(K-1)) {
      val[k,1] <- (1 + exp(-intcpt-exp(b[k])*dosescaled[k-1]))^{-1}
      val[k,2] <- (1 + exp(-intcpt-exp(b[k+1])*dosescaled[k+1]))^{-1}
    }
  }
  else { stop( " model specified not available."); }

  for (k in 1:K) { Hset[k,1] <- b[k]; Hset[k,2] <- b[k+1]; }
  if (model=="empiric") a <- list(Hset=Hset, iint=val, target=target, prior=prior, model=model,detail=detail,property="sens")
  if (model=="logistic") a <- list(Hset=Hset,iint=val, target=target,prior=prior, model=model,intcpt=intcpt,
        details=detail,property="sens")
  class(a) <- "dxcrm"
  a
}

nopt = function(apcs,target,K,psi,correction=TRUE) {
  pL = target / ( target + psi*(1-target) )
  pU = target * psi / (1 - target + target*psi)
  sigL = sqrt( target*(1-target) + pL*(1-pL) + 2*pL*(1-target) )
  sigU = sqrt( target*(1-target) + pU*(1-pU) + 2*target*(1-pU) )

  if (!correction) {
    delL = (target - pL)/sigL
    delU = (pU - target)/sigU
    del = (delL + delU)/2
    arg = 1 - K*(1-apcs)/2/(K-1)
    val = ( qnorm(arg) / del )^2
  }
  else {
    n = 1
    while (T) {
	delL = (target-pL+0.5/n)/sigL
	delU = (pU-target-0.5/n)/sigU
	del = (delL + delU)/2
	bindex = 1/K + (1-1/K) * (2 * pnorm(sqrt(n)*del) - 1)
	if (bindex >= apcs) { break; }
	n = n+1
    }
    delL = (target-pL+0.5/n)/sigL
    delU = (pU-target-0.5/n)/sigU
    del = (delL + delU)/2
    arg = 1 - K*(1-apcs)/2/(K-1)
    val = ( qnorm(arg) / del )^2
  }
  val
}

getn = function(apcs,target,nlevel,psi,correction=TRUE,detail=FALSE) {
  msg = "IMPORTANT: The calculated n is intended as a starting point, and is not recommended for use without validation by simulation"
  if (target<0.1 | target>0.3) msg = c(msg, "Warning: The sample size algorithm has been validated for target rate ranging from 0.1 to 0.3 only")
  if (nlevel<4 | nlevel >8) msg = c(msg, "Warning: The sample size algorithm has been validated for 4 to 8 dose levels only")
  if (psi<1.25 | psi>2.5) msg = c(msg, "Warning: The sample size algorithm has been validated for odds ratio (psi) ranging from 1.25 to 2.50 only")
  K = nlevel
  lterm = ( log(apcs/(1-apcs)) - 2.26 + 0.00235*K^2 + 0.7*psi + 1.903/psi ) / 0.854
  b = exp(lterm)/(1 + exp(lterm))
  na = nopt(b,target,K,psi,correction=correction)
  nb = nopt(apcs,target,K,psi,correction=correction)
  eff = nb / na
  if (na > 40)   msg = c(msg, "Warning: The calculation likely involves extrapolation - validate with simulation.")

  val = list(n = ceiling(na), astar=apcs, target=target, nlevel=K, psi=psi, bstar=b, efficiency = eff, correction=correction, na=na,nb=nb, messages=msg, detail=detail)
  class(val) = "crmsize"
  val
}

print.crmsize = function(x, detail=x$detail, ...) {
  cat(" Target rate:\t\t\t", x$target, "\n")
  cat(" Number of dose levels:\t\t",x$nlevel,"\n")
  cat(" Effect size (odds ratio):\t",x$psi,"\n")
  cat(" Required accuracy:\t\t",x$astar,"\n")
  cat(" Calculated sample size:\t",x$n,"\n\n")

  if (detail) {
  	if (x$correction) { cat(" Method: With continuity correction\n"); }
	else { cat(" Method: Without continuity correction\n"); }

	cat(" Required sample size (to first decimal place):\t\t",round(x$na,digits=1),"\n")
	cat(" Sample size lower bound (to first decimal place):\t",round(x$nb,digits=1),"\n")
      cat(" Efficiency relative to optimal benchmark:\t\t",x$efficiency,"\n\n")

      cat("Messages\n")
      print(x$messages)
  }

  cat("\n")
}
    

getprior <- function(halfwidth,target,nu,nlevel,model="empiric", intcpt=3) {
  dosescaled <- prior <- rep(NA,nlevel)
  b <- rep(NA,nlevel+1)
  b[1] <- -Inf
  b[(nlevel+1)] <- Inf
  if (model=="empiric") {
    dosescaled[nu] <- target
    for (k in nu:2) {
      b[k] <- log( log(target+halfwidth)/log(dosescaled[k]) )
      if (nu>1) { dosescaled[k-1] <- exp( log(target-halfwidth)/exp(b[k]) ); }
    }
    if (nu < nlevel) {
      for (k in nu:(nlevel-1)) {
        b[k+1] <- log( log(target-halfwidth)/log(dosescaled[k]) )
        dosescaled[k+1] <- exp( log(target+halfwidth)/exp(b[k+1]) )
      }
    }
    val <- dosescaled
  }
  else if (model=="logistic") {
    dosescaled[nu] <- log(target/(1-target)) - intcpt
    for (k in nu:2) {
      b[k] <- log( (log((target+halfwidth)/(1-target-halfwidth))-intcpt)/dosescaled[k] )
      if (nu>1) { dosescaled[k-1] <- (log((target-halfwidth)/(1-target+halfwidth))-intcpt) / exp(b[k]); }
    }
    if (nu<nlevel) {
      for (k in nu:(nlevel-1)) {
        b[k+1] <- log( (log((target-halfwidth)/(1-target+halfwidth))-intcpt) / dosescaled[k] )
        dosescaled[k+1] <- (log((target+halfwidth)/(1-target-halfwidth))-intcpt) / exp(b[k+1])
      }
    }
    val <- {1 + exp(-intcpt-dosescaled)}^{-1}
  }
  val
}
