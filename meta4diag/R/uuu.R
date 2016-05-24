
.logit <- function(x){return(log(x/(1-x)))}
.invlogit <- function(x){return(exp(x)/(1+exp(x)))}

.to <- function(x){return(log((1+x)/(1-x)))}
.from <- function(x){return(2*exp(x)/(1+exp(x))-1)}

########## summary marginal
.summary.marginal = function(m, level = level){
  m1 = INLA::inla.emarginal(function(x) x^1, m)
  m2 = INLA::inla.emarginal(function(x) x^2, m)
  stdev = sqrt(m2 - m1^2)
  q = INLA::inla.qmarginal(level, m)
  z = c(m1,stdev,q)
  names(z) = c("mean","sd",paste(level,"quant",sep=""))
  return(z)
}

.summary.samples = function(x, level=level){
  m1 = mean(x)
  stdev = sd(x)
  q = quantile(x, probs = level)
  z = c(m1,stdev,q)
  names(z) = c("mean","sd",paste(level,"quant",sep=""))
  return(z)
}

################################################################
################################################################
####
####             Prior Plot functions
####
################################################################
################################################################
.priorInvgamma <- function(a,b,xmax){
  x = seq(0,100,len=1000)
  tau = 1/x
  y = dgamma(tau,shape=a,rate=b)*x^(-2)
  plot(x,y,xlab=expression(sigma^2),ylab=expression(pi(sigma^2)),type="l",xlim=c(0,xmax))
}

.priorSigmaPC <- function(u,alpha,xmax){
  x = seq(0,100,len=1000)
  tau = 1/x
  theta = -log(alpha)/u
  y = 0.5*theta*tau^(-1.5)*exp(-theta/sqrt(tau))*x^(-2)
  plot(x,y,xlab=expression(sigma^2),ylab=expression(pi(sigma^2)),type="l",xlim=c(0,xmax))
}

.priorHalfCauchy = function(gamma,xmax){
  sig = c(seq(0,3,by=0.001),seq(3,10,by=0.1),seq(10,1000,by=1))
  sig = unique(sig)
  pi_sig = 2*gamma/(pi*(sig^2+gamma^2))
  x = sig^2
  y = pi_sig*abs(0.5/sig)
  plot(x,y,xlab=expression(sigma^2),ylab=expression(pi(sigma^2)),type="l",xlim=c(0,xmax))
}

.priorSigmaTnorm <- function(m, v, xmax){
  par.sig = sqrt(v)
  x = seq(0,7,len=1000)
  var.sig = sqrt(x)
  pi.var.sig = 1/par.sig*dnorm((var.sig-m)/par.sig)/(1-pnorm(-m/par.sig))
  y = 0.5*pi.var.sig/var.sig
  plot(x,y,xlab=expression(sigma^2),ylab=expression(pi(sigma^2)),type="l",xlim=c(0,xmax))
}

.priorUniformSig = function(xmax){
  sig = c(seq(0,3,by=0.001),seq(3,10,by=0.1),seq(10,1000,by=1))
  sig = unique(sig)
  pi_sig = rep(1/max(sig),length(sig))
  x = sig^2
  y = pi_sig*abs(0.5/sig)
  plot(x,y,xlab=expression(sigma^2),ylab=expression(pi(sigma^2)),type="l",xlim=c(0,xmax))
}

.priorTableSigma = function(xmax,mrv){
  x = mrv$priorfile$x
  if(any(x<=0)){
    var1dialog <-gtkMessageDialog(mrv$main_window,"destroy-with-parent","warning","ok",
                                    "Please check the input prior, variance should be in [0, infinity]!")
    if (var1dialog$run()==GtkResponseType["ok"]){
      stop("")
    }
    var1dialog$destroy()
  }
  y = mrv$priorfile$y
  plot(x,y,xlab=expression(sigma^2),ylab=expression(pi(sigma^2)),type="l",xlim=c(0,xmax))
}

.priorTableRho = function(xmax,mrv){
  x = mrv$priorfile$x
  if(any(x< -1 | x>1)){
    rhodialog <-gtkMessageDialog(mrv$main_window,"destroy-with-parent","warning","ok",
                                  "Please check the input prior, correlation should be in [-1, 1]!")
    if (rhodialog$run()==GtkResponseType["ok"]){
      stop("")
    }
    rhodialog$destroy()
  }
  y = mrv$priorfile$y
  plot(x,y,xlab=expression(rho),ylab=expression(pi(rho)),type="l",xlim=c(-1,xmax))
}

.priorInvWishart = function(nu,R11,R22,R12,xmax){
  R = matrix(c(R11,R12,R12,R22),2,2)
  detR = det(R)
  gamma_p = gamma(0.5*nu)*gamma(0.5*(nu-1))*sqrt(pi)
  const = detR^(0.5*nu)/((2^nu)*gamma_p)
  sig1 = seq(0.01,20,by=0.1)
  sig2 = seq(0.01,10,by=0.1)
  cor = seq(-0.99999,0.99999,by=0.01)
  M1 = matrix(0,length(sig1),length(sig2))
  for(i in 1:length(sig1)){
    M = lapply(sig2, function(y){
      S1 = matrix(c(sig1[i]^2,0,0,y^2),2,2)
      detS = det(S1)
      RS = R*solve(S1)
      pi_S1 = const*detS^(-0.5*(nu+3))*exp(-0.5*(RS[1,1]+RS[2,2]))
      return(pi_S1)
    })
    M1[i,] = unlist(M)
  }
  M2 = matrix(0,length(sig1),length(cor))
  for(i in 1:length(sig1)){
    M = lapply(cor, function(y){
      S1 = matrix(sig1[i]^2*c(1,y,y,1),2,2)
      detS = det(S1)
      RS = R*solve(S1)
      pi_S1 = const*detS^(-0.5*(nu+3))*exp(-0.5*(RS[1,1]+RS[2,2]))
      return(pi_S1)
    })
    M2[i,] = unlist(M)
  }
  par(mfrow=c(1,2))
  contour(sig1,sig2,M1, nlevel=10,xlim=c(0,xmax),ylim=c(0,xmax),xlab=expression(sigma[1]^2),ylab=expression(sigma[2]^2))
  contour(sig1,cor,M2, nlevel=10,xlim=c(0,xmax),ylim=c(-1,1),xlab=expression(sigma^2),ylab=expression(rho))
}

.priorRhoBetaPlot <- function(a,b){
  # get the density of old normal prior of correlation parameter "rho"
  rho=c(seq(-1,-0.9,len=500),seq(-0.9001,0.8999,len=500),seq(0.9,1,len=500))
  z = 0.5*(rho+1)
  dens = 0.5*dbeta(z,shape1=a,shape2=b)
  plot(-1,-1,xlim=c(-1,1), ylim=c(0,5),xlab=expression(rho),ylab=expression(pi(rho)))
  lines(rho,dens,lty=1,lwd=1)
}

.priorRhoNormalPlot <- function(mean,variance){
  rho=c(seq(-1,-0.9,len=500),seq(-0.9001,0.8999,len=500),seq(0.9,1,len=500))
  transf = function(rho){.logit(0.5*rho+0.5)}
  z = transf(rho)
  lndens = dnorm(z,mean=mean,sd=sqrt(variance))*abs(2/(1-rho^2))
  plot(-1,-1,xlim=c(-1,1), ylim=c(0,5),xlab=expression(rho),ylab=expression(pi(rho)))
  lines(rho,lndens,lty=1,lwd=1)
    
}

.int.density <- function(lambda,rho,rho.ref){
  kld = (1-rho*rho.ref)/(1-rho.ref^2)-0.5*log((1-rho^2)/(1-rho.ref^2))-1
  y = exp(-lambda*sqrt(2*kld))
  return(y)
}

.priorExpRhoS1 <- function(rho.ref,left.portion,Umin,alpha1){
  density.name = "exp"
  lambda1 = c(seq(0.0000000001,0.1,len=100),seq(0.10000001,0.5,len=400),seq(0.5000001,50,len=200))
  f2 = function(lam1,left){return(left/(1-left)*lam1)}
  f1 = left.portion*(.int.density(lambda1,rho.ref,rho.ref)
                     -.int.density(lambda1,Umin,rho.ref))-(left.portion-alpha1)
  g1 = splinefun(f1,lambda1)
  lam1 = g1(0)
  lam2 = f2(lam1=lam1,left=left.portion)
  drho = .rhoPDF(rho.ref=rho.ref, lambda1=lam1, lambda2=lam2, density.name=density.name)
  plot(drho$x,drho$y,xlab=expression(rho),ylab=expression(pi(rho)),type="l",ylim=c(0,drho$y[10000]+0.65))
}


.priorExpRhoS2 <- function(rho.ref,left.portion,Umax,alpha2){
  density.name = "exp"
  right.portion = 1 - left.portion
  lambda2 = c(seq(0.0000000001,0.1,len=100),seq(0.10000001,0.5,len=400),seq(0.5000001,50,len=200))
  f1 = function(lam2,left){return((1-left)/left*lam2)}
  f2 = right.portion*(-.int.density(lambda2,Umax,rho.ref)
                      +.int.density(lambda2,rho.ref,rho.ref))-right.portion+alpha2
  g2 = splinefun(f2,lambda2)
  lam2 = g2(0)
  lam1 = f1(lam2,left.portion)
  drho = .rhoPDF(rho.ref=rho.ref, lambda1=lam1, lambda2=lam2, density.name=density.name)
  plot(drho$x,drho$y,xlab=expression(rho),ylab=expression(pi(rho)),type="l",ylim=c(0,drho$y[10000]+0.65))
}

.priorExpRhoS3 <- function(rho.ref,Umin,alpha1,Umax,alpha2){
  density.name = "exp"
  lambda1 = seq(0.0000000001,50,len=200)
  lambda2 = seq(0.0000000001,50,len=200)
  n = length(lambda1)
  f1 = matrix(0,n,n)
  f2 = f1  
  f1_temp = .int.density(lambda1,rho.ref,rho.ref)-.int.density(lambda1,Umin,rho.ref)
  f2_temp = -.int.density(lambda2,Umax,rho.ref) + .int.density(lambda2,rho.ref,rho.ref)  
  for (i in c(1:n)){
    w1 = lambda2[i]/(lambda1+lambda2[i])
    w2 = lambda1[i]/(lambda1[i]+lambda2)
    f1[,i] = w1*f1_temp - w1 + alpha1
    f2[i,] = w2*f2_temp - w2 + alpha2
  }
  lines1 = contourLines(x = lambda1,y = lambda2,z = f1,levels = c(0))
  lines2 = contourLines(x = lambda1,y = lambda2,z = f2,levels = c(0))
  g1 = splinefun(lines1[[1]]$x,lines1[[1]]$y)
  g2 = splinefun(lines2[[1]]$x,lines2[[1]]$y)
  max.value = min(max(lines1[[1]]$x),max(lines2[[1]]$x))
  min.value = max(min(lines1[[1]]$x),min(lines2[[1]]$x))
  xx = seq(min.value,max.value,0.0001)
  yy = g1(xx)-g2(xx)
  g = splinefun(yy,xx)
  lam1 = g(0)
  lam2 = g1(lam1)
  drho = .rhoPDF(rho.ref=rho.ref, lambda1=lam1, lambda2=lam2, density.name=density.name)
  plot(drho$x,drho$y,xlab=expression(rho),ylab=expression(pi(rho)),type="l",ylim=c(0,drho$y[10000]+0.65))
}

.checkNumEntry <- function(entry){
  text <- entry$getText()
  if (nzchar(gsub("[0-9.-]", "", text))) {
    entry$setIconFromStock("primary", "gtk-no")
    entry$setIconTooltipText("primary","Only numbers are allowed")
  } else { 
    entry$setIconFromStock("primary", NULL)
    entry$setIconTooltipText("secondary", NULL)
  }
}

# #########################
# .makeResult <- function(x){
#   level = sort(unique(c(x$level,0.025,0.5,0.975)))
#   effect.length = length(level) + 2
#   
#   data = x$datafile
#   colnames(data) = tolower(colnames(data))
#   outdata = x$outdata
#   outpriors = x$outpriors
#   model = x$model
#   
#   nsample = x$nsample
#   model.type = x$model.type
#   var.prior = x$var.prior
#   var2.prior = x$var2.prior
#   cor.prior = x$cor.prior
#   var.par = x$var.par
#   var2.par = x$var2.par
#   cor.par = x$cor.par
#   init = c(0.01,0.01,0)
#   link=x$link
#   verbose = x$verbose
#   covariates = x$covariates
#   
#   studynames = as.character(outdata$studynames[seq(from=1,to=dim(outdata)[1],by=2)])
#   
#   variables.names = colnames(data)
#   
#   est = list()
#   est$data = data
#   est$outdata = outdata
#   est$priors.density = outpriors$density
#   est$priors.distance = outpriors$density.distance
#   
#   
#   ########## model type
#   if(x$model.type==1){
#     est$names.fitted = c("Se","Sp")
#     est$names.transf.fitted = c("1-Se","1-Sp")
#   }
#   if(x$model.type==2){
#     est$names.fitted = c("Se","1-Sp")
#     est$names.transf.fitted = c("1-Se","Sp")
#   }
#   if(x$model.type==3){
#     est$names.fitted = c("1-Se","Sp")
#     est$names.transf.fitted = c("Se","1-Sp")
#   }
#   if(x$model.type==4){
#     est$names.fitted = c("1-Se","1-Sp")
#     est$names.transf.fitted = c("Se","Sp")
#   }
#   
#   
#   n.marginal.point = dim(model$marginals.fixed[[1]])[1]
#   ########## model cpu & call
#   est$cpu.used = model[["cpu.used"]]
#   est$call = model[["call"]]
#   
#   ######## fixed
#   est$summary.fixed = model[["summary.fixed"]][,c(1:effect.length)]
#   est$marginals.fixed = model[["marginals.fixed"]]
#   
#   ######## summarized.fixed
#   names.summarized.fixed = paste("mean(logit(",est$names.fitted,"))",sep="")
#   est$summary.summarized.fixed = model[["summary.lincomb.derived"]][,c(2:(effect.length+1))]
#   est$marginals.summarized.fixed = model[["marginals.lincomb.derived"]]
#   
#   
#   ####### summarized.fitted
#   names.summarized.fitted = paste("mean(",est$names.fitted,")",sep="")
#   marginals.summarized.fitted.temp = lapply(names.summarized.fixed, function(y){
#     inla.tmarginal(function(x) .invlogit(x), model[["marginals.lincomb.derived"]][[y]], n=n.marginal.point)
#   })
#   names(marginals.summarized.fitted.temp) = names.summarized.fitted
#   # est$marginals.summarized.fitted = marginals.summarized.fitted.temp
#   
#   se.marginals = marginals.summarized.fitted.temp[[1]]
#   sp.marginals = marginals.summarized.fitted.temp[[2]]
#   suminfo1 = .summary.marginal(se.marginals,level=level)
#   suminfo2 = .summary.marginal(sp.marginals,level=level)
#   summary.summarized.fitted.temp = rbind(suminfo1, suminfo2)
#   summary.summarized.fitted.temp = as.matrix(summary.summarized.fitted.temp)
#   rownames(summary.summarized.fitted.temp) = names.summarized.fitted
#   
#   ####### summarized.transf.fitted
#   names.summarized.transf.fitted = paste("mean(",est$names.transf.fitted,")",sep="")
#   marginals.summarized.transf.fitted.temp = lapply(names.summarized.fixed, function(y){
#     inla.tmarginal(function(x) 1-.invlogit(x), model[["marginals.lincomb.derived"]][[y]], n=n.marginal.point)
#   })
#   names(marginals.summarized.transf.fitted.temp) = names.summarized.transf.fitted
#   marginals.summarized.transf.fitted = marginals.summarized.transf.fitted.temp
#   
#   se.marginals = marginals.summarized.transf.fitted[[1]]
#   sp.marginals = marginals.summarized.transf.fitted[[2]]
#   suminfo1 = .summary.marginal(se.marginals,level=level)
#   suminfo2 = .summary.marginal(sp.marginals,level=level)
#   summary.summarized.transf.fitted.temp = rbind(suminfo1, suminfo2)
#   summary.summarized.transf.fitted.temp = as.matrix(summary.summarized.transf.fitted.temp)
#   rownames(summary.summarized.transf.fitted.temp) = names.summarized.transf.fitted
#   
#   summary.summarized.fitted.total = rbind(summary.summarized.fitted.temp,summary.summarized.transf.fitted.temp)
#   est$summary.summarized.fitted = summary.summarized.fitted.total
#   
#   marginals.summarized.fitted.total = append(marginals.summarized.fitted.temp,marginals.summarized.transf.fitted)
#   est$marginals.summarized.fitted = marginals.summarized.fitted.total
#   
#   ############# hyperpar: need transformation
#   names.var = paste("var(logit(",est$names.fitted,"))",sep="")
#   names.cor = "cor(logits)"
#   tau1.marginals = model[["marginals.hyperpar"]][[1]]
#   var1.marginals = inla.tmarginal(function(x) 1/x, tau1.marginals, n=n.marginal.point)
#   tau2.marginals = model[["marginals.hyperpar"]][[2]]
#   var2.marginals = inla.tmarginal(function(x) 1/x, tau2.marginals, n=n.marginal.point)  
#   marginals.hyperpar.temp = list()
#   marginals.hyperpar.temp[[names.var[1]]] = var1.marginals
#   marginals.hyperpar.temp[[names.var[2]]] = var2.marginals
#   marginals.hyperpar.temp[[names.var[3]]] = model[["marginals.hyperpar"]][[3]]
#   est$marginals.hyperpar = marginals.hyperpar.temp
#   
#   var1.marginals = est$marginals.hyperpar[[1]]
#   var2.marginals = est$marginals.hyperpar[[2]]
#   suminfo1 = .summary.marginal(var1.marginals,level=level)
#   suminfo2 = .summary.marginal(var2.marginals,level=level)
#   
#   summary.hyperpar.temp = rbind(suminfo1, suminfo2, model[["summary.hyperpar"]][3,c(1:effect.length)])
#   summary.hyperpar.temp = as.matrix(summary.hyperpar.temp)
#   rownames(summary.hyperpar.temp) = c(names.var,names.cor)
#   est$summary.hyperpar = summary.hyperpar.temp
#   
#   #############
#   est$summarized.fixed.correlation.matrix = model[["misc"]]$lincomb.derived.correlation.matrix
#   est$summarized.fixed.covariance.matrix = model[["misc"]]$lincomb.derived.covariance.matrix
#   
#   fitted1.ind = seq(1,dim(outdata)[1],by=2)
#   fitted2.ind = seq(2,dim(outdata)[1],by=2)
#   ############# predict
#   summary.predict.temp1 = as.matrix(model[["summary.linear.predictor"]][fitted1.ind,c(1:effect.length)])
#   rownames(summary.predict.temp1) = studynames
#   summary.predict.temp2 = as.matrix(model[["summary.linear.predictor"]][fitted2.ind,c(1:effect.length)])
#   rownames(summary.predict.temp2) = studynames
#   est[[paste("summary.predict.(",est$names.fitted[1],")",sep="")]] = summary.predict.temp1
#   est[[paste("summary.predict.(",est$names.fitted[2],")",sep="")]] = summary.predict.temp2
#   
#   marginals.predict.temp1 = lapply(fitted1.ind, function(y){model[["marginals.linear.predictor"]][[y]]})
#   names(marginals.predict.temp1) = studynames
#   marginals.predict.temp2 = lapply(fitted2.ind, function(y){model[["marginals.linear.predictor"]][[y]]})
#   names(marginals.predict.temp2) = studynames
#   est[[paste("marginals.predict.(",est$names.fitted[1],")",sep="")]] = marginals.predict.temp1
#   est[[paste("marginals.predict.(",est$names.fitted[2],")",sep="")]] = marginals.predict.temp2
#   
#   ############# transform to other accurary that not fitted here directly
#   .transfunc = function(x) 1-.invlogit(x)
#   marginals.transf.fitted1 = lapply(marginals.predict.temp1, function(x){inla.tmarginal(.transfunc, x, n=n.marginal.point)})
#   marginals.transf.fitted2 = lapply(marginals.predict.temp2, function(x){inla.tmarginal(.transfunc, x, n=n.marginal.point)})
#   
#   suminfo1 = do.call(rbind,lapply(marginals.transf.fitted1, function(x) .summary.marginal(x,level=level)))
#   suminfo2 = do.call(rbind,lapply(marginals.transf.fitted2, function(x) .summary.marginal(x,level=level)))
#   
#   est[[paste("summary.fitted.(",est$names.transf.fitted[1],")",sep="")]] = suminfo1
#   est[[paste("summary.fitted.(",est$names.transf.fitted[2],")",sep="")]] = suminfo2
#   
#   ############# fitted
#   summary.fitted.temp1 = as.matrix(model[["summary.fitted.values"]][fitted1.ind,c(1:effect.length)])
#   rownames(summary.fitted.temp1) = studynames
#   summary.fitted.temp2 = as.matrix(model[["summary.fitted.values"]][fitted2.ind,c(1:effect.length)])
#   rownames(summary.fitted.temp2) = studynames
#   est[[paste("summary.fitted.(",est$names.fitted[1],")",sep="")]] = summary.fitted.temp1
#   est[[paste("summary.fitted.(",est$names.fitted[2],")",sep="")]] = summary.fitted.temp2
#   
#   marginals.fitted.temp1 = lapply(fitted1.ind, function(y){model[["marginals.fitted.values"]][[y]]})
#   names(marginals.fitted.temp1) = studynames
#   marginals.fitted.temp2 = lapply(fitted2.ind, function(y){model[["marginals.fitted.values"]][[y]]})
#   names(marginals.fitted.temp2) = studynames
#   est[[paste("marginals.fitted.(",est$names.fitted[1],")",sep="")]] = marginals.fitted.temp1
#   est[[paste("marginals.fitted.(",est$names.fitted[2],")",sep="")]] = marginals.fitted.temp2
#   
#   
#   ############# samples
#   if(is.logical(nsample)){
#     if(nsample==FALSE){
#       est$misc$sample.flag = FALSE
#     }else{
#       est$misc$sample.flag = TRUE
#       message("Argument \"nsample\" set to TRUE, we will give 5000 samples!")
#       est$misc$nsample = 5000
#     }
#   }else if(is.numeric(nsample)){
#     est$misc$sample.flag = TRUE
#     if(abs(nsample - round(nsample)) < .Machine$double.eps^0.5){
#       est$misc$nsample = nsample
#     }else{
#       message("Argument \"nsample\" should be a integer, we round the given number to interger.")
#       est$misc$nsample = round(nsample, digits = 0)
#     }
#   }else{
#     message("Argument \"nsample\" should be TRUE, FALSE or a integer, we set it to FALSE!")
#   }
#   if(est$misc$sample.flag){
#     options(warn=-1)
#     samples = inla.posterior.sample(n = nsample, model)
#     options(warn=0)
#     predictors.samples = do.call(cbind,lapply(c(1:nsample), function(x) samples[[x]]$latent[1:dim(outdata)[1]]))
#     fixed.samples = do.call(cbind,lapply(c(1:nsample), function(x){
#       length.latent = length(samples[[x]]$latent)
#       a = samples[[x]]$latent[(length.latent-dim(model$summary.fixed)[1]+1):length.latent]
#       return(a)
#     }))
#     fixed.names = rownames(model$summary.fixed)
#     rownames(fixed.samples) = fixed.names
#     ind.fitted1  = c(agrep("mu", fixed.names, max.distance=0), agrep("alpha", fixed.names, max.distance=0))
#     ind.fitted2  = c(agrep("nu", fixed.names, max.distance=0), agrep("beta", fixed.names, max.distance=0))
#     if(length(ind.fitted1)==1){
#       mean.logit.fitted1.samples = fixed.samples[ind.fitted1,]
#       mean.logit.fitted2.samples = fixed.samples[ind.fitted2,]
#     }else{
#       mean.logit.fitted1.samples = colSums(fixed.samples[ind.fitted1,])
#       mean.logit.fitted2.samples = colSums(fixed.samples[ind.fitted2,])
#     }
#     
#     mean.fitted1.samples = .invlogit(mean.logit.fitted1.samples)
#     mean.fitted2.samples = .invlogit(mean.logit.fitted2.samples)
#     
#     if(model.type==1){
#       mean.LRpos.samples = mean.fitted1.samples/(1-mean.fitted2.samples)
#       mean.LRneg.samples = (1-mean.fitted1.samples)/mean.fitted2.samples
#       mean.DOR.samples = mean.LRpos.samples/mean.LRneg.samples
#     }
#     if(model.type==2){
#       mean.LRpos.samples = mean.fitted1.samples/mean.fitted2.samples
#       mean.LRneg.samples = (1-mean.fitted1.samples)/(1-mean.fitted2.samples)
#       mean.DOR.samples = mean.LRpos.samples/mean.LRneg.samples
#     }
#     if(model.type==3){
#       mean.LRpos.samples = (1-mean.fitted1.samples)/(1-mean.fitted2.samples)
#       mean.LRneg.samples = mean.fitted1.samples/mean.fitted2.samples
#       mean.DOR.samples = mean.LRpos.samples/mean.LRneg.samples
#     }
#     if(model.type==4){
#       mean.LRpos.samples = (1-mean.fitted1.samples)/mean.fitted2.samples
#       mean.LRneg.samples = mean.fitted1.samples/(1-mean.fitted2.samples)
#       mean.DOR.samples = mean.LRpos.samples/mean.LRneg.samples
#     }
#     
#     summary.LRpos.temp = .summary.samples(mean.LRpos.samples, level=level)
#     summary.LRneg.temp = .summary.samples(mean.LRneg.samples, level=level)
#     summary.DOR.temp = .summary.samples(mean.DOR.samples, level=level)
#     
#     summary.summarized.statistics = rbind(summary.LRpos.temp, summary.LRneg.temp, summary.DOR.temp)
#     rownames(summary.summarized.statistics) = c("mean(LRpos)","mean(LRneg)","mean(DOR)")
#     est$summary.summarized.statistics = summary.summarized.statistics
#     
#     ############# fitted samples to calculate LRpos, LRneg, DOR for each study
#     fitted1.samples = .invlogit(predictors.samples[fitted1.ind,])
#     fitted2.samples = .invlogit(predictors.samples[fitted2.ind,])
#     
#     if(model.type==1){
#       LRpos.samples = fitted1.samples/(1-fitted2.samples)
#       LRneg.samples = (1-fitted1.samples)/fitted2.samples
#       DOR.samples = LRpos.samples/LRneg.samples
#     }
#     if(model.type==2){
#       LRpos.samples = fitted1.samples/fitted2.samples
#       LRneg.samples = (1-fitted1.samples)/(1-fitted2.samples)
#       DOR.samples = LRpos.samples/LRneg.samples
#     }
#     if(model.type==3){
#       LRpos.samples = (1-fitted1.samples)/(1-fitted2.samples)
#       LRneg.samples = fitted1.samples/fitted2.samples
#       DOR.samples = LRpos.samples/LRneg.samples
#     }
#     if(model.type==4){
#       LRpos.samples = (1-fitted1.samples)/fitted2.samples
#       LRneg.samples = fitted1.samples/(1-fitted2.samples)
#       DOR.samples = LRpos.samples/LRneg.samples
#     }
#     
#     summary.fitted.LRpos.temp = t(apply(LRpos.samples, 1, function(x) .summary.samples(x, level=level)))
#     rownames(summary.fitted.LRpos.temp) = studynames
#     summary.fitted.LRneg.temp = t(apply(LRneg.samples, 1, function(x) .summary.samples(x, level=level)))
#     rownames(summary.fitted.LRneg.temp) = studynames
#     summary.fitted.DOR.temp = t(apply(DOR.samples, 1, function(x) .summary.samples(x, level=level)))
#     rownames(summary.fitted.DOR.temp) = studynames
#     est$summary.fitted.LRpos = summary.fitted.LRpos.temp
#     est$summary.fitted.LRneg = summary.fitted.LRneg.temp
#     est$summary.fitted.DOR = summary.fitted.DOR.temp
#     
#     if(model.type==1){
#       est$Se.samples = fitted1.samples
#       est$Sp.samples = fitted2.samples
#       est[["mean(Se).samples"]] = mean.fitted1.samples
#       est[["mean(Sp).samples"]] = mean.fitted2.samples
#     }
#     if(model.type==2){
#       est$Se.samples = fitted1.samples
#       est$Sp.samples = 1-fitted2.samples
#       est[["mean(Se).samples"]] = mean.fitted1.samples
#       est[["mean(Sp).samples"]] = 1-mean.fitted2.samples
#     }
#     if(model.type==3){
#       est$Se.samples = 1-fitted1.samples
#       est$Sp.samples = fitted2.samples
#       est[["mean(Se).samples"]] = 1-mean.fitted1.samples
#       est[["mean(Sp).samples"]] = mean.fitted2.samples
#     }
#     if(model.type==4){
#       est$Se.samples = 1-fitted1.samples
#       est$Sp.samples = 1-fitted2.samples
#       est[["mean(Se).samples"]] = 1-mean.fitted1.samples
#       est[["mean(Sp).samples"]] = 1-mean.fitted2.samples
#     }
#   }
#   
#   ############# scores
#   est$waic = model[["waic"]]
#   est$mlik = model[["mlik"]]
#   est$cpo = model[["cpo"]]
#   est$dic = model[["dic"]]
#   
#   ############# save inla result
#   est$inla.result = model
#   ############# save general setting
#   est$misc$var.prior = var.prior
#   est$misc$cor.prior = cor.prior
#   est$misc$var.par = var.par
#   est$misc$cor.par = cor.par
#   est$misc$covariates = covariates
#   if(is.null(covariates) || covariates==FALSE){
#     est$misc$covariates.flag=FALSE
#   } else{
#     est$misc$covariates.flag=TRUE
#   }
#   if(is.null(data$modality)){
#     est$misc$modality.flag=FALSE
#   }else{
#     est$misc$modality.flag=TRUE
#     est$misc$modality.place.in.data = which(variables.names=="modality")
#     est$misc$modality.level = length(unique(data$modality))
#   }
#   est$misc$link = link
#   est$misc$model.type = model.type
#   
#   return(est)
# }

################################
####### Prior need
################################
.lprior = function(rho, rho.ref=0, lambda=1, density.name="exponential")
{
  kld = (1-rho*rho.ref)/(1-rho.ref^2)-0.5*log((1-rho^2)/(1-rho.ref^2))-1
  d = sqrt(2*kld)
  partial = (rho/(1-rho^2)-rho.ref/(1-rho.ref^2))/d
  
  if (density.name=="exp" || density.name=="exponential"){
    log.prior = log(lambda)-lambda*d + log(abs(partial))
  }
  if (density.name=="cauchy" || density.name=="half-cauchy"){
    log.prior = log(2*lambda)-log(pi)-log(d^2+lambda^2) + log(abs(partial))
  }
  if (density.name=="normal" || density.name=="half-normal"){
    log.prior = 0.5*(log(2)-log(pi))-log(lambda)-0.5*d^2/lambda^2 + log(abs(partial))
  }
  return(log.prior)
}

.rhoPDF = function(rho.ref=0, lambda1=1, lambda2=NULL, 
                   density.name="exponential", plot.flag=FALSE){
  if (is.null(lambda2) | lambda1==lambda2){
    lambda2 = lambda1
    Identity=TRUE
  } else{Identity=FALSE}
  rho.l = seq(-1,rho.ref,len=10000)
  rho.r = seq(rho.ref,1,len=10000)
  l.temp = .lprior(rho=rho.l[-length(rho.l)],rho.ref=rho.ref,lambda=lambda1,density.name=density.name)
  r.temp = .lprior(rho=rho.r[-1],rho.ref=rho.ref,lambda=lambda2,density.name=density.name)
  if (Identity==FALSE){
    sum.total = log(lambda1+lambda2)  
    if (density.name=="exp" || density.name=="exponential"){
      rate.l = log(lambda2)-sum.total # for left
      rate.r = log(lambda1)-sum.total # for right
    }
    if (density.name=="cauchy" || density.name=="half-cauchy" || density.name=="normal" || density.name=="half-normal"){
      rate.l = log(lambda1)-sum.total # for left
      rate.r = log(lambda2)-sum.total # for right
    }
    lp.l =  l.temp + rate.l
    lp.r =  r.temp + rate.r
  } else{
    lp.l = l.temp -log(2)
    lp.r = r.temp -log(2)
  }
  g1 = splinefun(rho.l[-length(rho.l)],lp.l)
  g2 = splinefun(rho.r[-1],lp.r)
  lp = c(g1(rho.l),g2(rho.r))
  pdens = exp(lp)
  rho = c(rho.l,rho.r)
  if (plot.flag){
    ylim = c(0,pdens[10000]+0.65)
    plot(rho,pdens,type="l",xlim=c(-1,1),ylim=ylim,ylab=expression(pi(rho)),xlab=expression(rho))
  }
  # ppd denotes prior probability density
  return(data.frame(x=rho,y=pdens))
}

.getRhoPDF = function(rho=0,rho.ref=0,
                      left.portion,
                      Umin,alpha1,
                      Umax,alpha2,
                      density.name="exponential",
                      plot.flag=FALSE,plot.lambda=TRUE){
  if(!missing(rho)){if(!is.numeric(rho)){stop("the interested point rho value should be numeric")}}
  if(!missing(rho.ref)){if(!is.numeric(rho.ref)){stop("the reference value of rho should be numeric and given")}}
  if(!missing(left.portion)){if(!is.numeric(left.portion)){stop("The precentile of left part of reference point should be numeric")}}
  if(!missing(Umin)){if(!is.numeric(Umin)){stop("Minimum U value should be numeric")}}
  if(!missing(Umax)){if(!is.numeric(Umax)){stop("Maximum U value should be numeric")}}
  if(!missing(alpha1)){if(!is.numeric(alpha1)){stop("alpha1 should be numeric")}}
  if(!missing(alpha2)){if(!is.numeric(alpha2)){stop("alpha2 should be numeric")}}
  if(!is.character(density.name)){stop("Distribution name should be character")}
  if(!is.logical(plot.flag)){stop("plot density should be TRUE or FALSE")}
  if(!is.logical(plot.lambda)){stop("plot lambda should be TRUE or FALSE")}
  
  # using the given values to obtain 2 lambdas
  res = .findLambdas(rho.ref=rho.ref,left.portion=left.portion,
                     Umin=Umin,alpha1=alpha1,Umax=Umax,alpha2=alpha2,
                     density.name=density.name,plot.flag=plot.lambda)
  # get the density
  data = .rhoPDF(rho.ref=rho.ref, lambda1=res$lambda1, lambda2=res$lambda2, 
                 density.name=density.name,plot.flag=plot.flag)
  func = splinefun(data$rho,data$rdv)
  if(!missing(rho)){value=func(rho)}else{value=NULL}
  result = list(density=data.frame(rho=data$rho,pdf=data$rdv),
                value=value,strategy=res$strategy,
                lambdas=c(res$lambda1,res$lambda2),
                rates=c(res$rate1,res$rate2),
                density.func=func)
  return(result)
}



.tlprior = function(theta, rho.ref=0, lambda=1, density.name="exponential")
{
  rho = .from(theta)
  kld = (1-rho*rho.ref)/(1-rho.ref^2)-0.5*log((1-rho^2)/(1-rho.ref^2))-1
  d = sqrt(2*kld)
  partial1 = (rho/(1-rho^2)-rho.ref/(1-rho.ref^2))/d
  partial2 = 2*exp(theta)/(1+exp(theta))^2
  partial = partial1*partial2
  if (density.name=="exp" || density.name=="exponential"){
    log.prior = log(lambda)-lambda*d + log(abs(partial))
  }
  if (density.name=="cauchy" || density.name=="half-cauchy"){
    log.prior = log(2*lambda)-log(pi)-log(d^2+lambda^2) + log(abs(partial))
  }
  if (density.name=="normal" || density.name=="half-normal"){
    log.prior = 0.5*(log(2)-log(pi))-log(lambda)-0.5*d^2/lambda^2 + log(abs(partial))
  }
  return(log.prior)
}

.thetaLPDF = function(rho.ref=0, lambda1=1, lambda2=NULL, 
                      density.name="exponential",ymax=3){  
  if (is.null(lambda2)==TRUE | lambda1==lambda2){
    lambda2 = lambda1
    Identity=TRUE
  } else{Identity=FALSE}
  theta.l = seq(-30,.to(rho.ref),len=10000)
  theta.r = seq(.to(rho.ref),30,len=10000)
  l.temp = .tlprior(theta=theta.l[-length(theta.l)],rho.ref=rho.ref,lambda=lambda1,density.name=density.name)
  r.temp = .tlprior(theta=theta.r[-1],rho.ref=rho.ref,lambda=lambda2,density.name=density.name)
  if (Identity==FALSE){
    sum.total = log(lambda1+lambda2)  
    if (density.name=="exp" || density.name=="exponential"){
      rate.l = log(lambda2)-sum.total # for left
      rate.r = log(lambda1)-sum.total # for right
    }
    if (density.name=="cauchy" || density.name=="half-cauchy" || density.name=="normal" || density.name=="half-normal"){
      rate.l = log(lambda1)-sum.total # for left
      rate.r = log(lambda2)-sum.total # for right
    }
    lp.l =  l.temp + rate.l
    lp.r =  r.temp + rate.r
  } else{
    lp.l = l.temp -log(2)
    lp.r = r.temp -log(2)
  }
  g1 = splinefun(theta.l[-length(theta.l)],lp.l)
  g2 = splinefun(theta.r[-1],lp.r)
  lp = c(g1(theta.l),g2(theta.r[-1]))
  #   pdens = exp(lp)
  theta = c(theta.l,theta.r[-1])
  return(data.frame(theta=theta,tldv=lp))
}

.findLambdas = function(rho.ref=0, left.portion = NULL,
                        Umin = NULL, alpha1 = NULL, Umax = NULL, alpha2 = NULL,
                        density.name="exp",
                        plot.flag=FALSE){
  if(is.null(density.name)){stop("Distributions in distance scale must be specified!!!")}
  if(!is.character(density.name)){stop("density.name must be character!!!")}
  
  density.name = tolower(density.name)
  
  # check.density = sapply(c("exp","exponential","cauchy","half-cauchy","normal","half-normal"), function(x){density.name==x})
  if (density.name!="exp" & density.name!="exponential" & density.name!="cauchy" & density.name!="half-cauchy" & density.name!="normal" & density.name!="half-normal"){
    stop("Distributions in distance scale must be set to one of the following: exponential, half-cauchy, half-normal or half-t.")
  }
  
  if (density.name=="exp" || density.name=="exponential"){
    density.name="e"
  }
  if (density.name=="cauchy" || density.name=="half-cauchy"){
    density.name="c"
  }
  if (density.name=="normal" || density.name=="half-normal"){
    density.name="n"
  }
  .erf <- function(x){2*pnorm(x*sqrt(2)) - 1}
  
  .CDFs <- function(lambda,d,density.name){
    if (density.name=="e"){
      y = 1 - exp(-lambda*d)
    }
    if (density.name=="c"){
      y = 2/pi*atan(d/lambda)
    }
    if (density.name=="n"){
      y = .erf(d/(sqrt(2)*lambda))
    }
    return(y)
  }
  
  .kld =function(rho, rho.ref) {(1-rho*rho.ref)/(1-rho.ref^2)-0.5*log((1-rho^2)/(1-rho.ref^2))-1}
  .d = function(rho, rho.ref) {sqrt(2*.kld(rho, rho.ref))}
  cc = c(0:14)
  vec = 1e-16*(10^cc)
  half_vec = 1-vec
  rho_vec = c(-half_vec,seq(-0.99,rho.ref,len=10000),seq(rho.ref,0.99,len=10000),rev(half_vec))
  dist = .d(rho_vec,rho.ref)*(rho_vec>=rho.ref)-.d(rho_vec,rho.ref)*(rho_vec<rho.ref)
  rho.to.d = splinefun(rho_vec,dist)
  d.to.rho = splinefun(dist,rho_vec)
  
  
  
  # if the left partition is not given
  if (is.null(left.portion)){
    # and also have missing values for 2 extremes, then set default values
    if (is.null(Umin) | is.null(alpha1) | is.null(Umax) | is.null(alpha2)){
      print("left.portion is not given and missing other input values, input set to default values: Umin=-0.99, alpha1=0.05, Umax=0.8, alpha2=0.2")
      Umin = -0.95
      alpha1 = 0.05
      Umax = 0.8
      alpha2 = 0.2
      strategy = 3
    }else {
      # if the 2 extremes are given, we say this is strategy 3
      strategy = 3
    }
    if(Umin>=rho.ref | Umax<=rho.ref | alpha1+alpha2>=1 | Umin==-1 | Umax==1){
      stop("The conditions must be satistied: -1<Umin<rho.ref, rho.ref<Umax<1, alpha1+alpha2<1")
    }
  } 
  else{ # if left partition is given
    if(left.portion>=1){
      stop("0<left.portion<1")
    }
    # check if left extrem is missing 
    if (is.null(Umin) | is.null(alpha1)){
      # if right extrem is also missing, set default value to left extrem 
      if (is.null(Umax) | is.null(alpha2)){
        print("missing limits probabilities, input set to default values: Umin=-0.95, alpha1=0.05")
        Umin = -0.95
        alpha1 = 0.05
        strategy = 1
        if(Umin>rho.ref | alpha1>left.portion){
          stop("The conditions must be satistied: Umin<rho.ref,alpha1<left.portion")}
      }else{ # if right exrem is not missing, we call it strategy 2
        right.portion = 1 - left.portion 
        strategy = 2
        if(Umax<rho.ref | alpha2>right.portion){
          stop("The conditions must be satistied: Umax>rho.ref,alpha2<right.portion")}
      }
    }else{ # otherwise, left extrem is not missing, (3 groups are all given),
      # then we focus on left partition and left extrem, which is strategy 1
      strategy = 1
      if(Umin>rho.ref | alpha1>left.portion){
        stop("The conditions must be satistied: Umin<rho.ref,alpha1<left.portion")}
    }
  }
  lambda1 = seq(0.0000000001,0.1,len=100)
  lambda1 = append(lambda1,seq(0.10001,0.5,len=400))
  lambda1 = append(lambda1,seq(0.5000001,10,len=2000))
  lambda2 = lambda1
  
  n = length(lambda1)
  if (strategy==1){
    # get lam1 from function f1, which is integration from Umin to rho.ref
    # we also know the relationship of lam1 and lam2 from given left partition
    # if density is exponential, left*lam1=right*lam2
    # if denisty is others (normal, cauchy), left/lam1 = right/lam2, which is left*lam2=right*lam1
    
    d_min = abs(rho.to.d(Umin))
    if (density.name=="e"){
      f2 = function(lam1,left){return(left/(1-left)*lam1)}
    }
    if (density.name=="c" || density.name=="n"){
      f2 = function(lam1,left){return((1-left)/left*lam1)}
    }
    
    f1 = left.portion*.CDFs(lambda1,d_min,density.name)-left.portion + alpha1
    g1 = splinefun(f1,lambda1)
    lam1 = g1(0)
    lam2 = f2(lam1=lam1,left=left.portion)
    if (plot.flag){
      ylim=c(-1,ceiling(lam2+1))
      xlim=c(0,ceiling(lam1+1))
      par(mar=c(5,4,4,5)+.1)
      plot(lambda1,f1,type="l",lwd=2,col="blue",xlab=expression(lambda[1]),ylab=expression(f(lambda[1])),ylim=ylim,xlim=xlim)
      lines(c(-1,lam1),c(0,0),lty=2,lwd=2,col="lightgray")
      lines(c(lam1,lam1),c(-10,0),lty=2,lwd=2,col="lightgray")
      points(lam1,0,pch="o")
      text(lam1, 0, toString(as.character(round(c(lam1,0),3))), cex=0.6, pos=3, col="black")
      par(new=TRUE)
      plot(lambda1,f2(lambda1,left.portion),type="l",lwd=2,col="black",,xaxt="n",yaxt="n",xlab="",ylab="",ylim=ylim,xlim=xlim)
      lines(c(lam1,max(lambda1)),c(lam2,lam2),lty=2,lwd=2,col="lightgray")
      lines(c(lam1,lam1),c(-10,lam2),lty=2,lwd=2,col="lightgray")
      points(lam1,lam2,pch=3)
      text(lam1, lam2, toString(as.character(round(c(lam1,lam2),3))), cex=0.6, pos=3, col="black")
      axis(4)
      mtext(expression(lambda[2]),side=4,line=3)
    }
  }
  if (strategy==2){
    # get lam2 from function f2, which is integration from Umax to 1
    # we also know the relationship of lam1 and lam2 from given left partition
    # left*lam1=right*lam2
    
    d_max = rho.to.d(Umax)
    if (density.name=="e"){
      f1 = function(lam2,left){return((1-left)/left*lam2)}
    }
    if (density.name=="c" || density.name=="n"){
      f1 = function(lam2,left){return(left/(1-left)*lam2)}
    }
    f2 = (1-left.portion)*.CDFs(lambda2,d_max,density.name) + alpha2 - (1-left.portion)
    
    g2 = splinefun(f2,lambda2)
    lam2 = g2(0)
    lam1 = f1(lam2,left.portion)
    if (plot.flag){
      ylim=c(-1,ceiling(lam1+1))
      xlim=c(0,ceiling(lam2+1))
      par(mar=c(5,4,4,5)+.1)
      plot(lambda2,f2,type="l",lwd=2,col="blue",xlab=expression(lambda[2]),ylab=expression(f(lambda[2])),ylim=ylim,xlim=xlim)
      lines(c(-1,lam2),c(0,0),lty=2,lwd=2,col="lightgray")
      lines(c(lam2,lam2),c(-10,0),lty=2,lwd=2,col="lightgray")
      points(lam2,0,pch="o")
      text(lam2, 0, toString(as.character(round(c(lam2,0),3))), cex=0.6, pos=3, col="black")
      par(new=TRUE)
      plot(lambda2,f1(lambda2,left.portion),type="l",lwd=2,col="black",xaxt="n",yaxt="n",xlab="",ylab="",ylim=ylim,xlim=xlim)
      lines(c(lam2,max(lambda2)),c(lam1,lam1),lty=2,lwd=2,col="lightgray")
      lines(c(lam2,lam2),c(-10,lam1),lty=2,lwd=2,col="lightgray")
      points(lam2,lam1,pch=3)
      text(lam2, lam1, toString(as.character(round(c(lam2,lam1),3))), cex=0.6, pos=3, col="black")
      axis(4)
      mtext(expression(lambda[1]),side=4,line=3)
    }
  }
  if (strategy==3){
    # using 2 integration: one is integrated from -1 to Umin and the other one is integrated from Umax to 1
    # get 2 integration surfaces of variables lambda1 and lambda2
    # we are interseted in the contour lines which f(lambda1,lambda2)=0
    # finding the intersection of 2 contour lines, resulting the value of lam1 and lam2
    
    d_min = abs(rho.to.d(Umin))
    d_max = rho.to.d(Umax)
    
    f1 = matrix(0,n,n)
    f2 = f1
    
    f1_temp = .CDFs(lambda1,d_min,density.name)
    f2_temp = .CDFs(lambda2,d_max,density.name) 
    
    for (i in c(1:n)){
      if (density.name=="e"){
        w1 = lambda2[i]/(lambda1+lambda2[i])
        w2 = lambda1[i]/(lambda1[i]+lambda2)
      }
      if (density.name=="c" || density.name=="n"){
        w1 = lambda1/(lambda1+lambda2[i])
        w2 = lambda2/(lambda1[i]+lambda2)
      }
      
      f1[,i] = w1*f1_temp - w1 + alpha1
      f2[i,] = w2*f2_temp - w2 + alpha2
    }
    lines1 = contourLines(x = lambda1,y = lambda2,z = f1,levels = c(0))
    lines2 = contourLines(x = lambda1,y = lambda2,z = f2,levels = c(0))
    g1 = splinefun(lines1[[1]]$x,lines1[[1]]$y)
    g2 = splinefun(lines2[[1]]$x,lines2[[1]]$y)
    max.value = min(max(lines1[[1]]$x),max(lines2[[1]]$x))
    min.value = max(min(lines1[[1]]$x),min(lines2[[1]]$x))
    xx = seq(min.value,max.value,0.0001)
    # plot(xx,g1(xx),type="l")
    # lines(xx,g2(xx),col="red")
    #     lines(lines1[[1]]$x,lines1[[1]]$y,col="red")
    #     lines(lines2[[1]]$x,lines2[[1]]$y,col="red")
    #     sum((lines1[[1]]$x==lines2[[1]]$x)*1)
    yy = g1(xx)-g2(xx)
    g = splinefun(yy,xx)
    lam1 = g(0)
    lam2 = g1(lam1)
    if(lam1<=0 || lam2<=0){
      cat("Internal Use: need to change the step size in function:findLambdas")
    }
    if(lam1>50 || lam2>50){
      cat("Internal Use: need to change the max limitation in function:findLambdas")
    }
    if(lam1<min.value || lam1>max.value){
      cat("Internal Use: restult maybe wrong, check strategy 3 in function:findLambdas")
    }
    
    if (plot.flag){
      plot(lines1[[1]]$x,lines1[[1]]$y,type="l",lwd=2,col="blue",xlab=expression(lambda[1]),ylab=expression(lambda[2]))
      lines(lines2[[1]]$x,lines2[[1]]$y,lwd=2,col="black")
      lines(c(lam1,lam1),c(0,lam2),lty=2,lwd=2,col="lightgray")
      lines(c(0,lam1),c(lam2,lam2),lty=2,lwd=2,col="lightgray")
      points(lam1,lam2,pch="o")
      text(lam1, lam2, toString(as.character(round(c(lam1,lam2),3))), cex=0.6, pos=3, col="black")
    }
  }
  # finally, function returns left portion and right portion and 2 lambdas
  if (density.name=="e"){
    return(data.frame(rate1=lam2/(lam1+lam2),rate2=lam1/(lam1+lam2),lambda1=lam1,lambda2=lam2,strategy=strategy))
  }
  if (density.name=="c" || density.name=="n"){
    return(data.frame(rate1=lam1/(lam1+lam2),rate2=lam2/(lam1+lam2),lambda1=lam1,lambda2=lam2,strategy=strategy))
  }
  
}               

