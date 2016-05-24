## Function to calculate the ecological log-likelihood w.r.t. beta (log-odds scale)

ec.log.likelihood.beta = function(beta,MM,NN,nn.var=NA){

  if(length(nn.var)==1){if(is.na(nn.var)){
    n.pos = enumerate.count(MM,NN)
    cap = 40*10^3
    if(n.pos <= cap ){ 
      nn.var = enumerate(MM,NN) 
      details = ec.log.likelihood.partial(beta,MM,NN,nn.var)
      hyblik = details$loglik

      return( list("hyblik" = hyblik) )
    }
    else{
      windowsize = ceiling( n.pos/ceiling(n.pos/cap) )
      full.windows = floor(n.pos/windowsize)
      remnants = n.pos%%windowsize

      loglik = rep(NA, full.windows + (remnants > 0) )
      for(wnd in 1:full.windows){
        startval = (wnd-1)*windowsize
        nn.var = enumerate.window(MM, NN, startval, windowsize)

        details = ec.log.likelihood.partial(beta,MM,NN,nn.var)
        loglik[wnd] = details$loglik

        rm('nn.var')
      }
      if(remnants > 0){
        nn.var = enumerate.window(MM, NN, full.windows*windowsize, remnants)
        details = ec.log.likelihood.partial(beta,MM,NN,nn.var)
        loglik[full.windows + (remnants > 0)] = details$loglik
      }
      rm('nn.var')

      baseline = which.max(loglik)
      non.base.lik = loglik[-baseline]-loglik[baseline]

      hyblik =  loglik[baseline] +
                log(1 + sum(exp(non.base.lik)))

      return(list("hyblik" = hyblik))
     }
  }}
  else{ 
    details = ec.log.likelihood.partial(beta,MM,NN,nn.var)
    hyblik = details$loglik

    return( list("hyblik" = hyblik) ) 
  }
}


ec.log.likelihood.partial = function(beta,MM,NN,nn.var){
  addbeta = beta[1] + c(0,beta[2:length(beta)])
  ebeta = exp(addbeta)

  log.like.part1 = sum(-MM*log(1+ebeta))

  n.rows = length(MM)

  columns = lchoose(MM, nn.var) + nn.var*addbeta

  eval.text = paste('columns[',1:n.rows,',]', collapse=' + ', sep='')
  log.prod = eval(parse(text = eval.text))

  baseline = which.max(log.prod)
  non.base.prod = log.prod[-baseline]-log.prod[baseline]

  log.like.part2 = log.prod[baseline] +
                   log(1 + sum(exp(non.base.prod)))
  
  log.likelihood = log.like.part1 + log.like.part2


  return(list("loglik" = log.likelihood)) 

}



## Estimation of Ecological Portion Only
log.ecportion.beta = function(beta,MM,NN,aprx=NA,nn.var=NA){
  addbeta = beta[1] + c(0,beta[2:length(beta)])
  ebeta = exp(addbeta)

  nparams = length(beta)

  px = ebeta/(1+ebeta)
  q = sum(px*MM/sum(NN))
  if(is.matrix(MM) == TRUE){q = apply(px*MM/sum(NN),2,sum)}
  else{q = sum(px*MM/sum(NN))}

  if(!is.na(aprx)){
    if(aprx=='binom'){
      ec.portion.log = dbinom(NN[2],sum(NN),q, log = TRUE)
    }
    else if(aprx=='norm'){
      px.qx.prod = ebeta/(1+ebeta)^2
      px.qx.prod[ebeta==Inf] = 0
      mu  = sum(NN)*q
      if(is.matrix(MM) == TRUE){
        sig = sqrt(apply(px.qx.prod*MM,2,sum))
        sig2 = apply(px.qx.prod*MM,2,sum)
      }
      else{
        sig = sqrt(sum(px.qx.prod*MM))
        sig2 = sum(px.qx.prod*MM)
      }

      ec.portion.log = dnorm(NN[2],mu,sig, log = TRUE)
    }


    else if(aprx=='pois'){
      lambda = sum(NN)*q
      ec.portion.log = dpois(NN[2],lambda, log = TRUE)
    }

  }
  else{
    if(length(nn.var)==1){if(is.na(nn.var)){nn.var = enumerate(MM,NN)}}
    ec.contribs = ec.log.likelihood.beta(beta,MM,NN,nn.var=nn.var)
    ec.portion.log = ec.contribs$hyblik
  }
  return(list("hyblik" = ec.portion.log))
}




## Function to calculate the hybrid likelihood, using the 
## alternative decomposition.  Likelihood calculation may be
## performed using aprximations (binomial, normal, poisson)

log.altdecomp.beta = function(beta,MM,NN,cc,aprx=NA,nn.var=NA,outvals=NA){
  mm = apply(cc,1,sum)
  nn = apply(cc,2,sum)

  MM.star = MM - mm
  NN.star = NN - nn
  N.star  = sum(NN.star)

  addbeta = beta[1] + c(0,beta[2:length(beta)])
  ebeta = exp(addbeta)

  nparams = length(beta)

  px = ebeta/(1+ebeta)
  px[ebeta==Inf] = 1

  if(is.matrix(MM) == TRUE){q = apply(px*MM.star/N.star,2,sum)}
  else{q = sum(px*MM.star/N.star)}

  if(!is.na(aprx)){
    if(aprx=='binom'){
      ec.portion.log = dbinom(NN.star[2],N.star,q, log = TRUE)
    }

    else if(aprx=='norm'){
      px.qx.prod = ebeta/(1+ebeta)^2
      px.qx.prod[ebeta==Inf] = 0
      mu  = N.star*q
      if(is.matrix(MM) == TRUE){
        sig = sqrt(apply(px.qx.prod*MM.star,2,sum))
        sig2 = apply(px.qx.prod*MM.star,2,sum)
      }
      else{
        sig = sqrt(sum(px.qx.prod*MM.star))
        sig2 = sum(px.qx.prod*MM.star)
      }
      ec.portion.log = dnorm(NN.star[2],mu,sig, log = TRUE)
    }


    else if(aprx=='pois'){
      lambda = N.star*q
      ec.portion.log = dpois(NN.star[2],lambda, log = TRUE)
    }


  }
  if(is.na(aprx)|(aprx=='NA')){
#    if(length(nn.var)==1){if(is.na(nn.var)){nn.var = enumerate(MM.star,NN.star)}}
#    ec.portion.log = ec.log.likelihood.beta(beta,MM.star,NN.star,nn.var=nn.var)
    ec.contribs = ec.log.likelihood.beta(beta,MM.star,NN.star)
    ec.portion.log = ec.contribs$hyblik
  }

  if(is.nan(ec.portion.log[1])){stop('Ecological portion calculation not recordable')}
#  else{if(ec.portion.log < -100){ec.portion.log = -100}}

  cc.portion.log.hyper = HG.log(mm,MM) - HG.log(nn,NN)

  cc.portion.log.binom = sum( lchoose(mm,cc[,2]) + cc[,2]*addbeta - mm*log(1 + ebeta) )

  cc.portion.log = cc.portion.log.hyper + cc.portion.log.binom
 
  return(list("hyblik" = ec.portion.log  + cc.portion.log))
}



## Estimation of Ecological Portion Only -- Pure ecological setting
log.ecportionECO.beta = function(beta,MM,NN,aprx=NA,nn.var=NA){
  addbeta = beta[1] + c(0,beta[2:length(beta)])
  ebeta = exp(addbeta)

  nparams = length(beta)

  px = ebeta/(1+ebeta)
  q = sum(px*MM/sum(NN))
  if(is.matrix(MM) == TRUE){q = apply(px*MM/sum(NN),2,sum)}
  else{q = sum(px*MM/sum(NN))}

  if(!is.na(aprx)){
    if(aprx=='binom'){
      ec.portion.log = dbinom(NN[2],sum(NN),q, log = TRUE)

      const1 = (NN[2]/q - NN[1]/(1-q))
      const2 = (NN[2]/q^2 + NN[1]/(1-q)^2)

      deriv.part1 = MM/sum(NN) * ebeta/(1+ebeta)^2
      if(is.matrix(MM) == TRUE){
        deriv.part2 = rbind(apply(deriv.part1,2,sum), deriv.part1[2:(nparams-1),] + rbind(deriv.part1[nparams,],deriv.part1[nparams,]) )
        hybgrad = matrix(rep(const1,(nparams-1)),byrow=TRUE,nrow=(nparams-1))*deriv.part2
      }
      else{
        deriv.part2 = c(sum(deriv.part1), deriv.part1[2:(nparams-1)] + deriv.part1[nparams] )
        hybgrad = const1*deriv.part2
      }

      if(is.matrix(MM) == TRUE){
        hybhess.part1 = apply(deriv.part2, 2, function(x){matrix(x,ncol=1)%*%matrix(x,nrow=1)})
  
        deriv2.part1 = ( MM/sum(NN) * ebeta*(1-ebeta)/(1+ebeta)^3 )
        deriv2.part2 = rbind(apply(deriv2.part1,2,sum), deriv2.part1[2:(nparams-1),] + rbind(deriv2.part1[nparams,],deriv2.part1[nparams,]) )

        hybhess.part2 = matrix(0, nrow=dim(hybhess.part1)[1], ncol = dim(hybhess.part1)[2])
        hybhess.part2[c(1,5,9),] = deriv2.part2
        hybhess.part2[c(2,3),] = hybhess.part2[c(4,7),] = deriv2.part2[-1,]
        hybhess.part2[6,] = hybhess.part2[8,] = deriv2.part1[4,]

        hybhess = -matrix(rep(const2,9),byrow=TRUE,nrow=9) * hybhess.part1 + matrix(rep(const1,9),byrow=TRUE,nrow=9)*hybhess.part2
      }
      else{
        hybhess.part1 = ( matrix(deriv.part2,ncol=1)%*%matrix(deriv.part2,nrow=1) )
  
        deriv2.part1 = ( MM/sum(NN) * ebeta*(1-ebeta)/(1+ebeta)^3 )
        deriv2.part2 = c(sum(deriv2.part1), deriv2.part1[2:(nparams-1)] + deriv2.part1[nparams] )
        hybhess.part2 = diag( deriv2.part2 )
        hybhess.part2[1,2:3] = hybhess.part2[2:3,1] = deriv2.part2[-1]
        hybhess.part2[2,3] = hybhess.part2[3,2] = deriv2.part1[4]
     
        hybhess = -const2 * hybhess.part1 + const1*hybhess.part2
      }
    }
    else if(aprx=='norm'){
      px.qx.prod = ebeta/(1+ebeta)^2
      px.qx.prod[ebeta==Inf] = 0
      mu  = sum(NN)*q
      if(is.matrix(MM) == TRUE){
        sig = sqrt(apply(px.qx.prod*MM,2,sum))
        sig2 = apply(px.qx.prod*MM,2,sum)
      }
      else{
        sig = sqrt(sum(px.qx.prod*MM))
        sig2 = sum(px.qx.prod*MM)
      }

      ec.portion.log = dnorm(NN[2],mu,sig, log = TRUE)

      d.mu = MM * ebeta/(1+ebeta)^2
      d.sig = d2.mu = MM * ebeta*(1-ebeta)/(1+ebeta)^3

      d2.sig = MM * ebeta*(1-4*ebeta+exp(2*addbeta))/(1+ebeta)^4

      const1 = (NN[2] - mu)/sig2
      const2 = (NN[2] - mu)^2/(2*sig2)^2 - 1/(2*sig2)


      if(is.matrix(MM) == TRUE){
        deriv.part1 = matrix(rep(const1,nparams),byrow=TRUE,nrow=nparams)*d.mu + 
                      matrix(rep(const2,nparams),byrow=TRUE,nrow=nparams)*d.sig
        hybgrad = rbind(apply(deriv.part1,2,sum), deriv.part1[2:(nparams-1),] + rbind(deriv.part1[nparams,],deriv.part1[nparams,]) )
      }
      else{
        deriv.part1 = const1*d.mu + const2*d.sig
        hybgrad = c(sum(deriv.part1), deriv.part1[2:(nparams-1)] + deriv.part1[nparams] )
      }

      if(is.matrix(MM) == TRUE){
        d.mu.beta = rbind(apply(d.mu,2,sum), d.mu[2:(nparams-1),] + rbind(d.mu[nparams,],d.mu[nparams,]) )
        d.sig.beta = rbind(apply(d.sig,2,sum), d.sig[2:(nparams-1),] + rbind(d.sig[nparams,],d.sig[nparams,]) )

        d2.mu.beta = matrix(0, nrow=9, ncol = dim(d.mu.beta)[2])
        d2.mu.beta[c(1,5,9),] = rbind(apply(d2.mu,2,sum), d2.mu[2:(nparams-1),] + rbind(d2.mu[nparams,],d2.mu[nparams,]) )
        d2.mu.beta[c(2,3),] = d2.mu.beta[c(4,7),] = rbind(d2.mu[2:(nparams-1),] + rbind(d2.mu[nparams,],d2.mu[nparams,]) )
        d2.mu.beta[6,] = d2.mu.beta[8,] = d2.mu[nparams,]
 
        d2.sig.beta = matrix(0, nrow=9, ncol = dim(d.sig.beta)[2]) 
        d2.sig.beta[c(1,5,9),] = rbind(apply(d2.sig,2,sum), d2.sig[2:(nparams-1),] + rbind(d2.sig[nparams,],d2.sig[nparams,]) )
        d2.sig.beta[c(2,3),] = d2.sig.beta[c(4,7),] = rbind(d2.sig[2:(nparams-1),] + rbind(d2.sig[nparams,],d2.sig[nparams,]) )
        d2.sig.beta[6,] = d2.sig.beta[8,] = d2.sig[nparams,]


        hybhess = 
        matrix(rep( (1/(2*sig2^2) - (NN[2] - mu)^2/sig2^3), 9), byrow=TRUE, nrow=9) * apply(d.sig.beta, 2, function(x){matrix(x,ncol=1)%*%matrix(x,nrow=1)}) -
        matrix(rep( ((NN[2] - mu)/sig2^2), 9), byrow=TRUE, nrow=9) * apply(rbind(d.mu.beta, d.sig.beta), 2, function(x){l = length(x); return(matrix(x[1:(l/2)],ncol=1)%*%matrix(x[(l/2+1):l],nrow=1))}) - 
        matrix(rep( ((NN[2] - mu)/sig2^2), 9), byrow=TRUE, nrow=9) * apply(rbind(d.sig.beta, d.mu.beta), 2, function(x){l = length(x); return(matrix(x[1:(l/2)],ncol=1)%*%matrix(x[(l/2+1):l],nrow=1))}) -
        matrix(rep( 1/sig2, 9), byrow=TRUE, nrow=9) * apply(d.mu.beta, 2, function(x){matrix(x,ncol=1)%*%matrix(x,nrow=1)}) + 
        matrix(rep( ((NN[2] - mu)^2/(2*sig2^2) - 1/(2*sig2)), 9), byrow=TRUE, nrow=9) * d2.sig.beta + 
        matrix(rep( (NN[2] - mu)/sig2, 9), byrow=TRUE, nrow=9) * d2.mu.beta
      }
      else{
        d.mu.beta = c(sum(d.mu), d.mu[2:(nparams-1)] + d.mu[nparams] )
        d.sig.beta = c(sum(d.sig), d.sig[2:(nparams-1)] + d.sig[nparams] )

        d2.mu.beta = diag(c(sum(d2.mu), d2.mu[2:(nparams-1)] + d2.mu[nparams] ))
        d2.mu.beta[1,2:3] = d2.mu.beta[2:3,1] = diag(d2.mu.beta)[-1]
        d2.mu.beta[2,3] = d2.mu.beta[3,2] = d2.mu[nparams]
 
        d2.sig.beta = diag(c(sum(d2.sig), d2.sig[2:(nparams-1)] + d2.sig[nparams] ))
        d2.sig.beta[1,2:3] = d2.sig.beta[2:3,1] = diag(d2.sig.beta)[-1]
        d2.sig.beta[2,3] = d2.sig.beta[3,2] = d2.sig[nparams]

        hybhess = 
        (1/(2*sig2^2) - (NN[2] - mu)^2/sig2^3) * (d.sig.beta) %*% t(d.sig.beta) -
        ((NN[2] - mu)/sig2^2) * (d.mu.beta) %*% t(d.sig.beta) - 
        ((NN[2] - mu)/sig2^2) * (d.sig.beta) %*% t(d.mu.beta) -
        1/sig2 * (d.mu.beta) %*% t(d.mu.beta) + 
        ((NN[2] - mu)^2/(2*sig2^2) - 1/(2*sig2)) * d2.sig.beta + 
        (NN[2] - mu)/sig2 * d2.mu.beta
      }

    }


    else if(aprx=='pois'){
      lambda = sum(NN)*q
      ec.portion.log = dpois(NN[2],lambda, log = TRUE)

      const1 = NN[2]/lambda -1
      const2 = NN[2]/lambda^2

      d.lam = MM * ebeta/(1+ebeta)^2
      d2.lam = MM * ebeta*(1-ebeta)/(1+ebeta)^3

      if(is.matrix(MM) == TRUE){
        hybgrad = matrix(rep(const1,(nparams-1)),byrow=TRUE,nrow=(nparams-1))*rbind(apply(d.lam,2,sum), d.lam[2:(nparams-1),] + rbind(d.lam[nparams,],d.lam[nparams,]) )
      }
      else{
        hybgrad = const1*c(sum(d.lam), d.lam[2:(nparams-1)] + d.lam[nparams] )
      }

      if(is.matrix(MM) == TRUE){
        d.lam.beta = rbind(apply(d.lam,2,sum), d.lam[2:(nparams-1),] + rbind(d.lam[nparams,],d.lam[nparams,]) )

        d2.lam.beta = matrix(0, nrow=9, ncol = dim(d.lam.beta)[2])
        d2.lam.beta[c(1,5,9),] = rbind(apply(d2.lam,2,sum), d2.lam[2:(nparams-1),] + rbind(d2.lam[nparams,],d2.lam[nparams,]) )
        d2.lam.beta[c(2,3),] = d2.lam.beta[c(4,7),] = rbind(d2.lam[2:(nparams-1),] + rbind(d2.lam[nparams,],d2.lam[nparams,]) )
        d2.lam.beta[6,] = d2.lam.beta[8,] = d2.lam[nparams,]

        hybhess = -matrix(rep(const2,9),byrow=TRUE,nrow=9) * apply(d.lam.beta, 2, function(x){matrix(x,ncol=1)%*%matrix(x,nrow=1)}) + matrix(rep(const1,9),byrow=TRUE,nrow=9) * d2.lam.beta
      }
      else{
        d.lam.beta = c(sum(d.lam), d.lam[2:(nparams-1)] + d.lam[nparams] )
        d2.lam.beta = diag(c(sum(d2.lam), d2.lam[2:(nparams-1)] + d2.lam[nparams] ))
        d2.lam.beta[1,2:3] = d2.lam.beta[2:3,1] = diag(d2.lam.beta)[-1]
        d2.lam.beta[2,3] = d2.lam.beta[3,2] = d2.lam[nparams]

        hybhess = -const2 * (d.lam.beta) %*% t(d.lam.beta) + const1 * d2.lam.beta
      }         
    }

  }
  else{
    if(length(nn.var)==1){if(is.na(nn.var)){nn.var = enumerate(MM,NN)}}
    ec.contribs = ec.log.likelihoodECO.beta(beta,MM,NN,nn.var=nn.var)
    ec.portion.log = ec.contribs$hyblik
    hybgrad = ec.contribs$hybgrad
    hybhess = ec.contribs$hybhess
  }

#  return(list("hyblik" = ec.portion.log, "hybgrad" = hybgrad))
  return(list("hyblik" = ec.portion.log, "hybgrad" = hybgrad, "hybhess" = hybhess))
}



## Function to calculate the hybrid likelihood, using the 
## alternative decomposition.  Likelihood calculation may be
## performed using aprximations (binomial, normal, poisson)

log.altdecompECO.beta = function(beta,MM,NN,cc,aprx=NA,nn.var=NA,outvals=NA){
  mm = apply(cc,1,sum)
  nn = apply(cc,2,sum)

  MM.star = MM - mm
  NN.star = NN - nn
  N.star  = sum(NN.star)

  addbeta = beta[1] + c(0,beta[2:length(beta)])
  ebeta = exp(addbeta)

  nparams = length(beta)

  px = ebeta/(1+ebeta)
  px[ebeta==Inf] = 1

  if(is.matrix(MM) == TRUE){q = apply(px*MM.star/N.star,2,sum)}
  else{q = sum(px*MM.star/N.star)}

  if(!is.na(aprx)){
    if(aprx=='binom'){
      ec.portion.log = dbinom(NN.star[2],N.star,q, log = TRUE)

      const1 = (NN.star[2]/q - NN.star[1]/(1-q))
      const2 = (NN.star[2]/q^2 + NN.star[1]/(1-q)^2)

      deriv.part1 = MM.star/sum(NN.star) * ebeta/(1+ebeta)^2
      if(is.matrix(MM) == TRUE){
        deriv.part2 = rbind(apply(deriv.part1,2,sum), deriv.part1[2:(nparams-1),] + rbind(deriv.part1[nparams,],deriv.part1[nparams,]) )
        ec.hybgrad = matrix(rep(const1,(nparams-1)),byrow=TRUE,nrow=(nparams-1))*deriv.part2
      }
      else{
        deriv.part2 = c(sum(deriv.part1), deriv.part1[2:(nparams-1)] + deriv.part1[nparams] )
        ec.hybgrad = const1*deriv.part2
      }

      if(is.matrix(MM) == TRUE){
        hybhess.part1 = apply(deriv.part2, 2, function(x){matrix(x,ncol=1)%*%matrix(x,nrow=1)})
  
        deriv2.part1 = ( MM.star/sum(NN.star) * ebeta*(1-ebeta)/(1+ebeta)^3 )
        deriv2.part2 = rbind(apply(deriv2.part1,2,sum), deriv2.part1[2:(nparams-1),] + rbind(deriv2.part1[nparams,],deriv2.part1[nparams,]) )

        hybhess.part2 = matrix(0, nrow=dim(hybhess.part1)[1], ncol = dim(hybhess.part1)[2])
        hybhess.part2[c(1,5,9),] = deriv2.part2
        hybhess.part2[c(2,3),] = hybhess.part2[c(4,7),] = deriv2.part2[-1,]
        hybhess.part2[6,] = hybhess.part2[8,] = deriv2.part1[4,]

        ec.hybhess = -matrix(rep(const2,9),byrow=TRUE,nrow=9) * hybhess.part1 + matrix(rep(const1,9),byrow=TRUE,nrow=9)*hybhess.part2
      }
      else{
        hybhess.part1 = ( matrix(deriv.part2,ncol=1)%*%matrix(deriv.part2,nrow=1) )
  
        deriv2.part1 = ( MM.star/sum(NN.star) * ebeta*(1-ebeta)/(1+ebeta)^3 )
        deriv2.part2 = c(sum(deriv2.part1), deriv2.part1[2:(nparams-1)] + deriv2.part1[nparams] )
        hybhess.part2 = diag( deriv2.part2 )
        hybhess.part2[1,2:3] = hybhess.part2[2:3,1] = deriv2.part2[-1]
        hybhess.part2[2,3] = hybhess.part2[3,2] = deriv2.part1[4]
     
        ec.hybhess = -const2 * hybhess.part1 + const1*hybhess.part2
      }
    }

    else if(aprx=='norm'){
      px.qx.prod = ebeta/(1+ebeta)^2
      px.qx.prod[ebeta==Inf] = 0
      mu  = N.star*q
      if(is.matrix(MM) == TRUE){
        sig = sqrt(apply(px.qx.prod*MM.star,2,sum))
        sig2 = apply(px.qx.prod*MM.star,2,sum)
      }
      else{
        sig = sqrt(sum(px.qx.prod*MM.star))
        sig2 = sum(px.qx.prod*MM.star)
      }
      ec.portion.log = dnorm(NN.star[2],mu,sig, log = TRUE)

      d.mu = MM.star * (ebeta/(1+ebeta)^2)
      d.sig = d2.mu = MM.star * (ebeta*(1-ebeta)/(1+ebeta)^3)

      d2.sig = MM.star * (ebeta*(1-4*ebeta+exp(2*addbeta))/(1+ebeta)^4)

      const1 = (NN.star[2] - mu)/sig2
      const2 = (NN.star[2] - mu)^2/(2*sig2)^2 - 1/(2*sig2)

      if(is.matrix(MM) == TRUE){
        deriv.part1 = matrix(rep(const1,nparams),byrow=TRUE,nrow=nparams)*d.mu + 
                      matrix(rep(const2,nparams),byrow=TRUE,nrow=nparams)*d.sig
        ec.hybgrad = rbind(apply(deriv.part1,2,sum), deriv.part1[2:(nparams-1),] + rbind(deriv.part1[nparams,],deriv.part1[nparams,]) )
      }
      else{
        deriv.part1 = const1*d.mu + const2*d.sig
        ec.hybgrad = c(sum(deriv.part1), deriv.part1[2:(nparams-1)] + deriv.part1[nparams] )
      }



      if(is.matrix(MM) == TRUE){
        d.mu.beta = rbind(apply(d.mu,2,sum), d.mu[2:(nparams-1),] + rbind(d.mu[nparams,],d.mu[nparams,]) )
        d.sig.beta = rbind(apply(d.sig,2,sum), d.sig[2:(nparams-1),] + rbind(d.sig[nparams,],d.sig[nparams,]) )

        d2.mu.beta = matrix(0, nrow=9, ncol = dim(d.mu.beta)[2])
        d2.mu.beta[c(1,5,9),] = rbind(apply(d2.mu,2,sum), d2.mu[2:(nparams-1),] + rbind(d2.mu[nparams,],d2.mu[nparams,]) )
        d2.mu.beta[c(2,3),] = d2.mu.beta[c(4,7),] = rbind(d2.mu[2:(nparams-1),] + rbind(d2.mu[nparams,],d2.mu[nparams,]) )
        d2.mu.beta[6,] = d2.mu.beta[8,] = d2.mu[nparams,]
 
        d2.sig.beta = matrix(0, nrow=9, ncol = dim(d.sig.beta)[2]) 
        d2.sig.beta[c(1,5,9),] = rbind(apply(d2.sig,2,sum), d2.sig[2:(nparams-1),] + rbind(d2.sig[nparams,],d2.sig[nparams,]) )
        d2.sig.beta[c(2,3),] = d2.sig.beta[c(4,7),] = rbind(d2.sig[2:(nparams-1),] + rbind(d2.sig[nparams,],d2.sig[nparams,]) )
        d2.sig.beta[6,] = d2.sig.beta[8,] = d2.sig[nparams,]


        ec.hybhess = 
        matrix(rep( (1/(2*sig2^2) - (NN.star[2] - mu)^2/sig2^3), 9), byrow=TRUE, nrow=9) * apply(d.sig.beta, 2, function(x){matrix(x,ncol=1)%*%matrix(x,nrow=1)}) -
        matrix(rep( ((NN.star[2] - mu)/sig2^2), 9), byrow=TRUE, nrow=9) * apply(rbind(d.mu.beta, d.sig.beta), 2, function(x){l = length(x); return(matrix(x[1:(l/2)],ncol=1)%*%matrix(x[(l/2+1):l],nrow=1))}) - 
        matrix(rep( ((NN.star[2] - mu)/sig2^2), 9), byrow=TRUE, nrow=9) * apply(rbind(d.sig.beta, d.mu.beta), 2, function(x){l = length(x); return(matrix(x[1:(l/2)],ncol=1)%*%matrix(x[(l/2+1):l],nrow=1))}) -
        matrix(rep( 1/sig2, 9), byrow=TRUE, nrow=9) * apply(d.mu.beta, 2, function(x){matrix(x,ncol=1)%*%matrix(x,nrow=1)}) + 
        matrix(rep( ((NN.star[2] - mu)^2/(2*sig2^2) - 1/(2*sig2)), 9), byrow=TRUE, nrow=9) * d2.sig.beta + 
        matrix(rep( (NN.star[2] - mu)/sig2, 9), byrow=TRUE, nrow=9) * d2.mu.beta
      }
      else{
        d.mu.beta = c(sum(d.mu), d.mu[2:(nparams-1)] + d.mu[nparams] )
        d.sig.beta = c(sum(d.sig), d.sig[2:(nparams-1)] + d.sig[nparams] )

        d2.mu.beta = diag(c(sum(d2.mu), d2.mu[2:(nparams-1)] + d2.mu[nparams] ))
        d2.mu.beta[1,2:3] = d2.mu.beta[2:3,1] = diag(d2.mu.beta)[-1]
        d2.mu.beta[2,3] = d2.mu.beta[3,2] = d2.mu[nparams]
 
        d2.sig.beta = diag(c(sum(d2.sig), d2.sig[2:(nparams-1)] + d2.sig[nparams] ))
        d2.sig.beta[1,2:3] = d2.sig.beta[2:3,1] = diag(d2.sig.beta)[-1]
        d2.sig.beta[2,3] = d2.sig.beta[3,2] = d2.sig[nparams]

        ec.hybhess = 
        (1/(2*sig2^2) - (NN.star[2] - mu)^2/sig2^3) * (d.sig.beta) %*% t(d.sig.beta) -
        ((NN.star[2] - mu)/sig2^2) * (d.mu.beta) %*% t(d.sig.beta) - 
        ((NN.star[2] - mu)/sig2^2) * (d.sig.beta) %*% t(d.mu.beta) -
        1/sig2 * (d.mu.beta) %*% t(d.mu.beta) + 
        ((NN.star[2] - mu)^2/(2*sig2^2) - 1/(2*sig2)) * d2.sig.beta + 
        (NN.star[2] - mu)/sig2 * d2.mu.beta
      }
    }


    else if(aprx=='pois'){
      lambda = N.star*q
      ec.portion.log = dpois(NN.star[2],lambda, log = TRUE)

      const1 = NN.star[2]/lambda -1
      const2 = NN.star[2]/lambda^2

      d.lam = MM.star * ebeta/(1+ebeta)^2
      d2.lam = MM.star * ebeta*(1-ebeta)/(1+ebeta)^3

      if(is.matrix(MM) == TRUE){
        ec.hybgrad = matrix(rep(const1,(nparams-1)),byrow=TRUE,nrow=(nparams-1))*rbind(apply(d.lam,2,sum), d.lam[2:(nparams-1),] + rbind(d.lam[nparams,],d.lam[nparams,]) )
      }
      else{
        ec.hybgrad = const1*c(sum(d.lam), d.lam[2:(nparams-1)] + d.lam[nparams] )
      }

      if(is.matrix(MM) == TRUE){
        d.lam.beta = rbind(apply(d.lam,2,sum), d.lam[2:(nparams-1),] + rbind(d.lam[nparams,],d.lam[nparams,]) )

        d2.lam.beta = matrix(0, nrow=9, ncol = dim(d.lam.beta)[2])
        d2.lam.beta[c(1,5,9),] = rbind(apply(d2.lam,2,sum), d2.lam[2:(nparams-1),] + rbind(d2.lam[nparams,],d2.lam[nparams,]) )
        d2.lam.beta[c(2,3),] = d2.lam.beta[c(4,7),] = rbind(d2.lam[2:(nparams-1),] + rbind(d2.lam[nparams,],d2.lam[nparams,]) )
        d2.lam.beta[6,] = d2.lam.beta[8,] = d2.lam[nparams,]

        ec.hybhess = -matrix(rep(const2,9),byrow=TRUE,nrow=9) * apply(d.lam.beta, 2, function(x){matrix(x,ncol=1)%*%matrix(x,nrow=1)}) + matrix(rep(const1,9),byrow=TRUE,nrow=9) * d2.lam.beta
      }
      else{
        d.lam.beta = c(sum(d.lam), d.lam[2:(nparams-1)] + d.lam[nparams] )
        d2.lam.beta = diag(c(sum(d2.lam), d2.lam[2:(nparams-1)] + d2.lam[nparams] ))
        d2.lam.beta[1,2:3] = d2.lam.beta[2:3,1] = diag(d2.lam.beta)[-1]
        d2.lam.beta[2,3] = d2.lam.beta[3,2] = d2.lam[nparams]

        ec.hybhess = -const2 * (d.lam.beta) %*% t(d.lam.beta) + const1 * d2.lam.beta
      }    
    }


  }
  if(is.na(aprx)|(aprx=='NA')){
#    if(length(nn.var)==1){if(is.na(nn.var)){nn.var = enumerate(MM.star,NN.star)}}
#    ec.portion.log = ec.log.likelihoodECO.beta(beta,MM.star,NN.star,nn.var=nn.var)
    ec.contribs = ec.log.likelihoodECO.beta(beta,MM.star,NN.star)
    ec.portion.log = ec.contribs$hyblik
    ec.hybgrad = ec.contribs$hybgrad
    ec.hybhess = ec.contribs$hybhess
  }

  if(is.nan(ec.portion.log[1])){stop('Ecological portion calculation not recordable')}
#  else{if(ec.portion.log < -100){ec.portion.log = -100}}

  cc.portion.log.hyper = HG.log(mm,MM) - HG.log(nn,NN)

  cc.portion.log.binom = sum( lchoose(mm,cc[,2]) + cc[,2]*addbeta - mm*log(1 + ebeta) )

  cc.portion.log = cc.portion.log.hyper + cc.portion.log.binom
 
  grad = cc[,2] - mm*ebeta/(1 + ebeta) 

  cc.db0 = sum(grad)
  cc.dbr = grad[2]+grad[4]
  cc.dbs = grad[3]+grad[4]

  cc.hybgrad = c(cc.db0, cc.dbr, cc.dbs)

  deriv2 = -mm*ebeta/(1+ebeta)^2

  cc.hybhess = diag( c(sum(deriv2), deriv2[2] + deriv2[4], deriv2[3] + deriv2[4]) )
  cc.hybhess[1,2] = cc.hybhess[2,1] = deriv2[2] + deriv2[4]
  cc.hybhess[1,3] = cc.hybhess[3,1] = deriv2[3] + deriv2[4]
  cc.hybhess[2,3] = cc.hybhess[3,2] = deriv2[4]

  if(is.matrix(MM) == TRUE){
    cc.hybhess = matrix(rep(as.vector(cc.hybhess), dim(ec.hybhess)[2]), ncol = dim(ec.hybhess)[2])
  }
 
#  return(list("hyblik" = ec.portion.log  + cc.portion.log, "hybgrad" = ec.hybgrad + cc.hybgrad))
  return(list("hyblik" = ec.portion.log  + cc.portion.log, "hybgrad" = ec.hybgrad + cc.hybgrad, "hybhess" = ec.hybhess + cc.hybhess))

#  return(list("hyblik" = cc.portion.log, "hybgrad" = cc.hybgrad, "hybhess" = cc.hybhess))
#  return(list("hyblik" = ec.portion.log, "hybgrad" = ec.hybgrad, "hybhess" = ec.hybhess))

}

## Function to calculate the ecological log-likelihood w.r.t. beta (log-odds scale)

ec.log.likelihoodECO.beta = function(beta,MM,NN,nn.var=NA){

  if(length(nn.var)==1){if(is.na(nn.var)){
    n.pos = enumerate.count(MM,NN)
    cap = 40*10^3
    if(n.pos <= cap ){ 
      nn.var = enumerate(MM,NN) 
      details = ec.log.likelihood.partial(beta,MM,NN,nn.var)
      hyblik = details$loglik
      hybgrad = details$grad / exp(hyblik)
   
      hybhess = -hybgrad %*% t(hybgrad) + details$hess/exp(hyblik)

      return( list("hyblik" = hyblik, "hybgrad" = hybgrad, "hybhess" = hybhess) )
    }
    else{
      windowsize = ceiling( n.pos/ceiling(n.pos/cap) )
      full.windows = floor(n.pos/windowsize)
      remnants = n.pos%%windowsize

      loglik = rep(NA, full.windows + (remnants > 0) )
      hybgrad.unscaled = rep(0, 3)
      hybhess.part2 = matrix(0,ncol=3,nrow=3)
      for(wnd in 1:full.windows){
        startval = (wnd-1)*windowsize
        nn.var = enumerate.window(MM, NN, startval, windowsize)

        details = ec.log.likelihood.partial(beta,MM,NN,nn.var)
        loglik[wnd] = details$loglik

        hybgrad.unscaled = hybgrad.unscaled + details$grad
   
        hybhess.part2 = hybhess.part2 + details$hess

        rm('nn.var')
      }
      if(remnants > 0){
        nn.var = enumerate.window(MM, NN, full.windows*windowsize, remnants)
        details = ec.log.likelihood.partial(beta,MM,NN,nn.var)
        loglik[full.windows + (remnants > 0)] = details$loglik

        hybgrad.unscaled = hybgrad.unscaled + details$grad
   
        hybhess.part2 = hybhess.part2 + details$hess
      }
      rm('nn.var')

      baseline = which.max(loglik)
      non.base.lik = loglik[-baseline]-loglik[baseline]

      hyblik =  loglik[baseline] +
                log(1 + sum(exp(non.base.lik)))

      hybgrad = hybgrad.unscaled / exp(hyblik)
   
      hybhess = -hybgrad %*% t(hybgrad) + hybhess.part2/exp(hyblik)

      return(list("hyblik" = hyblik, "hybgrad" = hybgrad, "hybhess" = hybhess))
     }
  }}
  else{ 
    details = ec.log.likelihood.partial(beta,MM,NN,nn.var)
    hyblik = details$loglik
    hybgrad = details$grad / exp(hyblik)
   
    hybhess = -hybgrad %*% t(hybgrad) + details$hess/exp(hyblik)

    return( list("hyblik" = hyblik, "hybgrad" = hybgrad, "hybhess" = hybhess) ) 
  }
}


ec.log.likelihood.partial = function(beta,MM,NN,nn.var){
  addbeta = beta[1] + c(0,beta[2:length(beta)])
  ebeta = exp(addbeta)

  log.like.part1 = sum(-MM*log(1+ebeta))

  n.rows = length(MM)

  columns = lchoose(MM, nn.var) + nn.var*addbeta

  eval.text = paste('columns[',1:n.rows,',]', collapse=' + ', sep='')
  log.prod = eval(parse(text = eval.text))

  baseline = which.max(log.prod)
  non.base.prod = log.prod[-baseline]-log.prod[baseline]

  log.like.part2 = log.prod[baseline] +
                   log(1 + sum(exp(non.base.prod)))
  
  log.likelihood = log.like.part1 + log.like.part2


  gradient.part1 = nn.var - MM*ebeta/(1+ebeta)
  dll.b0 = gradient.part1[1,] + gradient.part1[2,] + gradient.part1[3,] + gradient.part1[4,]
  dll.br = gradient.part1[2,] + gradient.part1[4,]
  dll.bs = gradient.part1[3,] + gradient.part1[4,]

  db0.unscaled = sum(exp(log.like.part1)*exp(log.prod)*dll.b0)
  dbr.unscaled = sum(exp(log.like.part1)*exp(log.prod)*dll.br)
  dbs.unscaled = sum(exp(log.like.part1)*exp(log.prod)*dll.bs)

  gradient = c(db0.unscaled, dbr.unscaled, dbs.unscaled)

  deriv2 = -MM*ebeta/(1+ebeta)^2
  hybhess.part2 = diag( c(sum(deriv2), deriv2[2] + deriv2[4], deriv2[3] + deriv2[4]) )
  hybhess.part2[1,2] = hybhess.part2[2,1] = deriv2[2] + deriv2[4]
  hybhess.part2[1,3] = hybhess.part2[3,1] = deriv2[3] + deriv2[4]
  hybhess.part2[2,3] = hybhess.part2[3,2] = deriv2[4]

  pcov.b0 = sum(exp(log.like.part1)*exp(log.prod)* (dll.b0*dll.b0 + hybhess.part2[1,1]) )
  pcov.br = sum(exp(log.like.part1)*exp(log.prod)* (dll.br*dll.br + hybhess.part2[2,2]) )
  pcov.bs = sum(exp(log.like.part1)*exp(log.prod)* (dll.bs*dll.bs + hybhess.part2[3,3]) )
  pcov.b0r = sum(exp(log.like.part1)*exp(log.prod)* (dll.b0*dll.br + hybhess.part2[1,2]) )
  pcov.b0s = sum(exp(log.like.part1)*exp(log.prod)* (dll.b0*dll.bs + hybhess.part2[1,3]) )
  pcov.brs = sum(exp(log.like.part1)*exp(log.prod)* (dll.br*dll.bs + hybhess.part2[2,3]) )

  partial.hessian = matrix(c(pcov.b0, pcov.b0r, pcov.b0s, pcov.b0r, pcov.br, 
    pcov.brs, pcov.b0s, pcov.brs, pcov.bs), byrow=TRUE, ncol=3)


  return(list("loglik" = log.likelihood, "grad" = gradient, "hess" = partial.hessian)) 

}




