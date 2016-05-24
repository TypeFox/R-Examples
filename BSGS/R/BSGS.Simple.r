#=======================================================================================================================================#
# Bayesian Sparse Group Selection                                                                                                       #
# Date: 01/15/2014                                                                                                                      #
# Maintainer : Kuo-Jung Lee & Ray-Bing Chen                                                                                             #                                                            #
# Description: Perform Bayesian sparse group selection                                                                                  #
#=======================================================================================================================================#
Compute.Z = function(k.vector, Y, X, beta.value, rho.value, tau2.value, sigma2.value)
{
  if(sum(k.vector)!=0){
    R.m = Y - X[,k.vector!=1]%*%cbind(beta.value[k.vector!=1])
    X.m = cbind(X[,k.vector==1])
    beta.m = cbind(beta.value[k.vector==1])

    #print(1/tau2.value[k.vector==1])
    T.m.inv = diag(1/tau2.value[k.vector==1], sum(k.vector))
    Q.value = t(X.m)%*%X.m/sigma2.value + T.m.inv
    Q.value.inv = solve(Q.value)
    #print((0.5*t(R.m)%*%X.m%*%Q.value.inv%*%t(X.m)%*%R.m / sigma2.value^2))
    Z.m = exp(0.5*t(R.m)%*%X.m%*%Q.value.inv%*%t(X.m)%*%R.m / sigma2.value^2)/(sqrt(prod(tau2.value[k.vector==1]))*sqrt(det(Q.value)))*
          prod(rho.value[k.vector==0])*prod((1-rho.value[k.vector==1]))
    #print(Z.m)
  }
  else
    Z.m =0
}

GroupZ.MH = function(Y, X, beta.value, rho, tau2, sigma2, theta)
{
  AllPossible= as.matrix(expand.grid(rep(list(0:1), length(tau2))))
  Z.sum = sum(as.vector(apply(AllPossible, 1, Compute.Z, Y, X, beta.value, rho, tau2, sigma2)))
  Z.m = prod(rho) + Z.sum
  prob.MH = (1-theta)/(1-theta + theta/Z.m)
  eta = ifelse(runif(1)<prob.MH, 1, 0)
}




BSGS.Simple = function(Y, X, Group.Index, r.value, eta.value, beta.value, tau2.value, rho.value, theta.value, sigma2.value, nu, lambda, Num.of.Iter.Inside.CompWise, Num.Of.Iteration=1000, MCSE.Sigma2.Given = 0.01)
{
  num.of.obs = length(Y)
  num.of.covariates = ncol(X)
  group.inf = (count(Group.Index))
  num.of.groups = dim(group.inf)[1]
  group.lables = group.inf[, 1]
  num.of.var.in.groups = group.inf[, 2]


  r.samples = eta.samples = sigma2.samples = beta.samples = numeric(0)
  #r.samples = matrix(0, Num.Of.Iteration, num.of.covariates)
  #eta.samples = matrix(0, Num.Of.Iteration, num.of.groups)
  #sigma2.samples = rep(0, Num.Of.Iteration)
  #beta.samples =  matrix(0, Num.Of.Iteration, num.of.covariates)

  r.sample = r.samples = r.value# r.true = as.numeric(beta.value!=0)
  beta.sample = beta.samples = beta.value
  sigma2.sample = sigma2.samples = sigma2.value
  eta.sample = eta.samples = eta.value #rep(0, num.of.group.var) 

  colnames(X) = names(r.sample) = names(beta.sample)= names(tau2.value) = names(rho.value) = Group.Index
  names(eta.sample) = names(theta.value) =  group.lables
  beta.sample = cbind(beta.sample)

  iter = 1
  start.run = proc.time()
  while(iter < Num.Of.Iteration){

  ## Step I
    group.selection = sample(group.lables, 1)

    Res = Y- (X[, Group.Index != group.selection]) %*% cbind(beta.sample[Group.Index != group.selection])

  ## Step II & III

    X.cws = cbind(X[, Group.Index == group.selection])
    beta.sample.cws = cbind(beta.sample[Group.Index == group.selection, ])
    r.sample.cws = r.sample[Group.Index == group.selection]
    rho.value.cws = rho.value[Group.Index == group.selection]
    theta.value.cws = theta.value[names(eta.sample) == group.selection]
    tau2.value.cws = tau2.value[Group.Index == group.selection]

    eta.cand = GroupZ.MH(Res, X.cws, beta.sample.cws, rho.value.cws, tau2.value.cws, sigma2.sample, theta.value.cws)
    if(eta.cand == 1){
      beta.r.cand = CompWiseGibbs(Res, X.cws, beta.sample.cws, r.sample.cws, tau2.value.cws, rho.value.cws, sigma2.sample, Num.of.Iter.Inside.CompWise)
      beta.sample[Group.Index == group.selection, ] = beta.r.cand$beta.sample #beta.true[group.index == group.selection]#
      r.sample[Group.Index == group.selection] = beta.r.cand$r.sample#r.true[group.index == group.selection]#
    }
    else{
      beta.sample[Group.Index == group.selection, ] = 0 #beta.true[group.index == group.selection]#
      r.sample[Group.Index == group.selection] = 0 #r.true[group.index == group.selection]
    }


   ## Step IV

    SSE = sum( (Y- X%*% beta.sample)^2 )
    sigma2.sample = 1/rgamma(1, ((num.of.obs+nu)/2), (SSE+nu*lambda)/2)#sigma2.true #
		
		eta.sample[names(eta.sample) == group.selection] = eta.cand
    
    
    r.samples = rbind(r.samples, r.sample)
    beta.samples = cbind(beta.samples, beta.sample)
    sigma2.samples = c(sigma2.samples, sigma2.sample)
    eta.samples = rbind(eta.samples, eta.sample)

    iter = iter + 1

    if(iter==Num.Of.Iteration)
      if(bm(sigma2.samples)$se > MCSE.Sigma2.Given){
        Num.Of.Iteration = Num.Of.Iteration + 100
        cat("Sigma2 = ", bm(sigma2.samples)$est, "\tMCSE = ", bm(sigma2.samples)$se, "\tNumber of Iterations = ", Num.Of.Iteration, "\n")
      }

  }
  end.run = proc.time()
  list(beta.samples = beta.samples, eta.samples = eta.samples, r.samples = r.samples, sigma2.samples = sigma2.samples, Iteration = Num.Of.Iteration, TimeElapsed = (end.run-start.run))
}


