genPBONN <- function (n, cmat.star, no.pois=0, no.bin=0, no.ord=0, no.nonn=0, 
                      pois.list=list(), bin.list=list(), ord.list=list(), is.ord.list.cum=FALSE, nonn.list=list()) {
  no.total = no.pois+no.bin+no.ord+no.nonn
  if(no.total <= 0) stop("There must be at least one variable!")
  nonn.list = lapply(nonn.list, .addMeanVar)
  .check.cor.mat(no.total, cmat.star)
  if(no.ord > 0 && !is.ord.list.cum) {
    ordbin.cum.list = lapply(c(bin.list, ord.list), .getCum)
  } else {
    ordbin.cum.list =  c(lapply(c(bin.list), .getCum), ord.list)
  }
  nonn.bcd.list = lapply(nonn.list, .getBCD)
  all.list = c(pois.list,ordbin.cum.list,nonn.bcd.list)
  z = mvrnorm(n, rep(0, no.total), cmat.star)
  x = z
  for(k in 1:no.total) {
    arg = all.list[[k]]
    gen = .n.to.nonn
    if(k <= no.pois) gen = .n.to.pois
    else if(k <= no.pois+no.bin+no.ord) gen = .n.to.ord
    x[,k] = sapply(z[,k], gen, arg) 
  }
  return(x)  
}

.addMeanVar <- function(mvector) {
  if(length(mvector) > 4) stop("Only the first four moments considered for continuous variables!")
  if(length(mvector) < 2) stop("At least two moments must be specified!")
  if(length(mvector)==4) return(mvector)
  if(length(mvector)==3) return(c(0,mvector))
  return(c(0,1,mvector))
}

check.params <- function(no.pois=0, no.ordbin=0, no.nonn=0, pois.list=list(), ordbin.list=list(), nonn.list=list()) {
  check.pois <- function(no.pois, pois.list) {
    if (length(pois.list) != no.pois) {
      stop("Dimension of the Poisson list does not match the number of Poisson variables!\n")
    }
    if (!all(pois.list > 0) ) {
      stop("Poisson rate parameters must be greater than zero!\n")
    }
    return(TRUE)
  }
  check.nonn <- function(no.nonn, nonn.list) {
    if (length(nonn.list) != no.nonn) {
      stop("Dimension of the continuous list does not match the number of variables!\n")
    }
    nonn.list = lapply(nonn.list,.addMeanVar)
    if(!all(unlist(lapply(nonn.list, .checkMoments)))) {
      stop("Non-normal skewness and excess kurtosis values must follow inequality: \nexcess kurtosis >= skew^2-2, and variance must be greater than 0!\n")
    }
  }
  check.ord <- function(no.ord, ord.list) {
    if (length(ord.list) != no.ord) {
      stop("Dimension of the ordinal and/or binary list does not match \nthe number of ordinal and/or binary variables!\n")
    }
    if(no.ord > 0) {
#Does not apply when list is cumulative
#      if(!all(lapply(ord.list, sum)==1)) {
#        stop("Ordinal rates must sum to 1!\n")
#      }
      if(!all(unlist(lapply(ord.list, min) >= 0)) || !all(unlist(lapply(ord.list,max)<=1))) {
        stop("Ordinal probabilities must be between 0 and 1!\n")
      }
      if(any(unlist(lapply(ord.list, is.unsorted)))) {
        stop("Each entry of list must be cumulative! \ncumulative ordinal probabilities must be nondecreasing!")
      }
    }
  }
  check.pois(no.pois, pois.list)
  check.ord(no.ordbin, ordbin.list)
  check.nonn(no.nonn, nonn.list)
}

.check.cor.mat <- function(no.total, cor.mat) {
  if(nrow(cor.mat)!= no.total) {
    stop("Dimension of cmat.star and number of variables do not match!\n")
  }
  if(!isSymmetric(cor.mat)) {
    stop("Correlation matrix is not square symmetric!\n")
  }
  if (!is.positive.definite(cor.mat)) {
    stop("Correlation matrix is not positive definite!\n")
  }
  if (!all(cor.mat <= 1)) stop("Correlation values cannot be greater than 1! \n")
  if (!all(cor.mat >= -1)) stop("Correlation values cannot be less than -1! \n")
  if (!all(diag(cor.mat) == 1)) stop("All diagonal elements of the correlation matrix must be 1! \n")  
}

lower.upper.cors <- function(no.pois=0, no.bin=0, no.ord=0, no.nonn=0, pois.list=list(), 
                             bin.list=list(), ord.list=list(), is.ord.list.cum=FALSE, nonn.list=list()) {
  no.total = no.pois+no.bin+no.ord+no.nonn
  if(no.total <= 0) stop("There must be at least one variable!")
  if(no.ord > 0 && !is.ord.list.cum) {
    ordbin.cum.list = lapply(c(bin.list, ord.list), .getCum)
  } else {
    ordbin.cum.list =  c(lapply(bin.list, .getCum), ord.list)
  }
  nonn.list = lapply(nonn.list, .addMeanVar)
  check.params(no.pois, no.bin+no.ord, no.nonn, pois.list, ordbin.cum.list, nonn.list)
  nonn.bcd.list = lapply(nonn.list, .getBCD)
  all.list = c(pois.list,ordbin.cum.list,nonn.bcd.list)
  
  allcombos = combn(no.total, 2)
  maxMins = apply(allcombos, 2, .getMaxMin, no.pois, no.bin, no.ord, no.nonn, all.list)
  mins = .make.sym.matrix(no.total, maxMins[2,])
  maxes = .make.sym.matrix(no.total, maxMins[1,])
  return(list(min=mins, max=maxes))
}

validate.cor.mat <- function(cor.mat, no.pois=0, no.bin=0, no.ord=0, no.nonn=0, pois.list=list(), 
                             bin.list=list(), ord.list=list(), is.ord.list.cum=FALSE, nonn.list=list()) {
  no.total = no.pois+no.bin+no.ord+no.nonn
  if(no.total <= 0) stop("There must be at least one variable!")
  .check.cor.mat(no.total, cor.mat)
  if(no.ord > 0 && !is.ord.list.cum) {
    ordbin.cum.list = lapply(c(bin.list, ord.list), .getCum)
  } else {
    ordbin.cum.list =  c(lapply(c(bin.list), .getCum), ord.list)
  }
  nonn.list = lapply(nonn.list, .addMeanVar)
  check.params(no.pois, no.bin+no.ord, no.nonn, pois.list, ordbin.cum.list, nonn.list)
  .check.cor.mat(no.total, cor.mat)
  nonn.bcd.list = lapply(nonn.list, .getBCD)
  all.list = c(pois.list,ordbin.cum.list,nonn.bcd.list)
  min.max.list = lower.upper.cors(no.pois, no.bin, no.ord, no.nonn, pois.list, bin.list, ord.list, is.ord.list.cum, nonn.list)
  tf = min.max.list[[1]] <= cor.mat & cor.mat <= min.max.list[[2]]
  if(all(tf)) return(TRUE)
  entries = which(!tf)-1
  cols = entries %/% nrow(tf)+1
  rows = entries %% ncol(tf)+1
  min = min.max.list[[1]][entries+1]
  max = min.max.list[[2]][entries+1]
  for(k in 1:sum(!tf)) {
    cat(paste(paste0("Cell [", rows[k], ",", cols[k], "] is invalid, must be between"), round(min[k],2), "and", round(max[k],2),"\n"))
  }
  return(FALSE)
}

find.cor.mat.star <- function(cor.mat, no.pois=0, no.bin=0, no.ord=0, no.nonn=0, pois.list=list(), 
                              bin.list=list(), ord.list=list(), is.ord.list.cum=FALSE, nonn.list=list()) {
  no.total = no.pois+no.bin+no.ord+no.nonn
  if(no.total <= 0) stop("There must be at least one variable!")
  if(validate.cor.mat(cor.mat, no.pois, no.bin, no.ord, no.nonn, pois.list, bin.list, ord.list, is.ord.list.cum, nonn.list) != TRUE) {
    stop("Cor.mat is not valid!")
  }
  if(no.ord > 0 && !is.ord.list.cum) {
    ordbin.cum.list = lapply(c(bin.list, ord.list), .getCum)
  } else {
    ordbin.cum.list =  c(lapply(c(bin.list), .getCum), ord.list)
  }
  nonn.list = lapply(nonn.list, .addMeanVar)
  nonn.bcd.list = lapply(nonn.list, .getBCD)
  all.list = c(pois.list,ordbin.cum.list,nonn.bcd.list)
  cor.mat.star = cor.mat
  if(no.pois > 1) {
    pp.combo = combn(no.pois, 2)
    for(k in 1:ncol(pp.combo)) {
      i=pp.combo[1,k]
      j=pp.combo[2,k]
      cor.mat.star[i,j] = .cor.P.P(all.list[[i]], all.list[[j]], cor.mat[i,j])
      cor.mat.star[j,i] = cor.mat.star[i,j]
    }
  }
  if(no.bin+no.ord > 1) {
    start = no.pois+1
    end = no.bin+no.ord+no.pois
    cor.mat.star[start:end,start:end] = ordcont(all.list[start:end], cor.mat[start:end,start:end])$SigmaC
  }
  if(no.nonn > 1) {
    nonn.combo = combn(no.nonn, 2)+no.pois+no.bin+no.ord
    for(k in 1:ncol(nonn.combo)) {
      i=nonn.combo[1,k]
      j=nonn.combo[2,k]
      cor.mat.star[i,j] = .cor.NN.NN(all.list[[i]], all.list[[j]], cor.mat[i,j])
      cor.mat.star[j,i] = cor.mat.star[i,j]
    }
  }

  if(no.pois > 0 && no.bin+no.ord > 0) {
    po.combo = expand.grid(1:no.pois, 1:(no.bin+no.ord)+no.pois)
    for(k in 1:nrow(po.combo)) {
      i=po.combo[k,1]
      j=po.combo[k,2]
      temp = .cor.P.N(all.list[[i]]) * .cor.O.N(all.list[[j]])
      cor.mat.star[i,j] = cor.mat[i,j] / temp
      cor.mat.star[j,i] = cor.mat.star[i,j]
    }
  }
  if(no.pois > 0 && no.nonn > 0) {
    pnn.combo = expand.grid(1:no.pois, 1:no.nonn+no.bin+no.ord+no.pois)
    for(k in 1:nrow(pnn.combo)) {
      i=pnn.combo[k,1]
      j=pnn.combo[k,2]
      temp = .cor.P.N(all.list[[i]]) * .cor.NN.N(all.list[[j]])
      cor.mat.star[i,j] = cor.mat[i,j]/temp
      cor.mat.star[j,i] = cor.mat.star[i,j]
    }
  }
  if(no.nonn > 0 && no.bin+no.ord > 0) {
    onn.combo = expand.grid(1:(no.bin+no.ord)+no.pois,1:no.nonn+no.bin+no.ord+no.pois)
    for(k in 1:nrow(onn.combo)) {
      i=onn.combo[k,1]
      j=onn.combo[k,2]
      temp = .cor.O.N(all.list[[i]]) * .cor.NN.N(all.list[[j]])
      cor.mat.star[i,j] = cor.mat[i,j] / temp
      cor.mat.star[j,i] = cor.mat.star[i,j]
    }
  }
  if(any(is.na(cor.mat.star))) {
    stop("Oops! Something went wrong!")
  }
  if(!is.positive.definite(cor.mat.star)) {
    cat("Warning! \nFinal correlation matrix had to be coerced to be positive definite!")
    cor.mat.star = nearPD(cor.mat.star, corr=TRUE)$mat
  }
  return(cor.mat.star)
}

#Follows the method of Demirtas and Hedeker, 2011.
.getMaxMinCorr <- function(gen1,arg1,gen2,arg2,n=100000) {
  x=gen1(n,arg1)
  y=gen2(n,arg2)
  xs=sort(x)
  ys=sort(y)
  yd=sort(y,decreasing=TRUE)
  maxCor=cor(xs,ys,use="all.obs")
  minCor=cor(xs,yd,use="all.obs")
  return(c(maxCor,minCor))
}

.checkMoments <- function(mvector) {
  return(mvector[2] > 0 && (mvector[4] >= mvector[3]^2-2))
}

.fleishman <- function(x, r1, r2) {
  b <- x[1]
  c <- x[2]
  d <- x[3]
  f <- rep(NA, 3)
  f[1] <- b^2 + 6 * b * d + 2 * c^2 + 15 * d^2 - 1
  f[2] <- 2*c * (b^2 + 24*b*d + 105*d^2 + 2) - r1
  f[3] <- b*d + c^2 * (1 + b^2 + 28 * b * d) + d^2 * (12 + 48 * b* d + 141 * c^2 + 225 * d^2) - r2/24
  f
}

.BBsolveWrapper <- function(arg1,arg2) {
  ans=BBsolve(par=runif(3),fn=.fleishman,r1=arg1,r2=arg2,quiet=TRUE)
  if(arg1*ans$par[2] < 0) ans$par = -1*ans$par
  return(ans)
}

#Solves system of equations to find (b,c,d) using BB package and tacks on mean and variance
.getBCD <- function(rs, maxIters=100) {
  rmean = rs[1]
  rvar = rs[2]
  r1 = rs[3]
  r2 = rs[4]
  ans=.BBsolveWrapper(r1, r2)
    iters = 0
    while(iters < maxIters && (ans$convergence != 0 || !.testC4(ans$par))) {
      iters = iters + 1
      ans = .BBsolveWrapper(r1,r2)
    }
  return(c(ans$par,rmean,rvar))
}

.testC4 <- function(bcd1) {
  return(bcd1[3] > sqrt((5+7*bcd1[1]^2)/75)-2*bcd1[1]/5)
}

.n.to.nonn <- function(z, bcd1) {
  bcd4 = c(-bcd1[2],bcd1[1:3])
  y = bcd4 %*% rbind(rep(1,length(z)),z,z^2,z^3)
  return(bcd1[4]+sqrt(bcd1[5])*as.vector(y))
}

.rnonn <- function(n, bcd1) {
  z = rnorm(n)
  bcd4 = c(-bcd1[2],bcd1[1:3])
  y = bcd4 %*% rbind(rep(1,n),z,z^2,z^3)
  return(bcd1[4]+sqrt(bcd1[5])*as.vector(y))
}

.getMaxMin <- function(x, no.pois, no.bin, no.ord, no.nonn, all.list) {
  n1 = x[1]
  n2 = x[2]
  arg1 = all.list[[n1]]
  arg2 = all.list[[n2]]
  gen1 = .rnonn
  gen2 = .rnonn
  if(n1 <= no.pois) gen1 = rpois
  else if(n1 <= no.pois+no.bin+no.ord) gen1 = .rord
  if(n2 <= no.pois) gen2 = rpois
  else if(n2 <= no.pois+no.bin+no.ord) gen2 = .rord
  return(.getMaxMinCorr(gen1,arg1,gen2,arg2))
}

.n.to.ord <- function(z, cum.vec) {
  u = pnorm(z)
  return(sapply(u, .numless, cum.vec))
}

.rord <- function(n, cum.vec) {
  u = runif(n)
  return(sapply(u, .numless, cum.vec))
}

.numless <- function(u, cum.vec) {
  return(sum(cum.vec<=u))
}

.cor.P.P <- function(lambda1, lambda2, cor12, n=10000) {
  u = runif(n)
  maxcor = cor(qpois(u, lambda1), qpois(u, lambda2))
  mincor = cor(qpois(u, lambda1), qpois(1 - u, lambda2))
  a = -maxcor * mincor/(maxcor + mincor)
  b = log((maxcor + a)/a, exp(1))
  c = -a
  corrected = log((cor12 + a)/a, exp(1))/b
  corrected = ifelse((corrected > 1 | corrected < (-1)), 
                     NA, corrected)
  return(corrected)
}

.cor.P.N <- function(lambda, n=10000) {
  y = rnorm(n)
  x = pnorm(y)
  p = qpois(x, lambda)
  return(cor(p,y))
}

.cor.O.N <- function(cum.vec, n=10000) {
  x = rnorm(n)
  u = pnorm(x)
  o = sapply(u, .numless, cum.vec)-1
  return(cor(x,o))
}

.cor.NN.NN <- function(bcd1, bcd2, rho) {
  coeff=c(-rho,bcd1[1]*bcd2[1]+3*bcd1[1]*bcd2[3]+3*bcd1[3]*bcd2[1]+9*bcd1[3]*bcd2[3],2*bcd1[2]*bcd2[2],6*bcd1[3]*bcd2[3])
  y=Re(polyroot(coeff))
  tf = y<=1 & y>=-1
  if(sum(tf) > 1) {
    b = coeff%*%rbind(rep(1,3),y,y^2,y^3)
    tf = tf & abs(b) < 10^-6
  }
  if(sum(tf) > 1) {
    tf = tf & abs(y-rho)==min(abs(y-rho)) & !duplicated(y)
  }
  rhoN = y[tf]
  return(rhoN)
}

.cor.NN.N <- function(bcd) {
  return(1/(bcd[1]+3*bcd[3]))
}

.getCum <- function(p.vec) {
  cum.vec = vector(length(p.vec)-1,mode="numeric")
  for (k in 1:length(cum.vec) )
    cum.vec[k] = sum(p.vec[1:k])
  return(cum.vec)
}

.n.to.pois <- function(z, lambda) {
  u = pnorm(z)
  return(qpois(u, lambda))
}

.make.sym.matrix <-function(numVars, rhos) {
  sigma = diag(numVars)
  sigma[upper.tri(sigma)]=rhos
  sigma[lower.tri(sigma)]=rhos
  return(sigma)
}