###########################################################################
## sample from the posterior distribution of a hierarchical beta binomial
## model in R using linked C++ code in Scythe
##
## y_{ij} ~ Binomial(s_{ij}, theta_{ij})
## theta_{ij} ~ Beta(alpha_j, beta_j)
## alpha_j ~ Pareto(1, a)
## beta_j  ~ Pareto(1, b)
##
## KQ 5/24/2011 - 6/25/2011
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################


"MCMChierBetaBinom" <- function(y, s, i.labels, j.labels, burnin=1000,
                                mcmc=10000, thin=1,
                                verbose=0, seed=NA, theta.start=NA,
                                alpha.start=NA, beta.start=NA,
                                a=0, b=0){

  ## checks
  check.mcmc.parameters(burnin, mcmc, thin)
  max.length <- max(length(y), length(s), length(i.labels), length(j.labels))
  min.length <- min(length(y), length(s), length(i.labels), length(j.labels))
  if (max.length != min.length){
    stop("y, s, i.labels, j.labels not all of same length\n")
  }

  na.indic <- is.na(y)
  na.indic <- na.indic + is.na(s)
  na.indic <- na.indic + is.na(i.labels)
  na.indic <- na.indic + is.na(j.labels)
  na.indic <- na.indic != 0
  y <- y[!na.indic]
  s <- s[!na.indic]
  i.labels <- i.labels[!na.indic]
  j.labels <- j.labels[!na.indic]

  if (min(y) < 0){
    stop("y cannot have an element less than 0\n")
  }
  if (max(y > s) > 0){
    stop("y[i] cannot be greater than s[i] for any i\n")
  }

  ## more checking needed


  ## seeds
  seeds <- form.seeds(seed) 
  lecuyer <- seeds[[1]]
  seed.array <- seeds[[2]]
  lecuyer.stream <- seeds[[3]]
  
  

  i.labels <- as.character(i.labels)
  j.labels <- as.character(j.labels)
  i.labels.unique <- unique(i.labels)
  j.labels.unique <- unique(j.labels)
  i.labels.num <- rep(NA, length(i.labels))
  j.labels.num <- rep(NA, length(j.labels))
  for (i in 1:length(i.labels)){
    i.labels.num[i] <- which(i.labels.unique == i.labels[i])
  }
  for (j in 1:length(j.labels)){
    j.labels.num[j] <- which(j.labels.unique == j.labels[j])
  }
  i.labels.num.unique <- unique(i.labels.num)
  j.labels.num.unique <- unique(j.labels.num)
                      


  ## starting values
  if (length(theta.start) == 1){
    if (is.numeric(theta.start)){
      theta.start <- rep(theta.start, length(y))
    }
    else if (is.na(theta.start)){
      theta.start <- (y + 0.01) / (s + 0.02)
    }
    else{
      theta.start <- rep(0.5, length(y))
    }
  }
  if (length(theta.start) != length(y)){
    stop("length(theta.start) != length(y)\n") 
  }
  if (max(theta.start) >= 1.0){
    stop("elements of theta.start must be less than 1.0\n")
  }
  if (min(theta.start) <= 0.0){
    stop("elements of theta.start must be greater than 0.0\n")
  }
  
  if (length(alpha.start) == 1){
    if (is.na(alpha.start)){
      alpha.start <- rep(1.001, length(j.labels.unique))
    }    
  }
  if (length(alpha.start) != length(j.labels.unique)){
    stop("length(alpha.start) != length(unique(j.labels))\n") 
  }
  if (length(beta.start) == 1){
    if (is.na(beta.start)){
      beta.start <- rep(1.001, length(j.labels.unique))
    }    
  }
  if (length(beta.start) != length(j.labels.unique)){
    stop("length(beta.start) != length(unique(j.labels))\n") 
  }


  accepts <- rep(0, length(j.labels.unique))



  

  ## get reasonable values for base.sigma
  base.sigma <- rep(1, length(j.labels.unique))

   

  ## call C++ code to draw the sample
  sample <- matrix(data=0, mcmc/thin, (length(y) +
                   2*length(j.labels.unique)) )
  
  
  posterior <- .C("hierBetaBinom",
                  samdata = as.double(sample),
                  samrow = as.integer(nrow(sample)),
                  samcol = as.integer(ncol(sample)),
                  y = as.integer(y),
                  s = as.integer(s),
                  theta = as.double(theta.start),
                  alpha = as.double(alpha.start),
                  beta = as.double(beta.start),
                  a = as.double(a),
                  b = as.double(b),
                  ilabels = as.integer(i.labels.num),
                  jlabels = as.integer(j.labels.num),
                  ilabelsunique = as.integer(i.labels.num.unique),
                  jlabelsunique = as.integer(j.labels.num.unique),
                  n = as.integer(length(y)),
                  ni = as.integer(length(i.labels.unique)),
                  nj = as.integer(length(j.labels.unique)),
                  burnin = as.integer(burnin),
                  mcmc = as.integer(mcmc),
                  thin = as.integer(thin),
                  lecuyer = as.integer(lecuyer),
                  seedarray = as.integer(seed.array),
                  lecuyerstream = as.integer(lecuyer.stream),
                  verbose = as.integer(verbose),
                  accepts = as.integer(accepts),
                  basesigma = as.double(base.sigma),
                  PACKAGE="MCMCpack")

  
  sample <- matrix(posterior$samdata,
                   posterior$samrow,
                   posterior$samcol,
                   byrow=FALSE)

  output <- mcmc(data=sample, start=burnin+1, end=burnin+mcmc, thin=thin)

  theta.names <- paste("theta.", i.labels, ".", j.labels, sep="")
  alpha.names <- paste("alpha.", j.labels.unique, sep="")
  beta.names <- paste("beta.", j.labels.unique, sep="")

  varnames(output) <- c(theta.names, alpha.names, beta.names)

  attr(output, "title") <- "MCMChierBetaBinom Posterior Sample"
  attr(output, "acceptance.rates") <- posterior$accepts / (posterior$mcmc + posterior$burnin)
  return(output)

  
} ## end MCMChierBetaBinom























                                                        










































