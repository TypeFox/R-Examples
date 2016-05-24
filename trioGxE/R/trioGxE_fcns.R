## functions used in trioGxE() function:
## pre.trioGxE(), trio.mustart(), link(), linkinv(), deriv.fcn(), triogam.loglkhd()
pre.trioGxE <- function(data,pgenos,cgeno,cenv,testGxE=FALSE) {
  #manipulate the data so that it can be in the appropriate form 
  # after cleaning up the data
  dat = data
  x <- dat[,cenv]
  pgenos <- dat[,pgenos]
  gm = pgenos[,1]; gf = pgenos[,2]
  cgeno <- dat[,cgeno]
  
  if(!testGxE){
    # data-checking
    # when testing for GxE, data-checking is not done
    # how many trios have complete genotypes and non-genetic covariates?
    available <- length(x)
    ind <- (!is.na(gm) & !is.na(gf)) & (!is.na(cgeno) & !is.na(cenv))
    complete <- sum(ind)
    ## Include only those trios with complete geno data
    if(sum(!ind)){
      gm <- gm[ind,]; gf <- gf[ind,]; cgeno <- cgeno[ind,]; x <- x[ind]
    }
    
    ## Genotype error checking
    if(!all(gm %in% c(0,1,2)) | !all(gf %in% c(0,1,2)))
      stop("Parental genenotypes ('pgenos') must be 0, 1, or 2.")
    if(!(all(cgeno %in% c(0,1,2))))
      stop("Child genotypes ('cgeno') must be 0, 1, or 2.")
    ## Checking whether too values exist for the non-genetic covariate 
    if(length(unique(x))<4 & (!is.factor(x)))
      stop("Nongenetic covariate ",cenv, "has fewer than 4 unique values.")
    
    ## Among trios with complete geno data, how many are informative?
    ind <- (gm==1 | gf==1)
    informative <- sum(ind)
    if(informative<20)
      warning(paste("Only ",informative, " trios with informative parental genotypes"))
    ## Include only informative trios
    if(sum(!ind)){
      gm <- gm[ind,]; gf <- gf[ind,]; cgeno <- cgeno[ind,]; x <- x[ind]
    }
    
    Gp <-rep(NA,length(x)) # represents parental mating types
    ind1 <- (gm+gf) == 1 #Gp=1: heterozygote + homozygote for non-index allele
    ind2 <- (gm+gf) == 3 #Gp=2: heterozygote + homozygote for index allele
    ind3 <- (gm+gf) == 2 #Gp=3: two heterozygotes
    Gp[ind1] <- 1; Gp[ind2] <- 2; Gp[ind3] <- 3
    
    ## Mendelian inconsistencies
    ind <- rep(FALSE,length(x))
    ind[Gp==1 & cgeno==2] <- TRUE #children with G=2 can't be from Gp=1
    ind[Gp==2 & cgeno==0] <- TRUE #children with G=0 can't be from Gp=2
    mendelian.consistent <- length(x) - sum(ind)
    
    ## include only those that are Mendelian-consistent
    if(sum(ind)){
      Gp <- Gp[!ind]; cgeno <- cgeno[!ind]; x <- x[!ind] 
      ind1 <- ind1[!ind]; ind2 <- ind2[!ind]; ind3 <- ind3[!ind]
    }
  }
  
  else{
    Gp <-rep(NA,length(x)) # represents parental mating types 
    ind1 <- (gm+gf) == 1 #Gp=1: heterozygote + homozygote for non-index allele
    ind2 <- (gm+gf) == 3 #Gp=2: heterozygote + homozygote for index allele
    ind3 <- (gm+gf) == 2 #Gp=3: two heterozygotes
    Gp[ind1] <- 1; Gp[ind2] <- 2; Gp[ind3] <- 3
  }
  
  ## Trios with Gp=3: response variable is a vector - need manipulation
  n<-length(x)
  id <- c(1:n)
  if(any(ind3)){
    id <- c(id,id[ind3])  
    mt = Gp
    mt[ind3] <- 3.1
    mt <- c(mt, rep(3.2, sum(ind3)))
    cgeno <- c(cgeno, cgeno[ind3])
    x <- c(x, x[ind3])
  }

  # MATING-TYPE-SPECIFIC RESPONSE
  y <- rep(NA, length(id))
  #mt1
  y[cgeno==0 & mt==1 ] = 0
  y[cgeno==1 & mt==1 ] = 1 
  #mt2
  y[cgeno==1 & mt==2 ] = 0
  y[cgeno==2 & mt==2 ] = 1
  #mt3
  #G=++
  y[cgeno==0 & mt==3.1 ] = 0 #y1
  y[cgeno==0 & mt==3.2 ] = 0 #y2
  
  #G=+-
  y[cgeno==1 & mt==3.1 ] = 1 #y1
  y[cgeno==1 & mt==3.2 ] = 0 #y2

  #G=--
  y[cgeno==2 & mt==3.1 ] = 0 #y1
  y[cgeno==2 & mt==3.2 ] = 1 #y2

  res <- data.frame(id=id, x=x, y=y, mt=mt)
  res <- res[order(id),]  
  attr(res,"Gp") = Gp
  res
}

## for fitting
trio.mustart<- function(y, mt){
  n <- length(y)
  
  ind= (mt==1|mt==2)
  gcats <- rep(NA,n)
  gcats[ind] <- 2

  ind= (floor(mt)==3)
  gcats[ind] <- 3
  
  trio.mustart <- (1+y)/(1+gcats)
  trio.mustart
}

link <- function(mu,mt){
  #link function for trios from MT1 or MT2 (i.e., eta)
  #NEEDED TO BE CHECKED
  n <- length(mu)
  link.val <- rep(NA,n)

  ind12 <- mt==1 | mt==2
  link.val[ind12] <- log(mu[ind12]/(1-mu[ind12]))

  #MT 3
  ind31 <- mt==3.1; ind32 <- mt==3.2
  link.val[ind31] <- log(mu[ind31]/(1-(mu[ind31]+mu[ind32])))
  link.val[ind32] <- log(mu[ind32]/mu[ind31]) ###
  link.val
}

linkinv <- function(eta, mt)
{##calculating g^(-1)(eta) = g^(-1)(g(mu))=mu
  
  n <- length(eta)
  linkinv.val <- rep(NA, n)
    
  ind12 <- mt==1 | mt==2
  ind31 <- mt==3.1; ind32 <- mt==3.2

  linkinv.val[ind12] <- exp(eta[ind12])/(1+exp(eta[ind12]))
  linkinv.val[ind31] <- exp(eta[ind31])/(1+exp(eta[ind31])+exp(eta[ind31]+eta[ind32]))
  linkinv.val[ind32] <- exp(eta[ind31]+eta[ind32])/(1+exp(eta[ind31])+exp(eta[ind31]+eta[ind32]))

  linkinv.val#vector
}

deriv.fcn <- function(y, mu, mt, w=1){
  #y and mu are assumed to have the same length
  n <- length(y)
  deriv.val <- rep(NA, n)
  ind12 <- mt==1 | mt==2
  deriv.val[ind12] <- y[ind12]-mu[ind12]
  ind31 <- mt==3.1; ind32 <- mt==3.2
  deriv.val[ind31]=(y[ind31]+y[ind32])-(mu[ind31]+mu[ind32])
  deriv.val[ind32]=y[ind32]-mu[ind32]

  deriv.val <- deriv.val*w

  deriv.val
}

triogam.loglkhd <- function(mt,y,eta,w=1){
  ind12 <- mt==1|mt==2; ind31 <- mt==3.1; ind32 <- mt==3.2

  ##n: number of trios + number of trios from MT3
  l <- rep(0,length(mt))
  if(any(ind12))
    l[ind12] <- y[ind12]*eta[ind12]-log(1+exp(eta[ind12]))

  ##leaving entries for mt==3.2 as zeor so that they won't contribute to
  ##calculation of log-likelihood
  if(any(ind31) | any(ind32)){
    if(any(ind31) & any(ind32)){
      l[ind31] <- y[ind31]*eta[ind31]+y[ind32]*(eta[ind31]+eta[ind32])-
        log(1+exp(eta[ind31])+exp(eta[ind31]+eta[ind32]))
    }
    
    else if(any(ind31) & !any(ind32)){
      l[ind31] <- y[ind31]*eta[ind31]+y[ind32]*(eta[ind31]-log(2))-
        log(1+exp(eta[ind31])+exp(eta[ind31]-log(2)))
    }
    
    else{ #(!any(ind31) & any(ind32))
      l[ind32] <- y[ind31]*log(2)+y[ind32]*(log(2)+eta[ind32])-
        log(1+exp(log(2))+exp(log(2)+eta[ind32]))
    }
  }#if(any(ind31) | any(ind32))
  
  l <- w*sum(l)
  l
}