"dirichlet" <-
  function(
                      cat.pen, # Category Penetration (catpen)
                      cat.buyrate, # Category Buyer's Average Purchase Rate (catbuyrate) in a given period.
                      brand.share, # Brand's Market Share
                      brand.pen.obs, # Brand Penetration

                      brand.name = NA,

                      #######  Other Tuning Parameters ##################

                      cat.pur.var = NA, # optionally use the VARIANCE of purchase rate to estimate K
                      nstar = 50,      # Max Purchase Rate to Approximate Infinite Prob. Sum
                      max.S = 30,       # maximum S 
                      max.K = 30,       # Maximum K
                      check = F         # Whether print diagnostic lines
                      )
{
  
  nbrand <- length(brand.pen.obs)             # number of brands considered
  if (is.na(brand.name)[1]) { brand.name <- paste("B",1:nbrand,sep="")}
  
  M0 <- M <- cat.pen * cat.buyrate  # NBD Parameter (base Period, Period=1)

  ## Estimate K
  if (is.na(cat.pur.var)) {
    cp <- log(1-cat.pen)
    eq1 <- function(K) { (K*log(1+M/K) + cp)^2 }  # fit by mean and zeros (if dist is reverse J)
    
    r <- optimize(eq1, c(0.0001,max.K))
    K <- r$minimum # NBD Parameter
  } else {
    K <- M^2 / (cat.pur.var - M)        # cat.pur.var > M , always! Wrong if otherwise.
  }

###########  Probability Functions for Estimating Dirichlet Model ##################
  
  ## For consumers buying the category n times in the analysis period, the conditional prob. of NOT
  ## buying brand j is pzeron(n,j,S), where S is the unknown parameter.
  ## This function is only used to find optimum S. Usually use "p.rj.n(rj,n,j)" below.
  pzeron <-  function(n, j, S) {
    alphaj <- S * brand.share[j]
    if (n==0) { r <- 1}
    else {
      a <- 0:(n-1)
      num <- log(S-alphaj+a)
      den <- log(S+a)
      r <- exp(sum(num-den))
    }
    r
   ## equiv to (but 10 times SLOW): r=1; for(i in 0:(n-1)) r=r * (S - alphaj + i) / (S + i)  
  }

  ## For consumers buying the category n times in the analysis period, the conditional prob. of
  ## buying brand j for exactly rj times.
  ## Note p.rj.n(0, n, j) == pzeron(n,j,S)
  ## last argument "j" can be a vector so that we can combine two brands j and k together
  ## => alpha=alpha_j + alpha_k
  p.rj.n <-  function(rj, n, j) {
    alphaj <- S * sum(sapply(j, function(x) brand.share[x]))
    choose(n,rj) * beta(alphaj + rj, S - alphaj + n - rj) / beta(alphaj, S-alphaj)
  }

## Prob. Dist. for Theoretical NBD of Category Purchases
Pn <- function(n) {
  if (n==0) { g <- 0 }
  else {
    a <- 0:(n-1)
    num <- log(K+a)
    den <- log(1+a)
    g <- sum(num-den)
  }
  exp((-K)*log(1+M/K) + g + n*log(M/(M+K)))
}

## Brand Penetration (b)
## "j":  Brand j
## "limit": Set of Category Buying Frequency
brand.pen <- function(j,limit=c(0:nstar)) {
   if (check == T) {cat("In brand.pen, nstar=", limit[length(limit)],"\n") }
  ## p(0): Overall Prob. of Not Buying Brand j                                        
  p0 <- sum(sapply(limit, function(i,j) {Pn(i) * p.rj.n(0,i,j)}, j=j))
  1 - p0                               # Brand Penetration
}


## The theoretical number of purchases of brand j per puyer of j (w)
brand.buyrate <- function(j,limit=c(1:nstar)){
  buyrate.n <- function(n,j){         # Given n Category Purchases, expected buying rate
    rate <- 1:n
    sum(rate * sapply(rate, p.rj.n, n=n, j=j))
  }
  if (check == T) {cat("In brand.buyrate, nstar=", limit[length(limit)],"\n") }
  ## w = \sum_{n=1}^\Infinity { Pn \sum_{r=1}^n [ r * p(r|n)] } / [1-p(0)]
  sum(sapply(limit, Pn) * sapply(limit, buyrate.n, j=j)) / brand.pen(j)
}

## and their average number of purchases of the Category (wp)
## (purchase rate of the category for the brand buyers)
wp <- function(j, limit=1:nstar){
  sum( sapply(limit,
              function(n, j) {
                n * Pn(n) * (1 - p.rj.n(0, n, j))
              },
              j=j)
      ) / brand.pen(j)
}

  
###########  END: Probability Functions for Estimating Dirichlet Model ##################3333
  
  ## For Estimating the Dirichlet Parameter S
  eq2 <- function(S,j) {
    ## Theoretical Penetraton
    t.pen <- 1 - sum(sapply(0:nstar, function(i,S,j) {Pn(i) * pzeron(i,j,S)}, S=S, j=j))
    o.pen <- brand.pen.obs[j]       # Observed Prob. of buying Brand j (prop of buyers)
                                        #  cat("S =",S, ", Pen(T) =",t.pen,", Pen (O) =",o.pen,"\n")
    (t.pen - o.pen)^2  # to be minimized to solve the equation.
  }

                                        # For Brand j, the solution of S is
  Sj <- function(j){
    r <- optimize(eq2, c(0, max.S),j=j)
    if (check==T) {
      cat("Objective Value is ", r$objective, "at S=", r$minimum, ", and Brand=",j, "\n")
    }
    r$minimum
  }

  if (check==T) {cat("Finding Optimum S for Each Brand ...\n")}
  Sall <- sapply(1:nbrand, Sj)            # S for each brand
  bp <- boxplot(Sall,plot=F)
  outlier <- bp$out   # delete the outlier of S
  outlier2 <- Sall[Sall>bp$conf[2]]     # also deem numbers outside
                                        # the upper "notch" as "outliers"
  outliers <- c(outlier,outlier2)
  schoose <- sapply(Sall, function(x) {! (x %in% outliers)})
  
  if (check==T && length(outliers)>0) {
    cat("Removing These Outliers of S:", outliers,", for brands:", c(1:nbrand)[!schoose],"\n")
  }    
  ## Final S, weighted by Brand's Market Share
  S <- weighted.mean(Sall[schoose], brand.share[schoose])  
  
########### All Parameters Are Estimated Now! ##################

  ## Check if "nstar" is sufficently large to cover the support of Pn
  ##  Category Prob should add up to 1;   Estimated mean from Pn should be equal to M0
  if (sum(sapply(0:nstar,Pn)) < 0.99 || abs(sum(c(0:nstar)* sapply(0:nstar,Pn))-M0)>0.1) {
    cat("nstar is too small! (nstar=",nstar,")\n")
    error=1
  }
  else {error=0}
  
  ## (Functional) Parameters for Dirichlet Model:
  dpar <- list(S=S, M=M, K=K,   # Estimated 
               nbrand=nbrand,
               nstar=nstar,

               ## Input Parameters
               cat.pen=cat.pen, # Category Penetration (catpen)
               cat.buyrate=cat.buyrate, # Category Buyer's Average Purchase Rate (catbuyrate) in a given period.
               brand.share=brand.share, # Brand's Market Share
               brand.pen.obs=brand.pen.obs, # Brand Penetration
               brand.name=brand.name,

               ## functions to be passed out
               period.set=function(t) {
                 M <<- M0 * t                   # change the time period in M
               },
               period.print=function(){
                 cat("Multiple of Base Time Period:",round(M/M0,2),", Current M =",M,"\n")
               },
               p.rj.n=p.rj.n,   # The conditional prob. of buying brand j for exactly rj times
                                # given n category purchases
               Pn=Pn,           # Prob. Dist. for Theoretical NBD of Category Purchases
               brand.pen=brand.pen,  # Brand Penetration
               brand.buyrate=brand.buyrate,  # Brand Buying Rate
               wp=wp,  # Purchase Rate of the Category for the Brand Buyers

               check=check,
               error=error
               
               )

  class(dpar) <- "dirichlet"

  dpar
}
