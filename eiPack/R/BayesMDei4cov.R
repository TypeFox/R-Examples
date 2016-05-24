BayesMDei4cov <- function(formula, covariate, total, data, lambda1 = 2,
                          lambda2 = 4, covariateprior = NULL,
                          tune.dr = NULL, tune.beta = NULL, tune.gamma
                          = NULL, tune.delta = NULL,
                          start.dr = NULL, start.betas = NULL,
                          start.gamma = NULL, start.delta = NULL,
                          sample = 1000,
                          thin = 1, burnin = 1000, verbose = 0, ret.beta =
                          'r', ret.mcmc = TRUE, usrfun = NULL, ...){

  if(thin < 1){stop('thin must be positive integer')}
  if(sample < 1){stop('thin must be positive integer')}
  if(burnin < 0){stop('burnin must be non-negative integer')}
  
  DD <- model.frame(formula, data)


  countParty <- countGroup <- propParty <- propGroup <- FALSE

  checkGroups <- round(apply(DD[[2]], 1, sum), 3)
  checkParties <- round(apply(DD[[1]], 1, sum), 3)
  if(all(DD[[1]] %% 1 == 0) & all(DD[[1]] >= 0)){countParty <- TRUE}
  else if(all(0 <= DD[[1]] & DD[[1]] <= 1)){
        if(all(checkParties == 1)){propParty <- TRUE}else{
        stop("column marginals are proportions that do not
                sum to 1 - please respecify data")}}
  else stop("column marginals are neither counts nor proportions - please
respecify data")

  if(all(DD[[2]] %% 1 == 0) & all(DD[[2]] >= 0)){countGroup <- TRUE}
  else if(all(0 <= DD[[2]] & DD[[2]] <= 1)){
        if(all(checkGroups == 1)){propGroup <- TRUE}else{
        stop("row marginals are proportions that do not sum to 1 - please
respecify data")}}
  else stop("row marginals are neither counts nor proportions - please
respecify data")

  if((propParty | propGroup) & is.null(total)){
    stop("one or both marginals are proportions - 'total' must be
provided")}

  if(propParty & !is.null(total)){
    DD[[1]] <- DD[[1]] * total
  warning("column margnials are proportions - multiplying by unit size")}
  if(propGroup & !is.null(total)){
    DD[[2]] <- DD[[2]] * total
   warning("row margnials are proportions - multiplying by unit size")}

  checkGroups <- round(apply(DD[[2]], 1, sum), 1)
  checkParties <- round(apply(DD[[1]], 1, sum), 1)
  if(identical(checkParties, checkGroups) ==  FALSE){
    stop("row and column totals unequal in some units - please
respecify data")}
  
  Groups <- DD[[2]]
  TT <- t(DD[[1]])
  XX <- t(Groups/apply(Groups,1,sum))
  group.names <- colnames(Groups)
  party.names <- rownames(TT)
  RR <- t(Groups)
  CC <- model.frame(covariate, data)
  ZZ <- as.matrix(CC)
  
  NG <- nrow(XX)
  NP <- nrow(TT)
  Precincts <- nrow(DD)

  if(is.null(start.dr)){
  start.dr <- matrix(rgamma(NG, lambda1, lambda2), NG)}
  if(min(start.dr) <= 0){stop("inadmissable starting values for dr")}

  
  if(is.null(start.betas)){
  start.betas <- array(NA, dim= c(NG, NP, Precincts))
  for(i in 1:Precincts){
    start.betas[,,i] <- rdirichlet(NG, rep(1,NP))}
  }
  if(identical(round(apply(start.betas, c(1,3), sum),10),matrix(1,NG, 
Precincts))!=TRUE){stop("inadmissable 
starting values
for beta")}

  
  if(is.null(start.gamma)){
  start.gamma <- cbind(matrix(rnorm(NG*(NP-1)), NG, NP-1),0)}
  if(identical(start.gamma[,NP], rep(0,NG))!=TRUE){stop("final column
of 'start.gamma' must be zero")}
  
  if(is.null(start.delta)){
  start.delta <- cbind(matrix(rnorm(NG*(NP-1)), NG, NP-1),0)}
  if(identical(start.delta[,NP], rep(0,NG))!=TRUE){stop("final column
of 'start.delta' must be zero")}
  
  
  usrenv <- environment(fun = usrfun)
  usrlen <- length(as.numeric(usrfun(list(start.dr, start.betas,
                                          start.gamma, start.delta, TT,
                                          RR))))
  
  if(is.null(tune.dr)){
    tune.dr <- rep(2,NG)}
  if(is.null(tune.beta)){
    tune.beta <- array(rep(.05, NG*(NP-1)*Precincts), c(NG, NP-1, Precincts))}
  if(is.null(tune.gamma)){
    tune.gamma <- matrix(.25, NG, NP-1)}
  if(is.null(tune.delta)){
    tune.delta <- matrix(.25, NG, NP-1)}
  
  if(identical(length(tune.dr), NG)!=TRUE) {stop("'tune.dr'
has incorrect dimensions")}
  
  if(identical(as.numeric(dim(tune.beta)), c(NG, NP-1, Precincts))!=TRUE) 
{stop("'tune.beta'
has incorrect dimensions")}

  if(identical(as.numeric(dim(tune.gamma)), c(NG, NP-1))!=TRUE) 
{stop("'tune.gamma'
has incorrect dimensions")}

  if(identical(as.numeric(dim(tune.delta)), c(NG, NP-1))!=TRUE) 
{stop("'tune.delta'
has incorrect dimensions")}

  if(is.null(covariateprior)){
    covprior <- 0
    delmean <- gammean <- rep(0, NG*(NP-1))
    delsd <- gamsd <- rep(1, NG*(NP-1))
  }else{
    covprior <- 1
    delmean <- covariateprior[[1]]
    delsd <- covariateprior[[2]]
    gammean <- covariateprior[[3]]
    gamsd <- covariateprior[[4]]
    if(identical(as.numeric(dim(delmean)), c(NG, NP-1))!=TRUE) 
      {stop("matrix of prior means for delta has incorrect dimensions")}
    if(identical(as.numeric(dim(delsd)), c(NG, NP-1))!=TRUE) 
      {stop("matrix of prior sd for delta has incorrect dimensions")}
    if(identical(as.numeric(dim(gammean)), c(NG, NP-1))!=TRUE) 
      {stop("matrix of prior means for gamma has incorrect dimensions")}
    if(identical(as.numeric(dim(gamsd)), c(NG, NP-1))!=TRUE) 
      {stop("matrix of prior sd for gamma has incorrect dimensions")}
    if(min(gamsd)<=0)
      {stop("prior sd for gamma must be > 0")}
    if(min(gamsd)<=0)
      {stop("prior sd for delta must be > 0")}
  }
  



  beta.names <- paste(paste(paste(group.names,matrix(rep(party.names,
                                                   NG),NG,NP, byrow=T)
                            ,sep="."),
                      matrix(rep(1:Precincts,NG*NP),NG*NP, Precincts,
                             byrow=TRUE),sep="."), ".txt.gz", sep="")

  if(ret.beta == 's'){touch.betas(beta.names)
                      ret.beta <- 2}
  if(ret.beta == 'd'){ret.beta <- 1}
  if(ret.beta == 'r'){ret.beta <- 0}
  if(is.numeric(ret.beta)==FALSE){stop("incorrect option for
ret.beta")}
  

  
  output <- .Call("rbycei_fcn4",
                  as.numeric(start.dr),
                  as.numeric(start.betas),
                  as.numeric(start.gamma),
                  as.numeric(start.delta),
                  as.numeric(TT),
                  as.numeric(XX),
                  as.numeric(ZZ),
                  as.numeric(tune.dr),
                  as.numeric(tune.beta),
                  as.numeric(tune.gamma),
                  as.numeric(tune.delta),
                  as.integer(NG),
                  as.integer(NP),
                  as.integer(Precincts),
                  as.numeric(lambda1),
                  as.numeric(lambda2),
                  as.integer(covprior),
                  as.numeric(delmean),
                  as.numeric(delsd),
                  as.numeric(gammean),
                  as.numeric(gamsd),
                  as.integer(sample),
                  as.integer(thin),
                  as.integer(burnin),
                  as.integer(verbose),
                  as.integer(ret.beta),
                  as.numeric(RR),
                  usrfun,
                  usrenv,
                  as.integer(usrlen),
                  as.character(beta.names)
                  )

  if(ret.beta==0){names(output) <- c("Dr", "Beta","Gamma","Delta",
                                  "dr.acc","beta.acc", "gamma.acc",
                                  "delta.acc","cell.count", "usrfun")}
  else{names(output) <- c("Dr","Gamma","Delta",
                                  "dr.acc","beta.acc", "gamma.acc",
                                  "delta.acc","cell.count", "usrfun")}
  
  if(ret.mcmc){
  colnames(output$Dr) <- paste("dr", group.names, sep=".")
  output$Dr <- coda::mcmc(output$Dr, thin=thin)
  colnames(output$cell.count) <- paste("ccount",matrix(rep(group.names,
                                                    NP),NG,NP)
                                 ,matrix(rep(party.names, NG),NG,NP,
                                         byrow=T) ,sep=".")
  output$cell.count <- coda::mcmc(output$cell.count, thin=thin)
  colnames(output$Gamma) <- paste("gamma",matrix(rep(group.names,
                                                    (NP-1)),NG,NP-1)
                                 ,matrix(rep(party.names[1:(NP-1)], NG),NG,NP-1,
                                         byrow=T) ,sep=".")
  output$Gamma <- coda::mcmc(output$Gamma, thin=thin)
  colnames(output$Delta) <- paste("delta",matrix(rep(group.names,
                                                    (NP-1)),NG,NP-1)
                                 ,matrix(rep(party.names[1:(NP-1)], NG),NG,NP-1,
                                         byrow=T) ,sep=".")
  output$Delta <- coda::mcmc(output$Delta, thin=thin)

  if(ret.beta==0){
    colnames(output$Beta) <- paste(paste("beta", group.names,matrix(rep(party.names, NG),NG,NP, byrow=T) ,sep="."), matrix(rep(1:Precincts,NG*NP),NG*NP, Precincts, byrow=TRUE),sep=".")
    output$Beta <- coda::mcmc(output$Beta, thin=thin)
  }
}else{

  output$Dr <- t(output$Dr) 
  dimnames(output$Dr) <- list(paste("dr", group.names, sep="."), 1:sample)
  output$cell.count <- array(t(output$cell.count), c(NG, NP, sample))
  dimnames(output$cell.count) <- list(group.names, party.names,
                                      1:sample)
  output$Gamma <- array(t(output$Gamma), c(NG, NP-1, sample))
  dimnames(output$Gamma) <- list(group.names, party.names[1:(NP-1)], 1:sample)
  output$Delta <- array(t(output$Delta), c(NG, NP-1,sample))
  dimnames(output$Delta) <- list(group.names, party.names[1:(NP-1)], 1:sample)
  
  if(ret.beta==0){
    output$Beta <- array(t(output$Beta), c(NG, NP, Precincts, sample))
    dimnames(output$Beta) <- list(group.names, party.names, 1:Precincts,
                                1:sample)
  }
}
    
  return(output)
}




