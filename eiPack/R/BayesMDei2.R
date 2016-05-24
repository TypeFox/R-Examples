BayesMDei2 <- function(formula, data, total, lambda1 = 4, lambda2 = 2,
                      tune.alpha = NULL, tune.beta = NULL,
                       start.alphas = NULL, start.betas = NULL,
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
  
  NG <- nrow(XX)
  NP <- nrow(TT)
  Precincts <- nrow(DD)

  usrenv <- environment(fun = usrfun)
  
  if(is.null(start.alphas)){
  start.alphas <- matrix(rgamma(NG*NP, lambda1, lambda2), NG, NP)}
  if(min(start.alphas) <= 0){stop("inadmissable starting values for alpha")}
  
  if(is.null(start.betas)){
  start.betas <- array(NA, dim= c(NG, NP, Precincts))
  for(i in 1:Precincts){
    start.betas[,,i] <- rdirichlet(NG, rep(1,NP))}
  }
  if(identical(round(apply(start.betas, c(1,3), sum),10),matrix(1,NG, 
Precincts))!=TRUE){stop("inadmissable 
starting values
for beta")}



  usrlen <- length(as.numeric(usrfun(list(start.alphas, start.betas, TT,
                                          RR))))
  
  if(is.null(tune.alpha)){
    tune.alpha <- matrix(rep(.25,NG*NP), NG, NP)}
  if(is.null(tune.beta)){
    tune.beta <- array(rep(.05, NG*(NP-1)*Precincts), c(NG, NP-1, Precincts))}

  tune.alpha <- as.matrix(tune.alpha)
  
  if(identical(dim(tune.alpha), c(NG, NP))!=TRUE) {stop("'tune.alpha'
has incorrect dimensions")}
  
  if(identical(as.numeric(dim(tune.beta)), c(NG, NP-1, Precincts))!=TRUE) 
{stop("'tune.beta'
has incorrect dimensions")}
  




  beta.names <- paste(paste(paste(group.names,matrix(rep(party.names,
                                                   NG),NG,NP, byrow=T)
                            ,sep="."),
                      matrix(rep(1:Precincts,NG*NP),NG*NP, Precincts,
                             byrow=TRUE),sep="."), ".txt.gz", sep="")
  if(ret.beta == 's'){touch.betas(beta.names)
                      ret.beta <- 2}
  if(ret.beta == 'd'){ret.beta <- 1}
  if(ret.beta == 'r'){ret.beta <- 0}
  if(is.numeric(ret.beta)==FALSE){stop('incorrect option for
ret.beta')}
  

  
  output <- .Call("rbycei_fcn2",
                  as.numeric(start.alphas),
                  as.numeric(start.betas),
                  as.numeric(TT),
                  as.numeric(XX),
                  as.numeric(tune.alpha),
                  as.numeric(tune.beta),
                  as.integer(NG),
                  as.integer(NP),
                  as.integer(Precincts),
                  as.numeric(lambda1),
                  as.numeric(lambda2),
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

  if(ret.beta==0){names(output) <- c("Alpha", "Beta", "alpha.acc",
                                  "beta.acc", "cell.count", "usrfun")}
  else{names(output) <- c("Alpha", "alpha.acc",
                          "beta.acc","cell.count", "usrfun")}
  
  if(ret.mcmc){
  colnames(output$Alpha) <- paste("alpha",matrix(rep(group.names,
                                                     NP),NG,NP)
                                 ,matrix(rep(party.names, NG),NG,NP,
                                         byrow=T) ,sep=".")
  output$Alpha <- coda::mcmc(output$Alpha, thin=thin)
  colnames(output$cell.count) <- paste("ccount",matrix(rep(group.names,
                                                    NP),NG,NP)
                                 ,matrix(rep(party.names, NG),NG,NP,
                                         byrow=T) ,sep=".")
  output$cell.count <- coda::mcmc(output$cell.count, thin=thin)
  colnames(output$usrfun) <- paste("usrfun", 1:ncol(output$usrfun), sep=".")
  output$usrfun <- coda::mcmc(output$usrfun, thin=thin)
  if(ret.beta==0){
    colnames(output$Beta) <- paste(paste("beta", group.names,matrix(rep(party.names, NG),NG,NP, byrow=T) ,sep="."), matrix(rep(1:Precincts,NG*NP),NG*NP, Precincts, byrow=TRUE),sep=".")
    output$Beta <- coda::mcmc(output$Beta, thin=thin)
  }
}else{
  output$Alpha <- array(t(output$Alpha), c(NG, NP, sample))
  dimnames(output$Alpha) <- list(group.names, party.names, 1:sample)
  output$cell.count <- array(t(output$cell.count), c(NG, NP, sample))
  dimnames(output$cell.count) <- list(group.names, party.names,
                                      1:sample)
  colnames(output$usrfun) <- paste("usrfun", 1:ncol(output$usrfun), sep=".")
  if(ret.beta==0){
    output$Beta <- array(t(output$Beta), c(NG, NP, Precincts, sample))
    dimnames(output$Beta) <- list(group.names, party.names, 1:Precincts,
                                1:sample)
  }
}

  return(output)
}




