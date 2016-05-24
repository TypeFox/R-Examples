##
## Functions for likelihood based inference in the geoR package
## ------------------------------------------------------------
##
## Includes:
##    - the main function for ML estimation - likfit()
##    - methods for class likGRF
##    - some other functions

".check.geoRparameters.values" <- function(list, messages = TRUE)
{
  if(!is.null(list$nugget))
    if(list$nugget < 0) stop("value for nugget must be non-negative")
  if(!is.null(list$psiA))
    if(list$psiA < 0) stop("value for psiA must be non-negative")
  if(!is.null(list$psiR))
    if(list$psiR < 1) stop("value for psiA must be >= 1")
  if(!is.null(list$kappa)){
    if(is.null(list$cov.model)) stop("cov.model must be provided when checking values of kappa")
    if(list$kappa[1] <= 0) stop("parameter kappa[1] must be greater than 0")
    if(any(list$cov.model == c("powered.exponential", "matern", "gneiting.matern", "gencauchy","cauchy"))){
      if(any(list$cov.model == c("gneiting.matern", "gencauchy"))){
        if(length(list$kappa) != 2)
          stop(paste("kappa must be of length 2 for the", list$cov.model, "correlation function"))
        if(list$cov.model == "gencauchy" && (list$kappa[2] <=0 | list$kappa[2] >2))
          stop("for the gencauchy model the kappa must be within (0,2]")          
      }
      else{
        if(list$cov.model == "powered.exponential" && (list$kappa <=0 | list$kappa >2))
          stop("for the powered.exponential model the kappa must be within (0,2]")          
        if(list$cov.model == "stable" && (list$kappa <=0 | list$kappa >2))
          stop("for the stable model the kappa must be within (0,2]")          
      }
    }
    else{
      if(messages && !is.null(list$kappa))
        cat(paste("kappa not used for the",list$cov.model, "correlation function\n"))
    }
  }
  return(invisible())
}

"likfit" <-
  function (geodata, coords=geodata$coords, data=geodata$data,
            trend = "cte", ini.cov.pars,
            fix.nugget = FALSE, nugget = 0, 
            fix.kappa = TRUE, kappa = 0.5, 
            fix.lambda = TRUE, lambda = 1, 
            fix.psiA = TRUE, psiA = 0, fix.psiR = TRUE, psiR = 1, 
            cov.model, realisations, lik.method = "ML",
            components = TRUE, nospatial = TRUE, limits = pars.limits(), 
            print.pars = FALSE, messages, ...) 
{
  name.geodata <- deparse(substitute(geodata))
  ##
  ## Checking input
  ##
  call.fc <- match.call()
  ldots <- list(...)
  temp.list <- list()
  temp.list$print.pars <- print.pars
  if(missing(messages))
    messages.screen <- as.logical(ifelse(is.null(getOption("geoR.messages")), TRUE, getOption("geoR.messages")))
  else messages.screen <- messages
  ##
  if(!missing(ini.cov.pars)){
    if(any(class(ini.cov.pars) == "eyefit")){
ini.cov.pars <- ini.cov.pars[[1]]
#      cov.model <- ini.cov.pars[[1]]$cov.model
#      kappa <- ini.cov.pars[[1]]$kappa
    }
    if(any(class(ini.cov.pars) == "variomodel")){
      cov.model <- ini.cov.pars$cov.model
      kappa <- ini.cov.pars$kappa
}
  }
  if(missing(cov.model)) cov.model <- "matern"
  cov.model <- match.arg(cov.model, choices = .geoR.cov.models)
  if(cov.model == "stable") cov.model <- "powered.exponential"
  if(any(cov.model == c("power", "gneiting.matern", "gencauchy")))
     stop(paste("parameter estimation for", cov.model, "is not yet implemented"))
  ##  if(any(cov.model == c("gneiting.matern", "gencauchy"))){
  ##    if(length(kappa != 2))
  ##      stop(paste(cov.model, "requires two values in the argument kappa"))
  ##    if(length(fix.kappa) == 1) fix.kappa <- rep(fix.kappa, 2) 
  ##    stop("parameter estimation for gneiting.matern model is not yet implemented")
  ##  }
  fixed.pars <- list(cov.model=cov.model)
  if(fix.nugget) fixed.pars$nugget <- nugget
  if(fix.kappa) fixed.pars$kappa <- kappa
  if(fix.psiA) fixed.pars$psiA <- psiA
  if(fix.psiR) fixed.pars$psiR <- psiR
  .check.geoRparameters.values(list=fixed.pars, messages = messages.screen)
  if(cov.model == "matern" & all(kappa == 0.5)) cov.model <- "exponential"
  temp.list$cov.model <- cov.model
  if(cov.model == "powered.exponential")
    if(limits$kappa["upper"] > 2) limits$kappa["upper"] <- 2
  if(cov.model == "gencauchy")
    if(limits$kappa2["upper"] > 2) limits$kappa2["upper"] <- 2
  ##
  ## Likelihood method
  ##
#####
##### temporary code back compatibility for argument "method"
  lik.MET <- c("ML", "ml", "RML", "REML", "rml", "reml")
  MET <- pmatch(names(ldots), "method") == 1
  if(!is.na(MET) && any(MET) && (ldots[[which(MET)]] %in% lik.MET)){
    warning("argument \"method\" has changed and is now used as an argument to be passed to optim(). Use \"lik.method\" to define the likelihood method")
    lik.method <- lik.MET[pmatch(ldots[[which(MET)]], lik.MET)]
    ldots[which(as.logical(pmatch(names(ldots), "method", nomatch=0)))] <- NULL
  }
#####
  method.lik <- lik.method
  if(method.lik %in% c("REML","reml","rml","RML"))  method.lik <- "RML"
  if(method.lik %in% c("ML", "ml")) method.lik <- "ML"
  if(method.lik == "ML" & cov.model == "power")
    stop("\n\"power\" model can only be used with method.lik=\"RML\".\nBe sure that what you want is not \"powered.exponential\"")
  temp.list$method.lik <- method.lik
  ##
  ## setting coordinates, data and covariate matrices
  ##
  coords <- as.matrix(coords)
  data <- as.vector(data)
  n <- length(data)
  if((nrow(coords) != n) | (2*n) != length(coords))
    stop("\nnumber of locations does not match with number of data")
  if(missing(geodata))
    xmat <- trend.spatial(trend=trend, geodata=list(coords = coords, data = data))
  else xmat <- unclass(trend.spatial(trend=trend, geodata=geodata))
  xmat.contrasts  <- attr(xmat,"contrasts")
  xmat <- unclass(xmat)
  if(nrow(xmat) != n)
    stop("trend matrix has dimension incompatible with the data")
  .solve.geoR(crossprod(xmat))
  beta.size <- temp.list$beta.size <- dim(xmat)[2]
  ##
  ## setting a factor to indicate different realisations
  ##
  if(missing(realisations))
    realisations <- as.factor(rep(1, n))
  else{
    if(!missing(geodata)){
        real.name <- deparse(substitute(realisations))
        if(all(isTRUE(as.logical(real.name))))
          if(is.null(geodata$realisations)) stop("element realisation not available in the geodata object")
          else realisations <- geodata$realisations
      else{
        if(!is.null(geodata[[real.name]]))
          realisations <- geodata[[real.name]]
      }
    }
    if(length(realisations) != n)
      stop("realisations must be a vector with the same length of the data")
    realisations <- as.factor(realisations)
  }
  temp.list$realisations <- realisations
  nrep <- temp.list$nrep <- length(levels(realisations))
  ind.rep <- split(1:n, realisations)
  vecdist <- function(x){as.vector(dist(x))}
  ##
  ## Initial values for parameters
  ##
  ## have to consider transformation, residuals from trend etc
#  var.data <- mean(tapply(data, realisations, var))
#  d.max <- max(by(ap$coords, ap$realisations, function(x) max(dist(x))))
#  if(missing(ini.cov.pars))
#    ini.cov.pars <- expand.grid(var.data/2, 3*var.data/4, var.data)
  if(any(class(ini.cov.pars) == "eyefit")){
    init <- nugget <- kappa <- NULL
    for(i in 1:length(ini.cov.pars)){
      init <- drop(rbind(init, ini.cov.pars[[i]]$cov.pars))
      nugget <- c(nugget, ini.cov.pars[[i]]$nugget)
      if(cov.model == "gneiting.matern")
        kappa <- drop(rbind(kappa, ini.cov.pars[[i]]$kappa))
      else
        kappa <- c(kappa, ini.cov.pars[[i]]$kappa)
    }
    ini.cov.pars <- init
  }
  if(any(class(ini.cov.pars) == "variomodel")){
    nugget <- ini.cov.pars$nugget
    kappa <- ini.cov.pars$kappa
    ini.cov.pars <- ini.cov.pars$cov.pars
  }
  if(is.matrix(ini.cov.pars) | is.data.frame(ini.cov.pars)){
    ini.cov.pars <- as.matrix(ini.cov.pars)
    if(nrow(ini.cov.pars) == 1)
      ini.cov.pars <- as.vector(ini.cov.pars)
    else{
      if((cov.model != "pure.nugget") & (ncol(ini.cov.pars) != 2))
        stop("\nini.cov.pars must be a matrix or data.frame with 2 components: \ninitial values for sigmasq and phi")
    }
  }
  if(is.vector(ini.cov.pars)){
    if((cov.model != "pure.nugget") & (length(ini.cov.pars) != 2))
      stop("\nini.cov.pars must be a vector with 2 components: \ninitial values for sigmasq and phi")
  }
  ##
  ## Checking for multiple initial values for preliminar search of   
  ## best initial value
  ##
  if(is.matrix(ini.cov.pars) | (length(nugget) > 1) | (length(kappa) > 1) | (length(lambda) > 1) | (length(psiR) > 1) | (length(psiA) > 1)){
    if(messages.screen) cat("likfit: searching for best initial value ...")
    ini.temp <- matrix(ini.cov.pars, ncol=2)
    grid.ini <- as.matrix(expand.grid(sigmasq=unique(ini.temp[,1]), phi=unique(ini.temp[,2]), tausq=unique(nugget), kappa=unique(kappa), lambda=unique(lambda), psiR=unique(psiR), psiA=unique(psiA)))
    assign(".likGRF.dists.vec",  lapply(split(as.data.frame(coords), realisations), vecdist), pos=1)
    temp.f <- function(parms, coords, data, temp.list)
      return(loglik.GRF(geodata = geodata,
                        coords = coords, data = as.vector(data),
                        cov.model=temp.list$cov.model,
                        cov.pars=parms[1:2],
                        nugget=parms["tausq"], kappa=parms["kappa"],
                        lambda=parms["lambda"], psiR=parms["psiR"],
                        psiA=parms["psiA"], trend= trend,
                        method.lik=temp.list$method.lik,
                        compute.dists=FALSE,
                        realisations = realisations))
    grid.lik <- apply(grid.ini, 1, temp.f, coords = coords,
                      data = data, temp.list = temp.list)
    grid.ini <- grid.ini[(grid.lik != Inf) & (grid.lik != -Inf) & !is.na(grid.lik) & !is.nan(grid.lik),, drop=FALSE] 
    grid.lik <- grid.lik[(grid.lik != Inf) & (grid.lik != -Inf) & !is.na(grid.lik) & !is.nan(grid.lik)] 
    ini.temp <- grid.ini[which(grid.lik == max(grid.lik)),, drop=FALSE]
    if(all(ini.temp[,"phi"] == 0)) ini.temp <- ini.temp[1,, drop=FALSE]
    rownames(ini.temp) <- "initial.value"
    if(messages.screen){
      cat(" selected values:\n")
      print(rbind(format(ini.temp, digits=2), status=ifelse(c(FALSE, FALSE, fix.nugget, fix.kappa, fix.lambda, fix.psiR, fix.psiA), "fix", "est")))
      cat(paste("likelihood value:", max(grid.lik), "\n"))
    }
    dimnames(ini.temp) <- NULL
    ini.cov.pars <- ini.temp[1:2]
    nugget <- ini.temp[3]
    kappa <- ini.temp[4]
    lambda <- ini.temp[5]
    psiR <- ini.temp[6]
    psiA <- ini.temp[7]
    grid.ini <- NULL
    remove(".likGRF.dists.vec", pos=1)
  }
  ##
  tausq <- nugget
  ##
  ## Box-Cox transformation for fixed lambda
  ##
  if(fix.lambda) {
    if(abs(lambda - 1) < 0.0001) {
      temp.list$log.jacobian <- 0
      temp.list$z <- as.vector(data)
    }
    else {
      if(any(data <= 0))
        stop("Transformation option not allowed when there are zeros or negative data")
      Jdata <- data^(lambda - 1)
      if(any(Jdata <= 0))
        temp.list$log.jacobian <- log(prod(Jdata))
      else temp.list$log.jacobian <- sum(log(Jdata))
      Jdata <- NULL
      if(abs(lambda) < 0.0001)
        temp.list$z <- log(data)
      else temp.list$z <- ((data^lambda) - 1)/lambda
    }
  }
  else{
    temp.list$z <- as.vector(data)
    temp.list$log.jacobian <- NULL
  }
  ##
  ## Coordinates transformation for fixed anisotropy parameters
  ##
  if(fix.psiR & fix.psiA){
    if(psiR != 1 | psiA != 0)
      coords <- coords.aniso(coords, aniso.pars=c(psiA, psiR))
      assign(".likGRF.dists.vec", lapply(split(as.data.frame(coords), realisations), vecdist), pos=1)
    range.dist <- range(get(".likGRF.dists.vec", pos=1))
    max.dist <- max(range.dist)
    min.dist <- min(range.dist)
  }
  ##
  ##
  ##
  ini <- ini.cov.pars[2]
  ##  fixed.pars <- NULL
  lower.optim <- c(limits$phi["lower"])
  upper.optim <- c(limits$phi["upper"])
  fixed.values <- list()
  if(fix.nugget) {
    ##    fixed.pars <- c(fixed.pars, 0)
    fixed.values$tausq <- nugget
  }
  else {
    ini <- c(ini, nugget/ini.cov.pars[1])
    lower.optim <- c(lower.optim, limits$tausq.rel["lower"])
    upper.optim <- c(upper.optim, limits$tausq.rel["upper"])
  }
  if(fix.kappa){
    ##    fixed.kappa <- c(fixed.pars, kappa)
    fixed.values$kappa <- kappa
  }
  else {
    ini <- c(ini, kappa)
    lower.optim <- c(lower.optim, limits$kappa["lower"])
    upper.optim <- c(upper.optim, limits$kappa["upper"])
  }
  if(fix.lambda){
    ##    fixed.pars <- c(fixed.pars, lambda)
    fixed.values$lambda <- lambda
  }
  else {
    ini <- c(ini, lambda)
    lower.optim <- c(lower.optim, limits$lambda["lower"])
    upper.optim <- c(upper.optim, limits$lambda["upper"])
  }
  if(fix.psiR){
    ##    fixed.pars <- c(fixed.pars, psiR)
    fixed.values$psiR <- psiR
  }
  else {
    ini <- c(ini, psiR)
    lower.optim <- c(lower.optim, limits$psiR["lower"])
    upper.optim <- c(upper.optim, limits$psiR["upper"])
  }
  if(fix.psiA){
    ##    fixed.pars <- c(fixed.pars, psiA)
    fixed.values$psiA <- psiA
  }
  else {
    ini <- c(ini, psiA)
    lower.optim <- c(lower.optim, limits$psiA["lower"])
    upper.optim <- c(upper.optim, limits$psiA["upper"])
  }
  ## This must be here, after the previous ones:
  if(fix.nugget & nugget > 0){
    ## Warning: Inverting order here, ini will be now: c(phi,sigmasg)
    ini <- c(ini, ini.cov.pars[1])
    lower.optim <- c(lower.optim, limits$sigmasq["lower"])
    upper.optim <- c(upper.optim, limits$sigmasq["upper"])
    ##    fixed.pars <- c(fixed.pars, ini.cov.pars[1])
    ##    fixed.values$sigmasq <- 0
  }
  ##
  names(ini) <- NULL
  if(length(ini) == 1) justone <- TRUE
  else justone <- FALSE
  ##
  ip <- list(f.tausq = fix.nugget, f.kappa = fix.kappa,
             f.lambda = fix.lambda,
             f.psiR = fix.psiR, f.psiA = fix.psiA)
  ##
  npars <- beta.size + 2 + sum(unlist(ip)==FALSE)
  temp.list$coords <- coords
  temp.list$xmat <- split(as.data.frame(unclass(xmat)), realisations)
  temp.list$xmat <- lapply(temp.list$xmat, as.matrix)
  temp.list$n <- as.vector(unlist(lapply(temp.list$xmat, nrow)))
  ##
  ## Constant term in the likelihood
  ##
  temp.list$loglik.cte <- rep(0, nrep)
  for(i in 1:nrep){
    if(method.lik == "ML"){
      if(ip$f.tausq & (tausq > 0))
        temp.list$loglik.cte[i] <-  (temp.list$n[i]/2)*(-log(2*pi))
      else
        temp.list$loglik.cte[i] <-  (temp.list$n[i]/2)*(-log(2*pi) +
                                                        log(temp.list$n[i]) -1)
    }
    if(method.lik == "RML"){
      xx.eigen <- eigen(crossprod(temp.list$xmat[[i]]),
                        symmetric = TRUE, only.values = TRUE)
      if(ip$f.tausq & (tausq > 0))
        temp.list$loglik.cte[i] <- - ((temp.list$n[i]-beta.size)/2)*(log(2*pi)) +
          0.5 * sum(log(xx.eigen$values))
      else
        temp.list$loglik.cte[i] <-  - ((temp.list$n[i]-beta.size)/2)*(log(2*pi)) +
          ((temp.list$n[i]-beta.size)/2)*(log(temp.list$n[i]-beta.size)) -
            ((temp.list$n[i]-beta.size)/2) + 0.5 * sum(log(xx.eigen$values))
    }
  }
  ##  
  if(messages.screen) {
    cat("---------------------------------------------------------------\n")
    cat("likfit: likelihood maximisation using the function ")
    if(is.R()){if(justone) cat("optimize.\n") else cat("optim.\n")} else cat("nlminb.\n")
    cat("likfit: Use control() to pass additional\n         arguments for the maximisation function.")
    cat("\n        For further details see documentation for ")
    if(is.R()){if(justone) cat("optimize.\n") else cat("optim.\n")} else cat("nlminb.\n")        
    cat("likfit: It is highly advisable to run this function several\n        times with different initial values for the parameters.\n")
    cat("likfit: WARNING: This step can be time demanding!\n")
    cat("---------------------------------------------------------------\n")
  }
  ##
  ## Numerical minimization of the -loglikelihood
  ##
  if(length(ini) == 1){
    if(upper.optim == Inf) upper.optim <- 50*max.dist
    lik.minim <- do.call("optimize", c(list(.negloglik.GRF,
                                            lower=lower.optim,
                                            upper=upper.optim,
                                            fp=fixed.values,
                                            ip=ip, temp.list = temp.list), ldots))
    lik.minim <- list(par = lik.minim$minimum,
                      value = lik.minim$objective,
                      convergence = 0,
                      message = "function optimize used")      
  }
  else{
    MET <- pmatch(names(ldots), names(formals(optim)))
    if(is.na(MET) || all(names(formals(optim))[MET] != "method"))
      ldots$method <- "L-BFGS-B"
    if(!is.null(names(ldots))){
      names(ldots)[which(as.logical(pmatch(names(ldots), "method", nomatch=0)))] <- "method"
    }
    if(!is.null(ldots$method) && ldots$method == "L-BFGS-B"){
      ldots$lower <- lower.optim
      ldots$upper <- upper.optim
    }
    lik.minim <- do.call("optim", c(list(par = ini, fn = .negloglik.GRF,
                                         fp=fixed.values, ip=ip, temp.list = temp.list), ldots))
    ##      lik.minim <- optim(par = ini, fn = .negloglik.GRF, method=optim.METHOD
    ##                         lower=lower.optim, upper=upper.optim,
    ##                         fp=fixed.values, ip=ip, temp.list = temp.list, ...)
  }
  ##
  if(messages.screen) cat("likfit: end of numerical maximisation.\n")
  par.est <- lik.minim$par
  if(any(par.est < 0)) par.est <- round(par.est, digits=12)
  phi <- par.est[1]
  ##
  ## Values of the maximised likelihood
  ##
  if(is.R())
    loglik.max <- - lik.minim$value
  else
    loglik.max <- - lik.minim$objective
  ##
  ## Assigning values for estimated parameters
  ##
  if(ip$f.tausq & ip$f.kappa & ip$f.lambda & ip$f.psiR & !ip$f.psiA){
    psiA <- par.est[2]
  }
  if(ip$f.tausq & ip$f.kappa & ip$f.lambda & !ip$f.psiR & ip$f.psiA){
    psiR <- par.est[2]
  }
  if(ip$f.tausq & ip$f.kappa & ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
    psiR <- par.est[2]
    psiA <- par.est[3]
  }
  if(ip$f.tausq & ip$f.kappa & !ip$f.lambda & ip$f.psiR & ip$f.psiA){
    lambda  <- par.est[2]
  }
  if(ip$f.tausq & ip$f.kappa & !ip$f.lambda & ip$f.psiR & !ip$f.psiA){
    lambda  <- par.est[2]
    psiA <- par.est[3]
  }
  if(ip$f.tausq & ip$f.kappa & !ip$f.lambda & !ip$f.psiR & ip$f.psiA){
    lambda  <- par.est[2]
    psiR <- par.est[3]
  }
  if(ip$f.tausq & ip$f.kappa & !ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
    lambda  <- par.est[2]
    psiR <- par.est[3]
    psiA <- par.est[4]
  }
  if(ip$f.tausq & !ip$f.kappa & ip$f.lambda & ip$f.psiR & ip$f.psiA){
    kappa  <-  par.est[2]
  }
  if(ip$f.tausq & !ip$f.kappa & ip$f.lambda & ip$f.psiR & !ip$f.psiA){
    kappa  <-  par.est[2]
    psiA <- par.est[3]
  }
  if(ip$f.tausq & !ip$f.kappa & ip$f.lambda & !ip$f.psiR & ip$f.psiA){
    kappa  <-  par.est[2]
    psiR <- par.est[3]
  }
  if(ip$f.tausq & !ip$f.kappa & ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
    kappa  <-  par.est[2]
    psiR <- par.est[3]
    psiA <- par.est[4]
  }
  if(ip$f.tausq & !ip$f.kappa & !ip$f.lambda & ip$f.psiR & ip$f.psiA){
    kappa <-  par.est[2]
    lambda <- par.est[3]
  }
  if(ip$f.tausq & !ip$f.kappa & !ip$f.lambda & ip$f.psiR & !ip$f.psiA){
    kappa <-  par.est[2]
    lambda <- par.est[3]
    psiA <- par.est[4]
  }
  if(ip$f.tausq & !ip$f.kappa & !ip$f.lambda & !ip$f.psiR & ip$f.psiA){
    kappa <-  par.est[2]
    lambda <- par.est[3]
    psiR<- par.est[4]
  }
  if(ip$f.tausq & !ip$f.kappa & !ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
    kappa <-  par.est[2]
    lambda <- par.est[3]
    psiR<- par.est[4]
    psiA<- par.est[5]
  }
  if(!ip$f.tausq & ip$f.kappa & ip$f.lambda & ip$f.psiR & ip$f.psiA){
    tausq <- par.est[2]
  }
  if(!ip$f.tausq & ip$f.kappa & ip$f.lambda & ip$f.psiR & !ip$f.psiA){
    tausq <- par.est[2]
    psiA<- par.est[3]
  }
  if(!ip$f.tausq & ip$f.kappa & ip$f.lambda & !ip$f.psiR & ip$f.psiA){
    tausq <- par.est[2]
    psiR<- par.est[3]
  }
  if(!ip$f.tausq & ip$f.kappa & ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
    tausq <- par.est[2]
    psiR<- par.est[3]
    psiA<- par.est[4]
  }
  if(!ip$f.tausq & ip$f.kappa & !ip$f.lambda & ip$f.psiR & ip$f.psiA){
    tausq <- par.est[2]
    lambda <- par.est[3]
  }
  if(!ip$f.tausq & ip$f.kappa & !ip$f.lambda & ip$f.psiR & !ip$f.psiA){
    tausq <- par.est[2]
    lambda <- par.est[3]
    psiA <- par.est[4]
  }
  if(!ip$f.tausq & ip$f.kappa & !ip$f.lambda & !ip$f.psiR & ip$f.psiA){
    tausq <- par.est[2]
    lambda <- par.est[3]
    psiR <- par.est[4]
  }
  if(!ip$f.tausq & ip$f.kappa & !ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
    tausq <- par.est[2]
    lambda <- par.est[3]
    psiR <- par.est[4]
    psiA <- par.est[5]
  }
  if(!ip$f.tausq & !ip$f.kappa & ip$f.lambda & ip$f.psiR & ip$f.psiA){
    tausq <- par.est[2]
    kappa <-  par.est[3]
  }
  if(!ip$f.tausq & !ip$f.kappa & ip$f.lambda & ip$f.psiR & !ip$f.psiA){
    tausq <- par.est[2]
    kappa <-  par.est[3]
    psiA <- par.est[4]
  }
  if(!ip$f.tausq & !ip$f.kappa & ip$f.lambda & !ip$f.psiR & ip$f.psiA){
    tausq <- par.est[2]
    kappa <-  par.est[3]
    psiR <- par.est[4]
  }
  if(!ip$f.tausq & !ip$f.kappa & ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
    tausq <- par.est[2]
    kappa <-  par.est[3]
    psiR <- par.est[4]
    psiA <- par.est[5]
  }
  if(!ip$f.tausq & !ip$f.kappa & !ip$f.lambda & ip$f.psiR & ip$f.psiA){
    tausq <- par.est[2]
    kappa <-  par.est[3]
    lambda <- par.est[4]
  }
  if(!ip$f.tausq & !ip$f.kappa & !ip$f.lambda & ip$f.psiR & !ip$f.psiA){
    tausq <- par.est[2]
    kappa <-  par.est[3]
    lambda <- par.est[4]
    psiA <- par.est[5]
  }
  if(!ip$f.tausq & !ip$f.kappa & !ip$f.lambda & !ip$f.psiR & ip$f.psiA){
    tausq <- par.est[2]
    kappa <-  par.est[3]
    lambda <- par.est[4]
    psiR <- par.est[5]
  }
  if(!ip$f.tausq & !ip$f.kappa & !ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
    tausq <- par.est[2]
    kappa <-  par.est[3]
    lambda <- par.est[4]
    psiR <- par.est[5]
    psiA <- par.est[6]
  }
  ##
  if(fix.nugget & nugget > 0){
    sigmasq <- par.est[length(par.est)]
    if(sigmasq > 1e-12) tausq <- nugget/sigmasq
    check.sigmasq <- TRUE
  }
  else check.sigmasq <- FALSE
  ##
  ##
  ## Transforming data according to the estimated lambda (Box-Cox) parameter
  ##
  if(!fix.lambda) {
    if(abs(lambda - 1) < 0.0001) {
      log.jacobian.max <- 0
    }
    else {
      if(any(data^(lambda - 1) <= 0))
        log.jacobian.max <- log(prod(data^(lambda - 1)))
      else log.jacobian.max <- sum(log(data^(lambda - 1)))
      temp.list$z <- ((data^lambda)-1)/lambda
    }
  }
  else{
    log.jacobian.max <- temp.list$log.jacobian
  }
  data.rep <- split(temp.list$z, realisations)
  coords.rep <- split(as.data.frame(coords), realisations)
  coords.rep <- lapply(coords.rep, as.matrix)
  ##
  ## Transforming coords for estimated anisotropy (if the case)
  ##
  if(fix.psiR & fix.psiA)
    remove(".likGRF.dists.vec", pos=1)
  else{
    if(round(psiR, digits=6) != 1 | round(psiA, digits=6) != 0)
      coords <- coords.aniso(coords, aniso.pars=c(psiA, psiR))
    rangevecdist <- function(x){range(as.vector(dist(x)))}
    range.dist <- lapply(split(as.data.frame(coords), realisations), rangevecdist)
    range.dist <- range(as.vector(unlist(range.dist)))
    max.dist <- max(range.dist)
    min.dist <- min(range.dist)
  }      
#  gc(verbose=FALSE)
  ##
  ## Computing estimated beta and tausq/sigmasq (if the case)
  ##
  xivx <- matrix(0, ncol=beta.size, nrow=beta.size)
  xivy <- matrix(0, ncol=1, nrow=beta.size)
  yivy <- 0
  for(i in 1:nrep){
    ni <- temp.list$n[i]
    if((phi < 1e-12))
      V <- diag(x=(1+tausq), ni)
    else{
      if(check.sigmasq){
        if(sigmasq < 1e-12){
          if(!fix.nugget)
            V <- diag(x=(1+tausq), ni)
          else
            V <- diag(x=sqrt(tausq), ni)          
        }
        else
          V <- varcov.spatial(coords = coords.rep[[i]],
                              cov.model = cov.model,
                              kappa = kappa, nugget = tausq,
                              cov.pars = c(1, phi))$varcov
      }
      else
        V <- varcov.spatial(coords = coords.rep[[i]],
                            cov.model = cov.model,
                            kappa = kappa, nugget = tausq,
                            cov.pars = c(1, phi))$varcov
    }
    ivyx <- solve(V,cbind(data.rep[[i]],temp.list$xmat[[i]]))
    xivx <- xivx + crossprod(ivyx[,-1],temp.list$xmat[[i]])
    xivy <- xivy + crossprod(ivyx[,-1],data.rep[[i]])
    yivy <- yivy + crossprod(data.rep[[i]],ivyx[,1])
  }
  betahat <- .solve.geoR(xivx, xivy)
  res <- as.vector(temp.list$z - xmat %*% betahat)
  if(!fix.nugget | (nugget < 1e-12)){
    ssres <- as.vector(yivy - 2*crossprod(betahat,xivy) +
                       crossprod(betahat,xivx) %*% betahat)  
    if(method.lik == "ML")
      sigmasq <- ssres/n
    else
      sigmasq <- ssres/(n - beta.size)
  }
  if(fix.nugget){
    if(nugget > 0)
      tausq <- nugget
  }
  else tausq <- tausq * sigmasq
  betahat.var <- .solve.geoR(xivx)
  if(sigmasq > 1e-12) betahat.var <- sigmasq * betahat.var
#  if(!fix.nugget & phi < 1e-16){
#    tausq <- sigmasq + tausq
#    sigmasq <- 0
#  }
  ##
  ## Preparing output
  ##
  if((phi < 0.001*min.dist)){
    tausq <- tausq + sigmasq
    sigmasq <- 0
  }
  if((sigmasq < 1e-12)) phi <- 0
  ##
  n.model.pars <- beta.size + 7
  par.su <- data.frame(status=rep(-9,n.model.pars))
  ind.par.su <- c(rep(0, beta.size), ip$f.tausq, 0, 0, ip$f.kappa,
                  ip$f.psiR, ip$f.psiA,ip$f.lambda)
  par.su$status <- ifelse(ind.par.su,"fixed", "estimated")
  par.su$values <- round(c(betahat, tausq, sigmasq, phi, kappa, psiR, psiA, lambda), digits=4)
  if(beta.size == 1) beta.name <- "beta"
  else beta.name <- paste("beta", 0:(beta.size-1), sep="")
  row.names(par.su) <- c(beta.name, "tausq", "sigmasq", "phi", "kappa",
                             "psiR", "psiA", "lambda")
  par.su <- par.su[c((1:(n.model.pars-3)), n.model.pars-1, n.model.pars-2, n.model.pars),] 
  ##
  lik.results <- list(cov.model = cov.model,
                      nugget = tausq,
                      cov.pars=c(sigmasq, phi),
                      sigmasq = sigmasq,
                      phi = phi,
                      kappa = kappa,
                      beta = as.vector(betahat),
                      beta.var = betahat.var,
                      lambda = lambda,
                      aniso.pars = c(psiA = psiA, psiR = psiR),
                      tausq = tausq,
                      practicalRange = practicalRange(cov.model=cov.model,
                        phi = phi, kappa = kappa),
                      method.lik = method.lik, trend = trend,
                      loglik = loglik.max,
                      npars = npars,
                      AIC = -2 * (loglik.max - npars),
                      BIC = -2 * (loglik.max - 0.5 * log(n) * npars),
#                      residuals = res,
                      parameters.summary = par.su,
                      info.minimisation.function = lik.minim,
                      max.dist = max.dist,
                      trend = trend,
                      trend.matrix= xmat,
                      transform.info = list(fix.lambda = fix.lambda,
                        log.jacobian = log.jacobian.max))
  ##
  ## Likelihood results for the model without spatial correlation
  ##
  if(nospatial){
    if(fix.lambda){
      beta.ns <- .solve.geoR(crossprod(xmat), crossprod(xmat, temp.list$z))
      ss.ns <- sum((as.vector(temp.list$z - xmat %*% beta.ns))^2)
      if(method.lik == "ML"){
        nugget.ns <- ss.ns/n
        loglik.ns <- (n/2)*((-log(2*pi)) - log(nugget.ns) - 1) + temp.list$log.jacobian
      }
      if(method.lik == "RML"){
        nugget.ns <- ss.ns/(n-beta.size)
        loglik.ns <- ((n-beta.size)/2)*((-log(2*pi)) - log(nugget.ns) -1) +
          temp.list$log.jacobian
      }
      npars.ns <- beta.size + 1 + !fix.lambda
      lambda.ns <- lambda
    }
    else{
      if(is.R())
        lik.lambda.ns <- optim(par=1, fn = .negloglik.boxcox,
                               method = "L-BFGS-B",
                               lower = limits$lambda["lower"],
                               upper = limits$lambda["upper"],
                               data = data, xmat = xmat,
                               lik.method = method.lik)
      else
        lik.lambda.ns <- nlminb(par=1, fn = .negloglik.boxcox,
                                lower=limits$lambda["lower"],
                                upper=limits$lambda["upper"],
                                data = data, xmat = xmat,
                                lik.method = method.lik)
      lambda.ns <- lik.lambda.ns$par
      if(abs(lambda) < 0.0001) tdata.ns <- log(data)
      else tdata.ns <- ((data^lambda.ns)-1)/lambda.ns
      beta.ns <- .solve.geoR(crossprod(xmat),crossprod(xmat,tdata.ns))
      ss.ns <- sum((as.vector(tdata.ns - xmat %*% beta.ns))^2)
      if(is.R())
        value.min.ns <- lik.lambda.ns$value
      else
        value.min.ns <- lik.lambda.ns$objective
      if(method.lik == "ML"){
        loglik.ns <- (- value.min.ns)+ (n/2)*((-log(2*pi)) + log(n) - 1)
        nugget.ns <- ss.ns/n
      }
      if(method.lik == "RML"){
        nugget.ns <- ss.ns/(n-beta.size)
        loglik.ns <- (- value.min.ns)+ ((n-beta.size)/2)*((-log(2*pi)) +
                                                          log(n-beta.size) - 1)
      }      
      npars.ns <- beta.size + 1 + !fix.lambda
    }
    lik.results$nospatial <- list(beta.ns = beta.ns, variance.ns = nugget.ns,
                                  loglik.ns = loglik.ns, npars.ns = npars.ns,
                                  lambda.ns = lambda.ns, AIC.ns = -2 * (loglik.ns - npars.ns),
                                  BIC.ns = -2 * (loglik.ns - 0.5 * log(n) * npars.ns))
  }
  ##
  ## Assigning names to the components of the mean vector beta
  ##
  if(length(lik.results$beta.var) == 1)
    lik.results$beta.var <- as.vector(lik.results$beta.var)
  if(length(lik.results$beta) > 1){
    ##    if(inherits(trend, "formula") || (!is.null(class(trend)) && any(class(trend) == "trend.spatial")))
    if(inherits(trend, "formula") || (length(class(trend)) > 0 && any(class(trend) == "trend.spatial")))
      beta.names <- c("intercept", paste("covar", 1:(ncol(xmat)-1), sep = ""))
    else
      if(trend == "1st")
        beta.names <- c("intercept", "x", "y")
      else
        if(trend == "2nd")
          beta.names <- c("intercept", "x", "y", "x2", "xy", "y2")
    names(lik.results$beta) <- beta.names
  }
  ##
  ## Computing residuals and predicted values
  ## (isolated components of the model)
  ##
  if(components) {
    if(!fix.psiR & !fix.psiA)
      if(psiR != 1 | psiA != 0)
        coords <- coords.aniso(coords, aniso.pars=c(psiA, psiR))
    #coords.rep <- split(as.data.frame(coords), realisations)
    #res.rep <- split(res, realisations)
    trend.comp <- temp.list$z - res
    spatial.comp <- list()
    for(i in 1:nrep){
#      invcov <- varcov.spatial(coords = coords[ind.rep[[i]],], cov.model = cov.model, 
#                               kappa = kappa, nugget = tausq,
#                               cov.pars = c(sigmasq, phi), inv=TRUE)$inverse 
#      covmat.signal <- varcov.spatial(coords = coords[ind.rep[[i]],],
#                                      cov.model = cov.model, 
#                                      kappa = kappa, nugget = 0,
#                                      cov.pars = c(sigmasq, phi))$varcov
      spatial.comp[[i]] <- as.vector(varcov.spatial(coords = coords[ind.rep[[i]],],
                                                    cov.model = cov.model, 
                                                    kappa = kappa, nugget = 0,
                                                    cov.pars = c(sigmasq, phi))$varcov %*%
                                     varcov.spatial(coords = coords[ind.rep[[i]],],
                                                    cov.model = cov.model, 
                                                    kappa = kappa, nugget = tausq,
                                                    cov.pars = c(sigmasq, phi), inv=TRUE)$inverse %*%
                                     res[ind.rep[[i]]]) 
    }
    spatial.comp <- as.vector(unlist(spatial.comp))[as.vector(unlist(ind.rep))]
    predict.comp <- trend.comp + spatial.comp
    residual.comp <- as.vector(temp.list$z - predict.comp)
#    residual.std <- as.vector(invcov %*% residual.comp)
#    residual.trend.std <- as.vector(invcov %*% res)
    lik.results$model.components <-
      data.frame(trend = trend.comp, spatial = spatial.comp, residuals = residual.comp)
#    lik.results$s2.random <- (crossprod(res,invcov) %*% res)/(n - beta.size)
#    lik.results$s2 <- (crossprod(residual.comp,invcov) %*% residual.comp)/(n - beta.size)
  }
  ##
  lik.results$contrasts <- xmat.contrasts
  lik.results$call <- call.fc
  ##
  ## Assigning classes
  ##
  oldClass(lik.results) <- c("likGRF", "variomodel")
  ##
  ## Some warning messages about particular possible results
  ##
  if(messages.screen){
    if((lik.results$cov.pars[1] < (0.01 * (lik.results$nugget + lik.results$cov.pars[1])))& lik.results$cov.pars[2] > 0)
      cat("\nWARNING: estimated sill is less than 1 hundredth of the total variance. Consider re-examine the model excluding spatial dependence\n" )      
    if((lik.results$cov.pars[2] > (10 * max.dist)) & lik.results$cov.pars[1] > 0 )
      cat("\nWARNING: estimated range is more than 10 times bigger than the biggest distance between two points. Consider re-examine the model:\n 1) excluding spatial dependence if estimated sill is too low and/or \n 2) taking trends (covariates) into account\n" ) 
    if(((lik.results$cov.pars[2] < (0.1 * min.dist)) & (lik.results$cov.pars[1] > 0)) & lik.results$cov.pars[2] > 0)
      cat("\nWARNING: estimated range is less than 1 tenth of the minimum distance between two points. Consider re-examine the model excluding spatial dependence\n" ) 
  }
  ##
  attr(lik.results, "geodata") <- name.geodata
  return(lik.results)
}

".negloglik.GRF" <-
  function(pars, fp, ip, temp.list)
### pars : values for the parameters to be estimated
  ## sequence is c(phi, tausq, kappa, lambda, psiR, psiA, sigmasq)
### fixed pars: parameters considered fixed
### ind.pars : list indicating which are fixed and which are to be estimated
  ##
  ## Warning:
  ##  if fix.nugget = TRUE and nugget > 0 ,
  ## sigmasq should be passed and fp$nugget is the value of the nugget
  ## otherwise the RELATIVE nugget should be passed
{
  p <- temp.list$beta.size
  log.jacobian <- temp.list$log.jacobian
  ## Obligatory parameter:
  phi <- pars[1]
  ## Others
  if(ip$f.tausq){
    if(fp$tausq > 0){
      npars.min <- length(pars)
      sigmasq <- pars[npars.min]
    }
    else sigmasq <- 1
  }
  else sigmasq <- 1
  if(ip$f.tausq & ip$f.kappa & ip$f.lambda & ip$f.psiR & ip$f.psiA){
    tausq <- fp$tausq
    kappa <- fp$kappa
    lambda <- fp$lambda
    psiR <- fp$psiR
    psiA <- fp$psiA
  }
  if(ip$f.tausq & ip$f.kappa & ip$f.lambda & ip$f.psiR & !ip$f.psiA){
    tausq <- fp$tausq
    kappa <- fp$kappa
    lambda <- fp$lambda
    psiR <- fp$psiR
    psiA <- pars[2]
  }
  if(ip$f.tausq & ip$f.kappa & ip$f.lambda & !ip$f.psiR & ip$f.psiA){
    tausq <- fp$tausq
    kappa <- fp$kappa
    lambda <- fp$lambda
    psiR <- pars[2]
    psiA <- fp$psiA
  }
  if(ip$f.tausq & ip$f.kappa & ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
    tausq <- fp$tausq
    kappa <- fp$kappa
    lambda <- fp$lambda
    psiR <- pars[2]
    psiA <- pars[3]
  }
  if(ip$f.tausq & ip$f.kappa & !ip$f.lambda & ip$f.psiR & ip$f.psiA){
    tausq <- fp$tausq
    kappa <- fp$kappa
    lambda <- pars[2]
    psiR <- fp$psiR
    psiA <- fp$psiA
  }
  if(ip$f.tausq & ip$f.kappa & !ip$f.lambda & ip$f.psiR & !ip$f.psiA){
    tausq <- fp$tausq
    kappa <- fp$kappa
    lambda <- pars[2]
    psiR <- fp$psiR
    psiA <- pars[3]
  }
  if(ip$f.tausq & ip$f.kappa & !ip$f.lambda & !ip$f.psiR & ip$f.psiA){
    tausq <- fp$tausq
    kappa <- fp$kappa
    lambda <- pars[2]
    psiR <- pars[3]
    psiA <- fp$psiA
  }
  if(ip$f.tausq & ip$f.kappa & !ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
    tausq <- fp$tausq
    kappa <- fp$kappa
    lambda <- pars[2]
    psiR <- pars[3]
    psiA <- pars[4]
  }
  if(ip$f.tausq & !ip$f.kappa & ip$f.lambda & ip$f.psiR & ip$f.psiA){
    tausq <- fp$tausq
    kappa <- pars[2]
    lambda <- fp$lambda
    psiR <- fp$psiR
    psiA <- fp$psiA
  }
  if(ip$f.tausq & !ip$f.kappa & ip$f.lambda & ip$f.psiR & !ip$f.psiA){
    tausq <- fp$tausq
    kappa <- pars[2]
    lambda <- fp$lambda
    psiR <- fp$psiR
    psiA <- pars[3]
  }
  if(ip$f.tausq & !ip$f.kappa & ip$f.lambda & !ip$f.psiR & ip$f.psiA){
    tausq <- fp$tausq
    kappa <- pars[2]
    lambda <- fp$lambda
    psiR <- pars[3]
    psiA <- fp$psiA
  }
  if(ip$f.tausq & !ip$f.kappa & ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
    tausq <- fp$tausq
    kappa <- pars[2]
    lambda <- fp$lambda
    psiR <- pars[3]
    psiA <- pars[4]
  }
  if(ip$f.tausq & !ip$f.kappa & !ip$f.lambda & ip$f.psiR & ip$f.psiA){
    tausq <- fp$tausq
    kappa <- pars[2]
    lambda <- pars[3]
    psiR <- fp$psiR
    psiA <- fp$psiA
  }
  if(ip$f.tausq & !ip$f.kappa & !ip$f.lambda & ip$f.psiR & !ip$f.psiA){
    tausq <- fp$tausq
    kappa <- pars[2]
    lambda <- pars[3]
    psiR <- fp$psiR
    psiA <- pars[4]
  }
  if(ip$f.tausq & !ip$f.kappa & !ip$f.lambda & !ip$f.psiR & ip$f.psiA){
    tausq <- fp$tausq
    kappa <- pars[2]
    lambda <- pars[3]
    psiR <- pars[4]
    psiA <- fp$psiA
  }
  if(ip$f.tausq & !ip$f.kappa & !ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
    tausq <- fp$tausq
    kappa <- pars[2]
    lambda <- pars[3]
    psiR <- pars[4]
    psiA <- pars[5]
  }
  if(!ip$f.tausq & ip$f.kappa & ip$f.lambda & ip$f.psiR & ip$f.psiA){
    tausq <- pars[2]
    kappa <- fp$kappa
    lambda <- fp$lambda
    psiR <- fp$psiR
    psiA <- fp$psiA
  }
  if(!ip$f.tausq & ip$f.kappa & ip$f.lambda & ip$f.psiR & !ip$f.psiA){
    tausq <- pars[2]
    kappa <- fp$kappa
    lambda <- fp$lambda
    psiR <- fp$psiR
    psiA <- pars[3]
  }
  if(!ip$f.tausq & ip$f.kappa & ip$f.lambda & !ip$f.psiR & ip$f.psiA){
    tausq <- pars[2]
    kappa <- fp$kappa
    lambda <- fp$lambda
    psiR <- pars[3]
    psiA <- fp$psiA
  }
  if(!ip$f.tausq & ip$f.kappa & ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
    tausq <- pars[2]
    kappa <- fp$kappa
    lambda <- fp$lambda
    psiR <- pars[3]
    psiA <- pars[4]
  }
  if(!ip$f.tausq & ip$f.kappa & !ip$f.lambda & ip$f.psiR & ip$f.psiA){
    tausq <- pars[2]
    kappa <- fp$kappa
    lambda <- pars[3]
    psiR <- fp$psiR
    psiA <- fp$psiA
  }
  if(!ip$f.tausq & ip$f.kappa & !ip$f.lambda & ip$f.psiR & !ip$f.psiA){
    tausq <- pars[2]
    kappa <- fp$kappa
    lambda <- pars[3]
    psiR <- fp$psiR
    psiA <- pars[4]
  }
  if(!ip$f.tausq & ip$f.kappa & !ip$f.lambda & !ip$f.psiR & ip$f.psiA){
    tausq <- pars[2]
    kappa <- fp$kappa
    lambda <- pars[3]
    psiR <- pars[4]
    psiA <- fp$psiA
  }
  if(!ip$f.tausq & ip$f.kappa & !ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
    tausq <- pars[2]
    kappa <- fp$kappa
    lambda <- pars[3]
    psiR <- pars[4]
    psiA <- pars[5]
  }
  if(!ip$f.tausq & !ip$f.kappa & ip$f.lambda & ip$f.psiR & ip$f.psiA){
    tausq <- pars[2]
    kappa <- pars[3]
    lambda <- fp$lambda
    psiR <- fp$psiR
    psiA <- fp$psiA
  }
  if(!ip$f.tausq & !ip$f.kappa & ip$f.lambda & ip$f.psiR & !ip$f.psiA){
    tausq <- pars[2]
    kappa <- pars[3]
    lambda <- fp$lambda
    psiR <- fp$psiR
    psiA <- pars[4]
  }
  if(!ip$f.tausq & !ip$f.kappa & ip$f.lambda & !ip$f.psiR & ip$f.psiA){
    tausq <- pars[2]
    kappa <- pars[3]
    lambda <- fp$lambda
    psiR <- pars[4]
    psiA <- fp$psiA
  }
  if(!ip$f.tausq & !ip$f.kappa & ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
    tausq <- pars[2]
    kappa <- pars[3]
    lambda <- fp$lambda
    psiR <- pars[4]
    psiA <- pars[5]
  }
  if(!ip$f.tausq & !ip$f.kappa & !ip$f.lambda & ip$f.psiR & ip$f.psiA){
    tausq <- pars[2]
    kappa <- pars[3]
    lambda <- pars[4]
    psiR <- fp$psiR
    psiA <- fp$psiA
  }
  if(!ip$f.tausq & !ip$f.kappa & !ip$f.lambda & ip$f.psiR & !ip$f.psiA){
    tausq <- pars[2]
    kappa <- pars[3]
    lambda <- pars[4]
    psiR <- fp$psiR
    psiA <- pars[5]
  }
  if(!ip$f.tausq & !ip$f.kappa & !ip$f.lambda & !ip$f.psiR & ip$f.psiA){
    tausq <- pars[2]
    kappa <- pars[3]
    lambda <- pars[4]
    psiR <- pars[5]
    psiA <- fp$psiA
  }
  if(!ip$f.tausq & !ip$f.kappa & !ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
    tausq <- pars[2]
    kappa <- pars[3]
    lambda <- pars[4]
    psiR <- pars[5]
    psiA <- pars[6]
  }
  ##
  if(temp.list$print.pars){
    running.pars <- c(phi = phi, tausq = tausq, kappa =kappa, psiA = psiA, psiR = psiR, lambda = lambda)
    if(ip$f.tausq && fp$tausq > 0)
      running.pars <- c(sigmasq=sigmasq, running.pars)
    print(running.pars)
  }
  ##
  ## Absurd values
  ##
  if(kappa < 1e-04 | (tausq+sigmasq) < (.Machine$double.eps^0.5) |
     any(c(phi, tausq, sigmasq, kappa) < 0))
    return(.Machine$double.xmax^0.5)
  ##
  ## Anisotropy
  ##
  if(!ip$f.psiR | !ip$f.psiA){
    coords.c <- coords.aniso(temp.list$coords, aniso.pars=c(psiA, psiR))
    vecdist <- function(x){as.vector(dist(x))}
    assign(".likGRF.dists.vec", lapply(split(as.data.frame(coords.c),
                                             temp.list$realisations), vecdist), pos=1)
  }
  ##
  ## Box-Cox transformation
  ##
  if(!ip$f.lambda){
    if(abs(lambda - 1) < 0.0001) {
      log.jacobian <- 0
    }
    else {
      if(any(temp.list$z <= 0))
        stop("Transformation not allowed for zero or negative data")
      data <- temp.list$z^(lambda - 1)
      if(any(data <= 0)) log.jacobian <- log(prod(data))
      else log.jacobian <- sum(log(data))
      data <- NULL
    }
    if(abs(lambda) < 0.0001)
      data <- log(temp.list$z)
    else data <- ((temp.list$z^lambda) - 1)/lambda
  }
  else data <- temp.list$z
  data <- split(data, as.factor(temp.list$realisations))
  ##
  ## Computing likelihood
  ##
  sumnegloglik <- 0
  for(i in 1:temp.list$nrep){
    ## NOTE: Likelihood for Independent observations 
    ##       arbitrary criteria used here:
    ##       (phi < 1-e16) or (sigmasq < 1-e16)  ==> independence
    ##
    n <- temp.list$n[i]
    xmat <- temp.list$xmat[[i]]
    z <- data[[i]]
    if((phi < 1e-16) | (sigmasq < 1e-16)){
      if(ip$f.tausq)
        v <- list(varcov = diag(x=(tausq+sigmasq), n),
                  log.det.to.half = (n/2) * log(tausq+sigmasq))
      else
        v <- list(varcov = diag(x=(1+tausq), n),
                  log.det.to.half = (n/2) * log(1+tausq))
    }
    else
      v <- varcov.spatial(dists.lowertri = get(".likGRF.dists.vec", pos=1)[[i]],
                          cov.model = temp.list$cov.model, kappa=kappa,
                          nugget = tausq, cov.pars=c(sigmasq, phi),
                          det = TRUE)
    if(!is.null(v$crash.parms)) return(.Machine$double.xmax^0.5)
    ivx <- solve(v$varcov,xmat)
    xivx <- crossprod(ivx,xmat)
    betahat <- try(.solve.geoR(xivx,crossprod(ivx,z)), silent=TRUE)
    if(inherits(betahat, "try-error")){
      t.ei <- eigen(xivx, symmetric = TRUE)
#      if(exists("trySilent"))
        betahat <- try(crossprod(t(t.ei$vec)/sqrt(t.ei$val)) %*% crossprod(ivx,z), silent=TRUE)
#      else{
#        error.now <- options()$show.error.message
#        options(show.error.messages = FALSE)
#        betahat <- try(crossprod(t(t.ei$vec)/sqrt(t.ei$val)) %*% crossprod(ivx,z))
#        if(is.null(error.now) || error.now) options(show.error.messages = TRUE)        
#      }
    }
    if(inherits(betahat, "try-error"))
      stop("Covariates have very different orders of magnitude. Try to multiply and/or divide them to bring them to similar orders of magnitude") 
    res <- z - xmat %*% betahat
    ssres <- drop(crossprod(res,solve(v$varcov,res)))
    if(temp.list$method.lik == "ML"){
      if(ip$f.tausq & (tausq > 0))
        negloglik <- v$log.det.to.half +  0.5 * ssres
      else
        negloglik <- (n/2) * log(ssres) +  v$log.det.to.half
    }
    if(temp.list$method.lik == "RML"){
      if(length(as.vector(xivx)) == 1) {
        choldet <- 0.5 * log(xivx)
      }
      else {
        chol.xivx <- chol(xivx)
        choldet <- sum(log(diag(chol.xivx)))
      }
      if(ip$f.tausq & (tausq > 0))
        negloglik <- v$log.det.to.half +  0.5 * ssres + choldet
      else
        negloglik <- ((n-p)/2) * log(ssres) +  v$log.det.to.half + choldet
    }  
    negloglik <- negloglik - temp.list$loglik.cte[i]
    sumnegloglik <- sumnegloglik + negloglik
  }
  sumnegloglik <- sumnegloglik - log.jacobian
  if(sumnegloglik > (.Machine$double.xmax^0.5) | sumnegloglik == Inf | sumnegloglik == -Inf)
    sumnegloglik <- .Machine$double.xmax^0.5
  if(temp.list$print.pars)
    cat(paste("log-likelihood = ", -sumnegloglik, "\n"))
  return(sumnegloglik) 
}

"likfit.limits" <- "pars.limits"

"print.likGRF" <-
  function (x, digits = max(3, getOption("digits") - 3), ...)
{
  est.pars <- as.vector(x$parameters.summary[x$parameters.summary[,1]
  == "estimated",2])
  names.est.pars <- dimnames(x$parameters.summary[x$parameters.summary[,1] == "estimated",])[[1]]
  names(est.pars) <- names.est.pars
  cat("likfit: estimated model parameters:\n")
  print.default(format(est.pars, digits = digits), ...)
  cat(paste("Practical Range with cor=0.05 for asymptotic range:", format(x$practicalRange, ...)))
  cat("\n")
  ##  print(round(est.pars, digits=digits))
  cat("\nlikfit: maximised log-likelihood = ")
  cat(format(x$loglik, digits = digits))
  ##  cat(round(x$loglik, digits=digits))
  cat("\n")
  return(invisible())
}  

"summary.likGRF" <-
  function(object, ...)
{
  names.pars <- dimnames(object$parameters.summary)[[1]]
  summ.lik <- list()
  if(object$method.lik == "ML")
    summ.lik$method.lik <- "maximum likelihood"
  if(object$method.lik == "RML")
    summ.lik$method.lik <- "restricted maximum likelihood"
  summ.lik$mean.component <- object$beta
  names(summ.lik$mean.component) <- names.pars[1:length(object$beta)]
  summ.lik$cov.model <- object$cov.model
  summ.lik$spatial.component <- object$parameters.summary[c("sigmasq", "phi"),]
  summ.lik$spatial.component.extra <- object$parameters.summary[c("kappa", "psiA", "psiR"),]
  summ.lik$nugget.component <- object$parameters.summary[c("tausq"),, drop=FALSE]
  summ.lik$transformation  <- object$parameters.summary[c("lambda"),, drop=FALSE]
  summ.lik$practicalRange <-  object$practicalRange
  summ.lik$likelihood <- list(log.L = object$loglik, n.params = as.integer(object$npars),
                               AIC = object$AIC, BIC = object$BIC)
  summ.lik$ns.likelihood <- list(log.L = object$nospatial$loglik.ns, n.params = object$nospatial$npars.ns, AIC=object$nospatial$AIC.ns, BIC=object$nospatial$BIC.ns)
  summ.lik$estimated.pars <- dimnames(object$parameters.summary[object$parameters.summary[,1] == "estimated",])[[1]]
  summ.lik$call <- object$call
  oldClass(summ.lik) <- "summary.likGRF"
  return(summ.lik)
}

"print.summary.likGRF" <-
  function (x, digits = max(3, getOption("digits") - 3), ...)
{
  if(length(class(x)) == 0 || all(class(x) != "summary.likGRF"))
    stop("object is not of the class \"summary.likGRF\"")
  cat("Summary of the parameter estimation\n")
  cat("-----------------------------------\n")
  cat(paste("Estimation method:", x$method.lik, "\n"))
  cat("\n")
  ##
  ## Estimates of the model components
  ## Model: Y(x) = X\beta + S(x) + e 
  ##
  cat("Parameters of the mean component (trend):")
  cat("\n")
  print.default(round(x$mean.component, digits = digits), ...)
  ##  print(round(x$mean.component, digits=digits))
  cat("\n")
  ##
  cat("Parameters of the spatial component:")
  cat("\n")
  cat(paste("   correlation function:", x$cov.model))
  cat(paste("\n      (estimated) variance parameter sigmasq (partial sill) = ", format(x$spatial.component[1,2], digits=digits)))
  cat(paste("\n      (estimated) cor. fct. parameter phi (range parameter)  = ", format(x$spatial.component[2,2], digits=digits)))
  if(any(x$cov.model == c("matern", "powered.exponential",
     "cauchy", "gencauchy", "gneiting.matern"))){
    kappa <- x$spatial.component.extra["kappa",2]
    if(x$spatial.component.extra["kappa",1] == "estimated")
      cat(paste("\n      (estimated) extra parameter kappa =", format(kappa, digits=digits)))
    else{
      cat(paste("\n      (fixed) extra parameter kappa = ", kappa))
      if(x$cov.model == "matern" & (format(kappa, digits=digits)  == 0.5))
      cat(" (exponential)")
    }
  }
  cat("\n")
  ##
  aniso <-  x$spatial.component.extra[c("psiA", "psiR"),]
  psiApsiR <- x$spatial.component.extra[c("psiA", "psiR"),2]
  cat("   anisotropy parameters:")
  if(aniso["psiA",1] == "estimated")
    cat(paste("\n      (estimated) anisotropy angle =",
              format(psiApsiR[1], digits=digits),
              " (",format((psiApsiR[1]*360)/(2*pi), dig=1), "degrees )"))
  else
    cat(paste("\n      (fixed) anisotropy angle =", psiApsiR[1],
              " (",(psiApsiR[1]*360)/(2*pi), "degrees )"))
  if(aniso["psiR",1] == "estimated")
    cat(paste("\n      (estimated) anisotropy ratio =",
              format(psiApsiR[2], digits=digits)))
  else
    cat(paste("\n      (fixed) anisotropy ratio =", psiApsiR[2]))
  cat("\n")
  cat("\n")  
  cat("Parameter of the error component:")
  if(x$nugget.component[,1] == "estimated")
    cat(paste("\n      (estimated) nugget = ", format(x$nugget.component[,2], digits=digits)))
  else
    cat(paste("\n      (fixed) nugget =", x$nugget.component[,2]))
  cat("\n")
  cat("\n")
  cat("Transformation parameter:")
  cat("\n")
  lambda <- x$transformation[,2]
  if(x$transformation[,1] == "estimated")
    cat(paste("      (estimated) Box-Cox parameter =", format(lambda, digits=digits)))
  else{
    cat(paste("      (fixed) Box-Cox parameter =", lambda))
    if(abs(lambda - 1) <  0.0001) cat(" (no transformation)")
    if(abs(lambda) < 0.0001) cat(" (log-transformation)")
  }
  cat("\n")
  cat("\n")
  cat("Practical Range with cor=0.05 for asymptotic range:",
      format(x$practicalRange, ...))
  cat("\n")
  cat("\n")
  cat("Maximised Likelihood:")
  cat("\n")
  print(format(x$likelihood, digits=digits), ...)
  cat("\n")
  cat("non spatial model:")
  cat("\n")
  print(format(x$ns.likelihood, digits=digits), ...)
  cat("\n")
  cat("Call:")
  cat("\n")
  print(x$call)
  cat("\n")
  invisible(x)
}

"loglik.GRF" <-
  function(geodata, coords=geodata$coords, data=geodata$data,
           obj.model = NULL,
           cov.model="exp", cov.pars,
           nugget=0, kappa=0.5, lambda=1, psiR=1, psiA=0,
           trend="cte", method.lik="ML",
           compute.dists = TRUE, realisations = NULL)
{
  if(!is.null(obj.model)){
    if(!is.null(obj.model$cov.model)) cov.model <- obj.model$cov.model
    if(!is.null(obj.model$cov.pars)) cov.pars <- obj.model$cov.pars
    if(!is.null(obj.model$nugget)) nugget <- obj.model$nugget
    if(!is.null(obj.model$kappa)) kappa <- obj.model$kappa
    if(!is.null(obj.model$lambda)) lambda <- obj.model$lambda
    if(!is.null(obj.model$psiR)) psiR <- obj.model$psiR
    if(!is.null(obj.model$psiA)) psiA <- obj.model$psiA      
   if(!is.null(obj.model$trend)) trend <- eval(obj.model$trend)
  ## a resolver: problema em passando  trend
  }
  sigmasq <- cov.pars[1]
  phi <- cov.pars[2]
  if(method.lik == "REML" | method.lik == "reml" | method.lik == "rml") 
    method.lik <- "RML"
  if(method.lik == "ML" | method.lik == "ml")
    method.lik <- "ML"
  if(is.null(realisations))
    realisations <- as.factor(rep(1, length(data)))
  else
    realisations <- as.factor(realisations)
  nrep <- length(levels(realisations))
  ##
  ## Absurd values
  ##
  if(kappa < 1e-04) return(-(.Machine$double.xmax^0.5))
  if((nugget+sigmasq) < 1e-16) return(-(.Machine$double.xmax^0.5))
  ##
  ## Trend matrix
  ##
  if(missing(geodata))
    xmat <- unclass(trend.spatial(trend=trend, geodata = list(coords = coords, data = data)))
  else
    xmat <- unclass(trend.spatial(trend=trend, geodata = geodata))
  if (nrow(xmat) != nrow(coords)) 
    stop("coords and trend have incompatible sizes")
  beta.size <- ncol(xmat)
  xmat <- split(as.data.frame(unclass(xmat)), realisations)
  xmat <- lapply(xmat, as.matrix)
  ##
  ## Anisotropy
  ##
  vecdist <- function(x){as.vector(dist(x))}
  if(psiR != 1 | psiA != 0){
    coords.c <- coords.aniso(coords, aniso.pars=c(psiA, psiR))
    .likGRF.dists.vec <- lapply(split(as.data.frame(coords.c),
                                      as.factor(realisations)), vecdist)
  }
  else if(compute.dists) .likGRF.dists.vec <- lapply(split(as.data.frame(coords),
                                                           as.factor(realisations)), vecdist)
  ##
  ## Box-Cox transformation
  ##
  z <- data
  if(abs(lambda - 1) < 0.0001)
    log.jacobian <- 0
  else {
    if(any(z <= 0))
      stop("Transformation not allowed for zero or negative data")
    data <- z^(lambda - 1)
    if(any(data <= 0)) log.jacobian <- log(prod(data))
    else log.jacobian <- sum(log(data))
    data <- NULL
    if(abs(lambda) < 0.0001)
      data <- log(z)
    else data <- ((z^lambda) - 1)/lambda
  }
  data <- split(data, as.factor(realisations))
  ##
  ## Computing likelihood
  ##
  sumnegloglik <- 0
  for(i in 1:nrep){
    ## NOTE: Likelihood for Independent observations 
    ##       arbitrary criteria used here:
    ##       (phi < 1-e16) or (sigmasq < 1-e16)  ==> independence
    ##
    n <- length(data[[1]])
    if((phi < 1e-16) | (sigmasq < 1e-16)){
      V <- list(varcov = diag(x=(nugget+sigmasq), n),
                log.det.to.half = (n/2) * log(nugget+sigmasq))
    }
    else{
      V <- varcov.spatial(dists.lowertri = .likGRF.dists.vec[[i]],
                          cov.model = cov.model, kappa=kappa,
                          nugget = nugget, cov.pars=c(sigmasq, phi),
                          det = TRUE)
    }
    if(!is.null(V$crash.parms)){
      cat("varcov.spatial: improper matrix for following the given parameters:")
      print(V$crash.parms)
      stop()
    }
    ivx <- solve(V$varcov,xmat[[i]])
    xivx <- crossprod(ivx,xmat[[i]])
    betahat <- .solve.geoR(xivx, crossprod(ivx,data[[i]]))
    res <- data[[i]] - xmat[[i]] %*% betahat
    ssres <- drop(crossprod(res, solve(V$varcov,res)))
    if(method.lik == "ML"){
      negloglik <- (n/2)*(log(2*pi)) + V$log.det.to.half +  0.5 * ssres
    }
    if(method.lik == "RML"){
      choldet <- sum(log(diag(chol(xivx))))
      negloglik <- V$log.det.to.half +  0.5 * ssres + choldet
      xx.eigen <- eigen(crossprod(xmat[[i]]), symmetric = TRUE, only.values = TRUE)
      negloglik <- negloglik + ((n-beta.size)/2)*(log(2*pi)) - 0.5 * sum(log(xx.eigen$values))
    }
    sumnegloglik <- sumnegloglik + negloglik
  }
  sumnegloglik <- sumnegloglik - log.jacobian
  if(sumnegloglik > (.Machine$double.xmax^0.5))
    sumnegloglik <- .Machine$double.xmax^0.5
  return(as.vector(-sumnegloglik))
}


##
## ----------------
## Other functions
## ----------------
##
##

"likfit.nospatial" <-
  function(temp.list, ...)
{
  results <- list()
  z <- temp.list$z
  n <- temp.list$n
  beta.size <- temp.list$beta.size
  xmat <- temp.list$xmat
  txmat <- temp.list$txmat
  ixx <- solve(crossprod(xmat))
  if(temp.list$fix.lambda == FALSE){
    if (temp.list$minimisation.function == "nlm"){
      assign(".temp.lower.lambda",-2, pos=1)
      assign(".temp.upper.lambda", 2, pos=1)
      results <- nlm(.proflik.lambda, 1, ...)
      if(exists(".temp.lambda")){
        results$lambda <- get(".temp.lambda", pos=1)
        remove(".temp.lambda", pos=1, inherits = TRUE)
      }
      else{
        results$lambda <- results$estimate
      }
      remove(".temp.lower.lambda", ".temp.upper.lambda", inherits = TRUE, pos=1)
    }
    if (temp.list$minimisation.function == "nlmP"){
      results <- .nlmP(.proflik.lambda, 1, lower=-2, upper=2,...)  
      results$lambda <- results$estimate
    }
    if (temp.list$minimisation.function == "optim"){
      results <- optim(1, .proflik.lambda, method="L-BFGS-B", lower=-2, upper=2,...)
      results$minimum <- results$value
      results$lambda <- results$par
    }
    if(results$lambda == 1) {
      temp.list$log.jacobian <- 0
    }
    else {
      if(any(z <= 0))
        stop("Transformation option not allowed when there are zeros or negative data")
      if(any(z^(results$lambda - 1) <= 0))
        temp.list$log.jacobian <- log(prod(z^(results$lambda - 1)))
      else temp.list$log.jacobian <- sum(log(z^(results$lambda - 1)))
      if(results$lambda == 0)
        z <- log(z)
      else z <- ((z^results$lambda) - 1)/results$lambda
    }
  }
  else{
    results$lambda <- temp.list$lambda
    results$code <- 1
    if (temp.list$minimisation.function == "optim") results$convergence <- 0
  }
  ssres <- (z %*% (diag(n) - xmat %*%
                   solve(crossprod(xmat)) %*% txmat) %*% z)
  if(temp.list$method == "ML"){
    results$tausqhat <- ssres/n
    if(temp.list$fix.lambda)
      results$minimum <- as.vector(((n/2) * log(2 * pi) +
                          (n/2) * log(results$tausqhat) +
                          (n/2)  -
                          temp.list$log.jacobian))
  }
  if(temp.list$method == "RML") {
    results$tausqhat  <- (ssres/(n-beta.size))
    if(temp.list$fix.lambda){
      results$minimum <- as.vector((((n - beta.size)/2) * log(2 * pi) +
                          ((n - beta.size)/2) * log(results$tausqhat) +
                          (n/2) -
                          temp.list$log.jacobian
                          ))
    }
  }
  if (temp.list$minimisation.function == "optim") results$value <- results$minimum    
  return(results)
}

"loglik.spatial" <- function(pars)
{
  tausq <- pars[1]
  sigmasq <- pars[2]
  sill.total <- tausq + sigmasq
  phi <- pars[3]
  lambda <- pars[4]
  z <- .temp.list$z
  n <- .temp.list$n
  if(!isTRUE(.temp.list$fix.lambda)) {
    if(lambda == 1) {
      .temp.list$log.jacobian <- 0
    }
    else {
      if(any(z < 0))
        stop("Transformation option not allowed when there are zeros or negative data"
             )
      if(any(z^(lambda - 1) <= 0))
        .temp.list$log.jacobian <- log(prod(z^(lambda -
                                               1)))
      else .temp.list$log.jacobian <- sum(log(z^(lambda - 1)))
      if(lambda == 0)
        z <- log(z)
      else z <- ((z^lambda) - 1)/lambda
    }
  }
  beta.size <- .temp.list$beta.size
  kappa <- .temp.list$kappa
#  covinf <- varcov.spatial(dists.lowertri = .temp.list$
#                           dists.lowertri, cov.model = .temp.list$cov.model,
#                           kappa = kappa, nugget = tausq,
#                           cov.pars = c(sigmasq, phi), scaled = FALSE,
#                           inv = TRUE, det = TRUE,
#                           only.inv.lower.diag = TRUE)
  COV <- varcov.spatial(dists.lowertri = .temp.list$
                        dists.lowertri, cov.model = .temp.list$cov.model,
                        kappa = kappa, nugget = tausq,
                        cov.pars = c(sigmasq, phi), scaled = FALSE,
                        det = TRUE)
#  xix <- as.double(rep(0, beta.size*beta.size))
#  xix <- .C("bilinearform_XAY",
#            as.double(covinf$lower.inverse),
#            as.double(covinf$diag.inverse),
#            as.double(as.vector(.temp.list$xmat)),
#            as.double(as.vector(.temp.list$xmat)),
#            as.integer(beta.size),
#            as.integer(beta.size),
#            as.integer(n),
#            res = xix, PACKAGE = "geoR")$res
#  attr(xix, "dim") <- c(beta.size, beta.size)
  xix <- crossprod(.temp.list$xmat, solve(COV$varcov,
                                           .temp.list$xmat))
  if(length(as.vector(xix)) == 1) {
    choldet <- 0.5 * log(xix)
  }
  else {
    choldet <- sum(log(diag(chol(xix))))
  }
  iy <- solve(COV$varcov, z)
  xiy <- crossprod(.temp.list$xmat, iy)
  beta.hat <- drop(solve(xix,xiy))
  yiy <- drop(crossprod(z, iy))
  ssresmat <- as.vector(yiy - 2*crossprod(beta.hat,xiy) +  beta.hat %*% xix %*% beta.hat)
  if(.temp.list$method == "ML") {
    loglik <- ( - (n/2) * log(2 * pi) -
               COV$log.det.to.half -
               0.5 * ssresmat + 
               .temp.list$log.jacobian)
  }
  if(.temp.list$method == "RML") {
    xx.eigen <- eigen(crossprod(.temp.list$xmat), symmetric = TRUE, only.values = TRUE)
    loglik <- ( - ((n - beta.size)/2) * log(2 * pi) +
               0.5 * sum(log(xx.eigen$values)) -
               COV$log.det.to.half -
               (0.5) * ssresmat -
               choldet + .temp.list$log.jacobian)
  }
  return(as.vector(loglik))
}

"logLik.likGRF" <- function(object, ...){
  val <- object$loglik
  attr(val, "df") <- object$npars
  attr(val, "method") <- object$method.lik
  class(val) <- "logLik"
  return(val)
}

"fitted.likGRF" <- 
  function(object, spatial = TRUE, ...)
{
  if(is.null(object$model.components)){
    object$call$components <- TRUE
    object$call$messages <- FALSE
    object$call$ini.cov.pars <- object$cov.pars
    object$call$nugget <- object$nugget
    object$call$lambda <- object$lambda
    cat("be patient ... this function currently require calling likfit again\n")
    object <- eval(object$call)
  }
  if(spatial)
    return(rowSums(object$model.components[,c("trend", "spatial")]))
  else
    return(object$model.components[,"trend"]) 
}

"resid.likGRF" <- "residuals.likGRF" <-
  function(object, spatial = FALSE, ...)
{
  if(is.null(object$model.components)){
    object$call$components <- TRUE
    object$call$messages <- FALSE
    object$call$ini.cov.pars <- object$cov.pars
    object$call$nugget <- object$nugget
    object$call$lambda <- object$lambda
    cat("be patient ... this function currently require calling likfit again\n")
    object <- eval(object$call)
  }
  if(spatial)
    return(rowSums(object$model.components[,c("spatial","residuals")]))
  else
    return(object$model.components[,"residuals"])
}

