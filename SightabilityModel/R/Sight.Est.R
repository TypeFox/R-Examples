Sight.Est <-
function(form, sdat=NULL, odat, sampinfo, method="Wong", logCI=TRUE, alpha=0.05, Vm.boot = FALSE, nboot = 1000, bet = NULL, varbet = NULL){

  # form = model formula (for fitting logistic regression model to sightability dataset)
  # sdat = sightability dataset containing response and sightability covariates
  # odat = dataset containing observed groups, sample unit ids, stratum identifiers, and sampling rates
  # sampinfo = data frame with sampling information (number of sample (nh) and population units (Nh) in each stratum
  # Method = method of variance estimation 
  # logCI = should the confidence interval be constructed under the assumption that tau^ is lognormally distributed
  # alpha = type I error rate for confidence interval construction
  # Vm.boot = should a bootstrap be used to estimate cov(theta[i,j],theta[i',j']), var/cov matrix of the expansion factors for detection
  # nboot = number of bootstraps if Vm.boot = TRUE
  # beta=regression parameters (can be passed, rather than using Sight.Est to fit the model)
  # varbet= var/cov matrix for regression parameter estimates (can be passed, rather than using Sight.Est to fit the model)

  # Check arguments
    if(is.null(sdat) & (is.null(bet) | is.null(varbet) )){
      stop("Need to specify sightability data set or beta^ and var(beta^)")
    }   
    if(method %in%c("Wong", "SS") != TRUE){
      stop("Method must be Wong or SS")
    }
    if(is.null(bet)){   
      sight.model <- glm(form, family = binomial(), data = sdat)
      beta <- coef(sight.model)
      varbet <- summary(sight.model)$cov.unscaled
    }else{
      sight.model <- NULL
      sight.model$bet <- bet 
      sight.model$varbet <- varbet
      sight.model$note <- "User supplied regression model"
      beta <- bet
      varbet <- varbet
    }  
    if(sum(c("stratum", "subunit", "total")%in%names(odat)) != 3){
      stop("Need to have variables stratum, subunit,  and total in observational survey dataset.  These are CASE-SENSITIVE")
    
    }  
    if(sum(odat$total == 0) > 0){
       print("Dropping records with 0 animals in the observational data set")
       odat <- odat[odat$total > 0, ]
    }   
    if(sum(c("nh","Nh", "stratum")%in%names(sampinfo)) != 3){
      stop("Need to have variables nh, Nh, and stratum in dataset containing sampling information")
    }  

 # Make sure all stratum names are in both observational and sampling information data sets
    s1 <- sort(unique(odat$stratum))
    s2 <- sort(unique(sampinfo$stratum))
    if(length(s1) != length(s2) | sum(s1%in%s2) != length(s2)){
      stop("The same list of stratum must be present in both both observational and sampling information data sets")
    }  
 # create sampling variables
    nh <- sampinfo$nh
    Nh <- sampinfo$Nh
  
 # Make sure strata are numbered 1:h
    sampinfo$stemp <- 1:nrow(sampinfo)
  
 # Want srates to be same legnth as total vector
    sampinfo$samp.rates <- sampinfo$nh/sampinfo$Nh
    odat <- merge(odat, sampinfo, by.x = "stratum", by.y = "stratum")  
  
 #pull off other vairables
    srates <- odat$samp.rates 
    stratum <- odat$stemp
    subunit <- odat$subunit
    total <- odat$total
   
    if(is.null(bet)){ 
      tempnm <- terms(form, data = sdat)
    }else{tempnm <- terms(form, data = odat)}
    tempnm2 <- attr(tempnm, "term.labels")
    if(sum(tempnm2%in%names(odat)) != length(tempnm2)){
       print("The exact same names need to be used for covariates in the sightability and operational survey datasets")
    }
    covars <- odat[, tempnm2]  
 
 #Bootstrap for Cov(theta)? 
    if(Vm.boot == TRUE){
      bets <- matrix(NA, nboot, length(beta))
      varbets <- array(NA,dim=c(nboot, length(beta), length(beta)))
      nsdat <- nrow(sdat)
      for(i in 1:nboot){
        bdat <- sdat[sample(1:nsdat, replace = TRUE), ] # sample observations with replacement (non-parametric bootstrap)
        fit <- glm(form, family = binomial, data = bdat)  # Fit model to bootstrap data
        bets[i,] <- coef(fit)  # pull of coefficients
        varbets[i,,] <- summary(fit)$cov.unscaled  # pull off variance/covariance matrix
      }
      smat <- covtheta(total, srates, stratum, subunit, covars, bets, varbets, nboot)
    }
    else{smat <- NULL}  
    if(method == "Wong"){
      est <- Wong.est(total, srates, nh, Nh, stratum, subunit, covars, beta, varbet, smat)
    }
    if(method == "SS"){
      est <- SS.est(total, srates, nh, Nh, stratum, subunit, covars, beta, varbet, smat)
    }    
    out <- NULL
    out$call <- match.call()
    out$sight.model <- sight.model
    out$est <- est$est
    out$var.method <- est$var.method 
    if(is.null(sdat) != TRUE){out$sdat <- sdat}
      out$odat <- odat
      out$samp <- sampinfo[,c("stratum", "nh", "Nh")]
      out$CI.method <- "lognormal"
      if(logCI == FALSE){  out$CI.method <- "normal"}
      out$alpha <- alpha
      class(out) <- "sightest"
      out
}
