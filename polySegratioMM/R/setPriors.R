`setPriors` <-
 function(model, type.prior=c("strong","vague","strong.tau","strong.s",
                   "specified"),
          mean.vague=1e-1, prec.vague=1e-1, 
          A.vague=1e-1, B.vague=1e-1,
          prec.strong=400, n.individuals=200,
          reffect.A=44, reffect.B=0.8,
          M.sd=0.025, STRONG.PREC=c(0.025, 0.975), UPPER=0.995, 
          PREC.INT=0.2, params=NULL, segRatio=NULL)
{

  ## NB: version 0.3-1  put in arguments for quartiles for strong priors for
  ## precision paraneters A, B : default set to c(0.025, c0.975)

  type <- match.arg(type.prior)

  n.comp <- model$n.components
  if (n.comp<2) stop("Number of components must be 2 or more")
  
  ##  vague priors 

  if (type=="vague") {
    M <- rep(mean.vague, n.comp)
    if (model$equal.variances) {
      prec <- prec.vague
      tau <- prec.vague
      A <- A.vague
      B <- B.vague
    } else {
      prec <- rep(prec.vague, n.comp)
      tau <- prec
      A <- rep(A.vague, n.comp)
      B <- rep(B.vague, n.comp)
    }
    if (model$random.effect){
      reffect.A <- A.vague
      reffect.B <- B.vague
    }
  }
  
  if (length(segRatio)==0 & (type=="strong" | type=="strong.tau" |
              type=="strong.s" )) {

    ## means of Normal components
    MM <- model$E.segRatio$ratio[1:n.comp]
    M <- gtools::logit(MM) # precision set below
    
    ## precisions of means of Normal components
    ##   NB: at least last one problematic so use UPPER as an upper bound
    limits <- cbind(low= MM - 2*M.sd, upp= MM + 2*M.sd)
    limits[limits>UPPER] <- UPPER  # for when things go over UPPER set the limit
    logit.limits <- gtools::logit(limits)
    logit.sd <- (logit.limits[,2]-logit.limits[,1])/4  # 4 sd's within the 95% CI 
    tau <- 1/logit.sd^2

    ## A and B for priors of precisions of Normal components
    ##    NB: at least last one problematic so limit to 
    ## sd <- (gtools::logit(pmin(MM + M.sd, UPPER))-gtools::logit(MM - M.sd))/4
    n.taus <- max(model$ploidy/2 - 1, n.comp) # ugh??? just n.comp?

    sd.low <- qbinom(p=0.025, size=n.individuals, prob=MM)/n.individuals
    sd.upp <- qbinom(p=0.975, size=n.individuals, prob=MM)/n.individuals
    sd.logit <- cbind(low=sd.low, upp=sd.upp)
    sd.logit[sd.logit>UPPER] <- UPPER
    sd.logit <- gtools::logit(sd.logit)

    sd.0 <- (sd.logit[,2]-sd.logit[,1])/4
    tau.0 <- 1/sd.0^2

    if (type=="strong") {
      B <-  rep((1+PREC.INT)^2 * (1-PREC.INT)^2, length(tau.0))
      A <- tau.0*B
    } 

    if (type=="strong.tau") {
      A <- rep(4/PREC.INT^2, length(tau.0))
      B <- A/tau.0
    } 
    
    if (type=="strong.s") {
      C <- PREC.INT^2 / (1+PREC.INT)^4 / (1-PREC.INT)^4
      B <- tau.0^3 * C
      A <- tau.0 * B
    } 

    if (model$equal.variances) {
      tau <- mean(tau)
      prec <- prec.strong                  # default: 400 = 1/0.05^2 NOT USED
      A <- A[1]
      B <- B[1]
    } else {
      ## prec <- rep(prec.strong, n.comp)     # default: 400 = 1/0.05^2
      prec <- rep(prec.strong, n.comp)
    }
  }

  if (type=="specified"){  # rudimentary error checking - watch out!
    if (length(params)==0) stop("params must be specified")
    M <- params$M; prec <- params$prec; A <- params$A; B <- params$B
    if (length(params$reffect.A)==1) reffect.A <- params$reffect.A
    if (length(params$reffect.B)==1) reffect.B <- params$reffect.B
  }
  
  ## bugs code for priors

  bc <- c(paste("## Priors:",toupper(type)),
          "P[]        ~ ddirch(alpha[]);    # prior for mixing proportion")
  bc <- c(bc, paste("for (j in 1:",model$n.components,"){",sep=""))
  bc <- c(bc, "  alpha[j]  <- 1;                  # vague  uniform priors",
          "}")
  bc <- c(bc,paste("mu[ 1 ]  ~ dnorm(0,",tau[1],");        #",type,"priors"))
  bc <- c(bc, paste("mu[",2:n.comp,"] <- mu[",1:(n.comp-1),"] + theta[",
                    1:(n.comp-1),"];"))
  if (model$equal.variances) {
    bc <- c(bc, paste("theta[", 1:(n.comp-1),"] ~ dnorm(",diff(M),",",
                      tau,")  T(0,) ;"))
    bc <- c(bc, paste("tau ~ dgamma(",A,",",B,");"))
  } else {
    bc <- c(bc, paste("theta[", 1:(n.comp-1),"] ~ dnorm(",diff(M),",",
                      tau[2:n.comp],")  T(0,) ;"))
    bc <- c(bc, paste("tau[",1:n.comp,"] ~ dgamma(",A,",",B,");") )
  }
  if (model$random.effect) {
    bc <- c(bc, paste("taub ~ dgamma(",reffect.A,",",reffect.B,");"))
    bc <- c(bc, "sigmab <- 1/sqrt(taub);")
  }
  bc <- c(bc, "}")

  ## return S3 class 

  res <- list(type=type, bugs.code=bc, random.effect=model$random.effect,
              equal.variances=model$equal.variances,
              params=list(logit.means=M, logit.prec=tau, A=A, B=B),
              call=match.call())
  if (model$random.effect){
    res$params$reffect.A=reffect.A
    res$params$reffect.B=reffect.B
  }
  oldClass(res) <- "priorsSegratioMM"

  ## if we've got this far - it would be nice to relabel the Priors: comment in
  ## model$bugs.code to be mode informative but I can't get it to work

  ## do this in calling environment???
  ## repline <- grep("Priors: ",model$bugs.code)
  ## replacement <-  paste((strsplit(model$bugs.code[repline],"Priors: "))[[1]][1],
  ##        "Priors: ", toupper(type))
###  parent.frame(n)
#  get(deparse(substitute(model)), envir=parent.frame())$bugs.code[repline] <- replacement # BUT see warning Venables pp 62
##  (eval(model,envir=sys.frame(parent.frame()), enclos))$bugs.code[repline] <- replacement

##  assign(model, envir=sys.frame(sys.parent()))$bugs.code[repline] <- replacement
 ## eval(model, envir=sys.frame(sys.parent())$bugs.code[repline] <- replacement
  
  return(res)
}

