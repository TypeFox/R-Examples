trtsel <-
function(event, trt, marker = NULL, data, fittedrisk.t0 = NULL, fittedrisk.t1 = NULL, thresh=0, study.design = "randomized cohort", cohort.attributes = NULL, marker.bounds = NULL, link = "logit", default.trt = "trt all" ){

  if(!is.data.frame(data)){stop('data must be a data.frame')}
  
  tmpnames <- c(event, trt, marker, fittedrisk.t0, fittedrisk.t1)
  if(!all(is.element(tmpnames, names(data)))) stop(paste("'", tmpnames[which(!is.element(tmpnames, names(data)))], "' cannot be found in data.frame provided", sep = ""))
    
  mycomplete <- complete.cases(data[,tmpnames]); 
  
  #check for missing data and throw it out, print a warning
  if(nrow(data)!=sum(mycomplete)){
    warning(paste(nrow(data)-sum(mycomplete), "observation(s) were removed due to missing data \n New sample size is now:", sum(mycomplete)))
    data <- data[mycomplete,]
    
  }
  ## Error Checking
  
  #check event
  if(!is.numeric(data[[event]]) | !all(is.element(unique(data[[event]]), c(0,1)))) stop( "event must be a numeric vector with elements 1 or 0")
 
  if(!is.element(link, c("logit", "probit", "cauchit", "log", "cloglog"))) stop("link must be one of ''logit'', ''probit'', ''cauchit'', ''log'', ''cloglog''")

  if(!is.numeric(data[[trt]]) | !all(is.element(unique(data[[trt]]), c(0,1)))) stop( "trt must be a numeric vector with elements 1 or 0") 
  
 # if(!is.numeric(data[[marker]])) stop( "marker must be a numeric") 
  if(length(marker) >1){ stop("only a single marker is allowed")}
  if(!is.element(default.trt, c("trt all", "trt none"))){ stop( "default.trt must be either 'trt all' or 'trt none'")}
  ## End Error Checking

  d <- thresh
#find out which bootstrapping functions to use based on type
  if( substr(study.design,1,4)  == "rand" ) { 
   
    boot.sample <- boot.sample.cohort 
    get.F <- get.F.cohort
    get.summary.measures <- get.summary.measures.cohort

    if(length(cohort.attributes) >0) warning("study.design = ''randomized cohort'', but cohort.attributes assigned: cohort.attributes will be ignored"); cohort.attributes=NULL;
     #just passing null value here 
    

  }else if( substr(study.design, 1, 4) =="nest") { 
    
    if(link!="logit") warning("when study.design is ''nested case-control'' only link=''logit'' is allowed, setting link = ''logit''"); link = "logit"
    boot.sample <- boot.sample.case.control 
    get.F <- get.F.case.control
    get.summary.measures <- get.summary.measures.case.control

    if(length(cohort.attributes) != 4){ 
    
    stop("cohort.attributes not specified correctly, when study.design=''nested case-control'': \n  
          cohort.attributes = c(cohort sample size, 
          proportion treated in cohort (Pr(trt==1)),
          adverse event prevalance in cohort (Pr(event==1)),
          fraction of cases sampled from cohort)\n")
    }

    if(length(trt) > cohort.attributes[1]) stop("Sub-cohort sample size larger than input cohort sample size, please check cohort.attributes[1]")

    if( any(cohort.attributes[-1]<0 ) | any(cohort.attributes[-1] > 1)) stop("Probabilities in cohort.attributes are not in (0,1), please check cohort.attributes")

   # if(length(cohort.attributes) == 3) cohort.attributes <- c(cohort.attributes, 1) # adding f = 1 by default
    #cohort attributes is c(N, Pr(trt = 1), Pr(event = 1) ) 

  }
  else if( substr(study.design, 1, 5) =="strat") { 
    if(link!="logit") warning("when study.design is ''stratified nested case-control'' only link=''logit'' is allowed, setting link = ''logit''"); link = "logit"
    boot.sample <- boot.sample.stratified.case.control
    get.F <- get.F.stratified.case.control  
    get.summary.measures <- get.summary.measures.stratified.case.control
   # if(length(cohort.attributes) == 5) cohort.attributes <- c(cohort.attributes, 1, 1)
   if(length(cohort.attributes) != 6){ 
    
    stop("cohort.attributes not specified correctly, when study.design=''nested case-control'': \n  
          cohort.attributes = c(cohort sample size, 
                              Pr(trt==0 & event==0) in cohort, 
                              Pr(trt==0 & event==1) in cohort, 
                              Pr(trt==1 & event==0) in cohort, 
                              fraction of cases with trt == 0 sampled from cohort, 
                              fraction of cases with trt == 1 sampled from cohort )\n ")
    }
    
    ca <- cohort.attributes
    cohort.attributes = c(ca[1], ca[2], ca[3], ca[4], 1-(ca[2]+ca[3]+ca[4]), ca[5], ca[6])


    if(length(trt) > ca[1]) stop("Sub-cohort sample size larger than input cohort sample size, please check cohort.attributes[1]")
    if( any(cohort.attributes[-1]<0) | any(cohort.attributes[-1] > 1)) stop("Probabilities in cohort.attributes are not in (0,1), please check cohort.attributes. \n Is Pr(trt==0 & event==0) + Pr(trt==0 & event==1) + Pr(trt==1 & event==0) > 1?")

  }
  else { stop("study.design not specified correctly, must be one of ''randomized cohort'', ''nested case-control'', or ''stratified nested case-control''") }

  
  functions <- list("boot.sample" = boot.sample, "get.F" = get.F, "get.summary.measures" = get.summary.measures)

  rho = cohort.attributes

  #make sure Pr(event = 1 | trt = 1) > Pr(event = 1 | trt = 0), if not, use T.star = 1-trt and give a warning


  # if( mean(event[trt==1]) > mean(event[trt==0]) & allow.switch ) {
  #  warning( "   Function assumes Pr(event = 1 | trt = 1) < Pr(event = 1 | trt = 0)\n   Redefining trt <- (1-trt)\n")
  #  trt <- 1-trt 
  # }
  
  event = data[[event]]
  trt = data[[trt]]
  if(!is.null(marker)) marker = data[[marker]]
  
  if(!is.null(fittedrisk.t0)) {
    fitted_risk_t0 = data[[fittedrisk.t0]]
    if(!is.null(marker)) warning("fitted risks provided: marker data will be ignored")
    marker <- NULL
    link <- "risks_provided"
    if(is.null(fittedrisk.t1)) stop("must provide fitted risk for trt = 1 as well")
    if(any(fitted_risk_t0 > 1) | any(fitted_risk_t0 <0)) stop("fitted risks for trt = 0 are outside of bounds (0,1)")
    
  }
  
  if(!is.null(fittedrisk.t1)){
    
   fitted_risk_t1 = data[[fittedrisk.t1]]
   if(is.null(fitted_risk_t0)) stop("must provide fitted risk for trt = 0 as well")
   if(any(fitted_risk_t1 > 1) | any(fitted_risk_t1 <0)) stop("fitted risks for trt = 1 are outside of bounds (0,1)")
   
  }  

  
  # model.fit
  #returns null if risks are provided
  coef <- get.coef( event, trt, marker, study.design, rho, link)

  model.fit <- list( "coefficients" = coef, "cohort.attributes" = rho, "study.design" = study.design, "marker.bounds" = marker.bounds, "link" = link, "thresh" = d)
  
  
  # derived.data

  #now that we allow for different link functions for "randomized cohorts", we need to get the fitted risks directly from the model

  if(link == "risks_provided"){
    linkinvfun <- NULL
  }else{
    linkinvfun <- binomial(link = link)$linkinv
    fitted_risk_t0 <- get.risk.t0(coef, marker, linkinvfun)
    fitted_risk_t1 <- get.risk.t1(coef, marker, linkinvfun)
    
  }
  

  trt.effect <- fitted_risk_t0 - fitted_risk_t1
  if(length(unique(marker))==2){ 
    #find which value of the marker is marker negative

    Meantrteff.y1 <- mean(trt.effect[marker == unique(marker)[1]])
    Meantrteff.y2 <- mean(trt.effect[marker == unique(marker)[2]])
  
   # if(Meantrteff.y1 < d & Meantrteff.y2 <d){ stop()}
    if(Meantrteff.y1 > Meantrteff.y2){ 
      marker.neg <- as.numeric(marker==unique(marker)[2])
      model.fit$disc.marker.neg = unique(marker)[2]
      }else{
      marker.neg <- as.numeric(marker==unique(marker)[1])
      model.fit$disc.marker.neg = unique(marker)[1]
    }
  }else{
    marker.neg <- ifelse( trt.effect < d, 1, 0) # indicator of being marker negative
  }

  
  ## if we dont use marker; we use fitted risks
  if(is.null(marker))
  {
    if(default.trt =="trt all"){
      derived.data <- data.frame( trt = trt, event = event, 
                                  fittedrisk.t0 = fitted_risk_t0, 
                                  fittedrisk.t1 = fitted_risk_t1,
                                  trt.effect = trt.effect,
                                  marker.neg = marker.neg)
      
    }else{
      marker.pos <- 1-marker.neg # indicator of being marker negative
      derived.data <- data.frame( trt = trt, event = event,  
                                  fittedrisk.t0 = fitted_risk_t0, 
                                  fittedrisk.t1 = fitted_risk_t1,
                                  trt.effect = trt.effect,
                                  marker.pos = marker.pos)
      
    }
    
  }else # if we do use marker
  {
     if(default.trt =="trt all"){
       derived.data <- data.frame( trt = trt, event = event, marker = marker, 
                                     fittedrisk.t0 = fitted_risk_t0, 
                                     fittedrisk.t1 = fitted_risk_t1,
                                     trt.effect = trt.effect,
                                     marker.neg = marker.neg)

     }else{
       marker.pos <- 1-marker.neg # indicator of being marker negative
       derived.data <- data.frame( trt = trt, event = event, marker = marker, 
                                fittedrisk.t0 = fitted_risk_t0, 
                                fittedrisk.t1 = fitted_risk_t1,
                                trt.effect = trt.effect,
                                marker.pos = marker.pos)
    
     }
  }
  out <- list(derived.data=derived.data, model.fit = model.fit, functions  = functions)
  class(out) = "trtsel"
  
  out

}
