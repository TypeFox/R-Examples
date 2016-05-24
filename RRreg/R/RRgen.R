# simulation of randomization procedure:
# input: model parameters
# output: data frame with randomized response
#' @export
#' @title Generate randomized response data
#' @description The method \code{RRgen} generates data according to a specified RR model, 
#' e.g., \code{"Warner"}. True states are either provided by a vector \code{trueState} or drawn randomly from a Bernoulli distribution. Useful for simulation and testing purposes, e.g., power analysis.
#' @param n sample size of generated data
#' @param pi.true true proportion in population (a vector for m-categorical \code{"FR"} or \code{"custom"} model)
#' @param model specifes the RR model, one of: \code{"Warner"}, \code{"UQTknown"}, \code{"UQTunknown"}, \code{"Mangat"}, \code{"Kuk"}, \code{"FR"}, \code{"Crosswise"}, \code{"CDM"}, \code{"CDMsym"}, \code{"SLD"},  \code{"mix.norm"},  \code{"mix.exp"}, \code{"custom"}. See \code{vignette("RRreg")} for details. 
#' @param p randomization probability (depending on model, see \code{\link{RRuni}} for details)
#' @param complyRates vector with two values giving the proportions of carriers and non-carriers who adhere to the instructions, respectively
#' @param sysBias probability of responding 'yes' (coded as 1) in case of non-compliance for carriers and non-carriers of the sensitive attribute, respectively. If \code{sysBias=c(0,0)}, carriers and non-carriers systematically give the nonsensitive response 'no' (also known as self-protective(SP)-'no' responses). If \code{sysBias=c(0,0.5)}, carriers always respond 'no' whereas non-carriers randomly select a response category. Note that \code{sysBias = c(0.5,0.5)} might be the best choice for \code{Kuk} and \code{Crosswise}. For the m-categorical \code{"FR"} or \code{"custom"} model, \code{sysBias} can be given as a probability vector for categories 0 to (m-1).
#' @param groupRatio proportion of participants in group 1. Only required for two-group models, e.g., \code{SLD} and \code{CDM}
#' @param Kukrep Number of repetitions of Kuk's procedure (how often red and black cards are drawn)
#' @param trueState optional vector containing true states of participants (i.e., 1 for carriers and 0 for noncarriers of sensitive attribute; for \code{FR}: values from 0,1,...,M-1 (M = number of response categories) which will be randomized according to the defined procedure (if specified, \code{n} and \code{pi.true} are ignored)
#' @return \code{data.frame} including the variables \code{true} and \code{response} (and for \code{SLD} and \code{CDM} a third variable \code{group}) 
#'@details If \code{trueState} is specified, the randomized response procedure will be simulated for this vector, otherwise a random vector of length \code{n} with true proportion \code{pi.true} is drawn. Respondents answer biases can be simulated by adjusting the compliance rates: if \code{complyRates} is set to \code{c(1,1)}, all respondents adhere to the randomization procedure. If one or both rates are smaller than 1, \code{sysBias} determines whether noncompliant respondents systematically choose the nonsensitive category or whether they answer randomly. 
#'
#'\code{SLD} - to generate data according to the stochastic lie detector with the proportion \code{t} of honest carriers, parameters are set to \code{complyRates=c(t,1)} and \code{sysBias=c(0,0)}
#'
#'\code{CDM} - to generate data according to the cheating detection model with the proportion \code{gamma} of cheaters, parameters are set to \code{complyRates=c(1-gamma,1-gamma)} and \code{sysBias=c(0,0)}
#'         
#' @seealso see \code{vignette('RRreg')} for a detailed description of the models and \code{\link{RRlog}}, \code{\link{RRlin}} and \code{\link{RRcor}} for the multivariate analysis of RR data
#' @examples 
#' # Generate responses of 1000 people according to Warner's model,
#' # every participant complies to the RR procedure
#' genData <- RRgen(n=1000, pi.true=.3, model="Warner", p=.7)
#' colMeans(genData)
#' 
#' # use Kuk's model with two decks of cards, 
#' # p gives the proportions of red cards for carriers/noncarriers
#' genData <- RRgen(n=1000, pi.true=.4, model="Kuk", p=c(.4,.7))
#' colMeans(genData)
#' 
#' # Stochastic Lie Detector (SLD):
#' # Only 80% of carriers answer according to the RR procedure
#' genData <- RRgen(n=1000, pi.true=.2, model="SLD", p=c(.2,.8),
#'                  complyRates=c(.8,1),sysBias=c(0,0))
#' colMeans(genData)

RRgen <- function(n,pi.true, model, p, 
                  complyRates = c(1,1), sysBias = c(0,0),
                  groupRatio=.5, Kukrep=1,trueState=NULL){
  model <- match.arg(model, c("Warner","UQTknown","UQTunknown","Mangat",
                              "Kuk","FR","Crosswise","CDM","CDMsym","SLD",
                              "mix.norm", "mix.exp", "custom"))
  true <- NULL
  if(!is.null(trueState)){
    trueState <- as.numeric(trueState)
    RRcheck.xp(model,trueState,p,"trueState")
    n <- length(trueState)
    true <- trueState
    pi.true <- mean(trueState)
    if (model %in% c("FR", "custom") && nrow(as.matrix(p))>2){
      pi.true <- table(true)/length(true)
    }
    if (model %in% c("mix.norm")){
      pi.true <- c(mean(trueState), sd(trueState))
    }
    if (model %in% c("mix.exp")){
      pi.true <- c(mean(trueState))
    }
  }
  if (model %in% c("custom", "FR") && length(pi.true)==1) 
    pi.true <- c(1-pi.true,pi.true) 
  
  # check input  
  RRcheck.p(model,p)
  RRcheck.pi(model,pi.true,n)
  RRcheck.rate(complyRates[1])
  RRcheck.rate(complyRates[2])
  RRcheck.groupRatio(groupRatio)
  if ( any(complyRates != 1 ) ){
    if (model %in% c("custom", "FR") &&
        (length(sysBias)!=nrow(as.matrix(p)) || sum(sysBias)!=1 || any(sysBias<0) || any(sysBias>1))  ){
      warning("For the m-categorical FR/custom model, the argument 'sysBias' must be a probability vector 
              of the same length as 'p'. sysBias is set to equal guessing across categories automatically.")
      sysBias <- rep(1,length(pi.true))/length(pi.true)
    }else if(!(model %in% c("custom", "FR")) && (min(sysBias)<0 || max(sysBias)>1 || length(sysBias) !=2))
        stop("The argument 'sysBias' gives the probabilities of 'no'-responses in case of non-compliance 
           for carriers and non-carriers, respectively (e.g., sysBias = c(0, 0.5).")
  }
  # initialisiere
  response <- rep(0,n)
  comply <- rep(1,n)
  
  randNum <- runif(n) #random numbers for data generation
  
  
  ##################
  # continuous mixture RR models
  if (model %in% c("mix.norm","mix.exp")){
    comply <- rep(1,n) # always comply
    if (model =="mix.norm"){
      if (is.null(true)) 
        true <- rnorm(n, mean=pi.true[1], sd=pi.true[2])
      mask <- rnorm(n, mean=p[2], sd=p[3])
    }else if (model =="mix.exp"){
      if (is.null(true)) 
        true <- rexp(n, rate=1/pi.true[1])
      mask <- rexp(n, rate=1/p[2])
    }    
    response <- ifelse(randNum<p[1], true, mask)  
  }
  # Forced response: multinomial response categories possible
  else if (model %in% c("FR", "custom")){ # && length(pi.true>1)){
    numCat <- nrow(as.matrix(p))
    if (numCat != length(pi.true)){
      stop("The length of vector 'pi.true' and 'p' has to match in the 'FR'/'custom' model")
    }
    # distribute values across categories acording to 'pi.true'
    if (is.null(trueState)){
      true <- findInterval(runif(n),cumsum(pi.true))
    } 
    
    if(model == "FR"){
      # return response 'i' with probability 'p[i]', otherwise true response x
      chooseCat <- findInterval(runif(n),cumsum(p))
      responseComply <- ifelse(chooseCat==length(p),
                               true,
                               chooseCat) 
    }else if (model == "custom"){
      pcumsum <- apply(p[,true+1], 2, cumsum)
      responseComply <- apply( rbind(pcumsum, runif(n)), 2, function(xx) findInterval(xx[numCat+1], xx[1:numCat]))
    }
    
    if(missing(complyRates) || all(complyRates ==1))
      complyRates <- rep(1, numCat)
    if ( any(complyRates!=1) && length(complyRates) != numCat){
      warning("For the polytomous forced response ('FR'/'custom') model,'complyRates' 
                must have the same length as 'p', defining the compliance rate for 
                each of the true states separately.")
      
    }
    responseNonComply <- findInterval(runif(n),cumsum(sysBias))
    
    # do participants follow the instructions:
    comply <- rep(NA,n)
    for (i in 1:numCat){
      comply[true == i-1] <- ifelse(randNum[true == i-1]<complyRates[i] ,1,0)
    }
    response <- ifelse(comply == 1 , 
                       responseComply ,     # participants comply
                       responseNonComply)   # participants don't comply
    
    
    ##############################
    # for dichotomous models
  }else{
    if (model %in% c("CDM","CDMsym")){
      # adjustment, so pi.true will fit to the estimation
      pi.true <- 2*pi.true / (complyRates[1]+complyRates[2])
      if(pi.true>1) stop("For CDM and CDMsym, change arguments complyRates and/or pi in order to generate valid data.")
    }
    # sensitive attribute is binomiallly distributed with probabiliy pi.true
    if (is.null(trueState)){
      true <- rbinom(n,1,pi.true[1])
    }
    
    # unbiased response according to instructions:
    switch(model,
           "FR" = {
             chooseCat <- findInterval(runif(n),cumsum(p));
             responseComply <- ifelse(chooseCat==length(p),
                                      true,      # true response
                                      chooseCat) # forced response
           },
           "Warner" = {
             responseComply <- ifelse(randNum<p,
                                      true,  # normal question
                                      1-true)  # reversed question
           },
           "Mangat" = {
             responseComply <- ifelse(true==1,
                                      1,          # carriers answer honestly
                                      ifelse(randNum<p,
                                             0,   # true answer with prob 'p'
                                             1))  # noncarriers forced to answer 1 with prob '1-p'
           },
           "Kuk" = {          # p[1], p[2] give proportion of red cards
             responseComply <- ifelse(true==1, 
                                      rbinom(n,Kukrep,p[1]),    #card deck for carriers
                                      rbinom(n,Kukrep,p[2]) )  # for noncarriers
           },
           "UQTknown" = {
             responseComply <- ifelse(randNum<p[1], 
                                      true, # answer to relevant question with probability p[1]
                                      ifelse(runif(n)<p[2],1,0)) 
             # answer to irrelevant question "Yes" with probability p[2]
           },
           "UQTunknown" = {
             if (length(pi.true)==1){
               pi.true <- c(pi.true,runif(1,.2,.8))
               #                warning(paste0("The prevalence of the unrelated question was randomly set to ", round(pi.true[2],3),". To explicitly choose the prevalence, set 'pi.true=c(pi.sensitive, pi.unrelated)'"),call.=F)
             }
             randUQ <- runif(n)
             nn1 <- round(n*groupRatio)
             split <- c(rep(1,nn1),rep(2,n-nn1))
             responseComply <- ifelse(split==1, 
                                      ifelse( randNum<p[1],     # group 1: answer sens question with prob p[1]
                                              true,  
                                              ifelse (randUQ<pi.true[2],1,0)),
                                      ifelse( randNum<p[2],     # group 2: answer sens question with prob p[2]
                                              true,  
                                              ifelse (randUQ<pi.true[2],1,0)) )
           },
           "Crosswise"={
             responseComply <- ifelse(randNum<p,    # see Warner
                                      true,
                                      1-true)             
           },
           "CDM" = {
             nn1 <- round(n*groupRatio)
             split <- c(rep(1,nn1),rep(2,n-nn1))
             responseComply <- ifelse(true==1, 
                                      1,   # honest-carriers always answer "yes"
                                      ifelse( split ==1,     # honest non-carriers
                                              ifelse (randNum<p[1],1,0),  # group 1
                                              ifelse (randNum<p[2],1,0)   # group 2
                                      ))
           },
           "CDMsym" = {
             nn1 <- round(n*groupRatio)
             split <- c(rep(1,nn1),rep(2,n-nn1))
             responseComply <- ifelse(true==1, 
                                      ifelse( split == 1,     # honest-carriers 
                                              ifelse(randNum<p[2],0,1),   # forced No  group1
                                              ifelse(randNum<p[4],0,1)),  # forced No  group2
                                      ifelse( split == 1,     # honest non-carriers
                                              ifelse (randNum<p[1],1,0),  # forced Yes   group1
                                              ifelse (randNum<p[3],1,0)   # forced Yes   group2
                                      ))
           },
           "SLD" = {
             nn1 <- round(n*groupRatio)
             split <- c(rep(1,nn1),rep(2,n-nn1))
             responseComply <- ifelse(true==1,
                                      1, #carriers are not randomized
                                      ifelse( split ==1,     # honest non-carriers
                                              ifelse (randNum<p[1],0,1),  
                                              ifelse (randNum<p[2],0,1)
                                      )) 
           })
    
    # biased answer (nonComply)
    #     if (model %in% c("Kuk" ,"Crosswise") && sysBias!=0.5){
    #       sysBias <- 0.5
    #       #warning("Parameter 'sysBias' is ignored since it is not meaningful for Kuk's or Crosswise model.")
    #     }
    if (model =="Kuk"){
      responseNonComplyC <- rbinom(n,Kukrep,sysBias[1])    # Kuk: always random answer in case of  non-compliance
      responseNonComplyNC <- rbinom(n,Kukrep,sysBias[2])
    } else{
      # cheaters always give the answer "no" (sysBias=0)
      responseNonComplyC <- rbinom(n,1,sysBias[1])   # random answers
      responseNonComplyNC <- rbinom(n,1,sysBias[2]) 
    }
    randNum <- runif(n) #random numbers for compliance
    # do participants follow the instructions:
    response[true==1] <- ifelse(randNum<complyRates[1] , 
                                responseComply ,              # carriers comply
                                responseNonComplyC)[true==1]   # carriers don't comply
    response[true==0] <- ifelse(randNum<complyRates[2] , 
                                responseComply ,              # noncarriers comply
                                responseNonComplyNC)[true==0]   # noncarriers don't comply
    comply[true==1] <- ifelse(randNum<complyRates[1] ,1,0)[true==1]
    comply[true==0] <- ifelse(randNum<complyRates[2] ,1,0)[true==0] 
  }
  
  
  # construct data frame
  data <- data.frame(true,comply,response)
  if (model %in% c("UQTunknown","SLD","CDM","CDMsym") ){
    data$group <- split
  }
  return(data)
}
