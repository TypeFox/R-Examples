stratasamp <-
function(n, Nh, Sh = NULL, Ch = NULL, type = 'prop')
{
### input handling
  type <- match.arg(type,c('prop', 'opt', 'costopt'))  

  # strata size
  if(!(is.numeric(Nh) && is.vector(Nh))){
    stop("Invalid input: ", sQuote("Nh"), " has to be a numeric vector.")
  }
  if(n < length(Nh)){
    stop("Invalid input: the sample size ", sQuote("n"), " has to be larger than number of strata.")
  }
  if(any(Nh < 1)){
    stop("Invalid input: population in ", sQuote("Nh"), " has to be positiv.")
  }
  # number of strata
  M <- length(Nh)
  # population size
  N <- sum(Nh)
 
  # strata deviation 
  if(type=='opt' || type=='costopt'){
    if (!(is.numeric(Sh) && is.vector(Sh))){
      stop("Invalid input: ", sQuote("Sh"), " has to be a numeric vector.")
    }
    if( length(Sh) != M){
      stop("Invalid input: ", sQuote("Sh"), " and ", sQuote("Nh"), " must have the same length.")
    }
    if(any(Sh < 0)){
      stop("Invalid input: standard diviation in ", sQuote("Sh"), " can not be negativ.")
    }
  }
  # strata costs
  if(type=='costopt'){
    if (!(is.numeric(Ch) && is.vector(Ch))){
      stop("Invalid input: ", sQuote("Ch"), " has to be a numeric vector.")
    }
    if ( length(Ch) != M ){
      stop("Invalid input: ", sQuote("Ch"), " and ", sQuote("Nh"), " must have the same length.")
    }
    if (any(Ch < 1)){
      stop("Invalid input: cost in ", sQuote("Ch"), " can not be negativ.")
    }
  }

### compute sample weights for each strata
  wh <- nh  <- integer(M)
  # proportional
  if(type=='prop'){
    wh <- Nh/N
  }
  # optimal
  if (type == 'opt'){
    for(i in 1:M){
      wh[i] <- (Nh[i] * Sh[i]) / sum(Nh*Sh)
    }
  } 
  # costoptimal 
  if(type == 'costopt'){
    for (i in 1:M){
      wh[i] <- (Nh[i] * Sh[i] / sqrt(Ch[i])) / sum(Nh*Sh/sqrt(Ch))
    }
  }

### compute sample size for each strata
  for (i in 1:M){
    nh[i] <- round(n * wh[i])
  }
  if(any(nh < 2)){
    warning("Warning: Sampling of less than 2 observations in a stratum is not recommended!") 
  }

### result object
  res <- rbind(1:length(Nh), nh)
  rownames(res) <- c("Stratum","Size")
  colnames(res) <- c(rep("",length(Nh)))
  return(res)
}
