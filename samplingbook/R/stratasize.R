stratasize <-
function(e, Nh, Sh, level = 0.95, type = "prop")
{ 
### error handling

  type <- match.arg(type, c("prop", "opt"))
  # precision
  if(e < 0){
    stop("Wrong input: precision ", sQuote("e"), " has to be positive number.")
  }
  # strata size
  if(!(is.numeric(Nh) && is.vector(Nh))){
    stop("Invalid input: ", sQuote("Nh"), " has to be a numeric vector.")
  }
  if(any(Nh < 1)){
    stop("Invalid input: population in ", sQuote("Nh"), " has to be positiv.") 
  }
  # strata deviation
  if(!(is.numeric(Sh) && is.vector(Sh))){
    stop("Invalid input: ", sQuote("Sh"), " has to be a numeric vector.")
  }
  if(any(Sh < 0)){
    stop("Wrong input: any number from vector of standard deviation of all strata ", sQuote("Sh"), " has to be positive.")
  }
  # level
  if (level < 0 | level > 1){
    stop("Wrong input: ", sQuote("level"), " has to be probability between 0 and 1.")
  }else{
    q <- qnorm((1 + level)/2)
  }
  
  # computation
  N <- sum(Nh)
  S <- sum(Sh)
  wh <- Nh/N
  if(type == "prop"){
    n = sum(wh^2 * N/(Nh-1) * Sh^2) / (e^2/q^2 + sum(wh^2 * Sh^2 / (Nh-1)))  
  } 
  else{
    Nhss = sum(Nh*Sh)
    n = sum(wh^2 * Nhss/(Nh-1) * Sh) / (e^2/q^2 + sum(wh^2 * Sh^2 / (Nh-1))) 
  }
  n = ceiling(n)

  # return argument
  ret <- list()
  ret$call <- list(e = e, Nh = Nh, Sh = Sh, level = level, type = type)
  ret$n <- n
  structure(ret, class = "stratasize")
}
