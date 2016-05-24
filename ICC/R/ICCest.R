ICCest <- function(x, y, data = NULL, alpha = 0.05, CI.type = c("THD", "Smith")){
  square <- function(z){z^2}
  icall <- list(y = substitute(y), x = substitute(x))

  if(is.character(icall$y)){
    warning("passing a character string to 'y' is deprecated since ICC vesion 2.3.0 and will not be supported in future versions. The argument to 'y' should either be an unquoted column name of 'data' or an object")
    if(missing(data)) stop("Supply either the unquoted name of the object containing 'y' or supply both 'data' and then 'y' as an unquoted column name to 'data'")
    icall$y <- eval(as.name(y), data, parent.frame())
  } 
  if(is.name(icall$y)) icall$y <- eval(icall$y, data, parent.frame())
  if(is.call(icall$y)) icall$y <- eval(icall$y, data, parent.frame())
  if(is.character(icall$y)) icall$y <- eval(as.name(icall$y), data, parent.frame())


  if(is.character(icall$x)){
    warning("passing a character string to 'x' is deprecated since ICC vesion 2.3.0 and will not be supported in future versions. The argument to 'x' should either be an unquoted column name of 'data' or an object")
    if(missing(data)) stop("Supply either the unquoted name of the object containing 'x' or supply both 'data' and then 'x' as an unquoted column name to 'data'")
    icall$x <- eval(as.name(x), data, parent.frame())
  } 
  if(is.name(icall$x)) icall$x <- eval(icall$x, data, parent.frame())
  if(is.call(icall$x)) icall$x <- eval(icall$x, data, parent.frame())
  if(is.character(icall$x) && length(icall$x) == 1) icall$x <- eval(as.name(icall$x), data, parent.frame())


  tdata <- data.frame(icall)
  tdata <- na.omit(tdata)
  a <- length(unique(tdata$x))

  if(!is.null(attributes(tdata)$na.action)){
     warning(cat("NAs removed from rows:\n", unclass(attributes(tdata)$na.action), "\n"))
  } 
  if(!is.factor(tdata$x)){
     warning("'x' has been coerced to a factor")
     tdata$x <- as.factor(tdata$x)
  } else{
       if(length(levels(tdata$x)) > a){
          tdata$x <- factor(as.character(tdata$x), levels = unique(tdata$x))
          warning("Missing levels of 'x' have been removed")
       } 
    } 

  tmpbb <- anova(aov(y ~ x, data = tdata))
  num.df <- tmpbb$Df[1]
  denom.df <- tmpbb$Df[2]
  MSa <- tmpbb$'Mean Sq'[1]
  MSw <- var.w <- tmpbb$'Mean Sq'[2]
  tmp.outj <- aggregate(y ~ x, data = tdata, FUN = length)$y
  k <- (1/(a-1)) * (sum(tmp.outj) - (sum(square(tmp.outj)) / sum(tmp.outj)))
  var.a <- (MSa - MSw) / k
  r <- var.a / (var.w + var.a)

  low.F <- qf(alpha/2, num.df, denom.df, lower.tail = FALSE)
  N <- nrow(tdata)
  n.bar <- N/a
  n.not <- n.bar - sum(square(tmp.outj - n.bar) / ((a - 1) * N))	
    type <- match.arg(CI.type)
      if(type == "THD"){
	up.F <- qf(alpha/2, denom.df, num.df, lower.tail = FALSE)	
	FL <- (MSa/MSw) / low.F
	FU <- (MSa/MSw) * up.F
	low.CI <- (FL - 1) / (FL + n.not - 1)
	up.CI <- (FU - 1) / (FU + n.not - 1)
       }
      if(type == "Smith"){
	z.not <- qnorm(alpha/2)
	Vr <- (2*square(1-r) / square(n.not)) * ((square((1+r*(n.not-1))) / (N-a)) + ((a-1)*(1-r)*(1+r*(2*n.not-1))+square(r)*(sum(square(tmp.outj))-2*(1/N)*sum((tmp.outj^3))+(1/square(N))*square(sum(square(tmp.outj)))))/ square(a-1))
	low.CI <- r + z.not * sqrt(Vr)
	up.CI <- r - z.not * sqrt(Vr) 
      }

list(ICC = r, LowerCI = low.CI, UpperCI = up.CI, N = a, k = k, varw = var.w, vara = var.a)
}

