###  Copyright (C) 2011-2012 Berwin A. Turlach <Berwin.Turlach@gmail.com>
###
###  This program is free software; you can redistribute it and/or modify
###  it under the terms of the GNU General Public License as published by
###  the Free Software Foundation; either version 2 of the License, or
###  (at your option) any later version.
###
###  This program is distributed in the hope that it will be useful,
###  but WITHOUT ANY WARRANTY; without even the implied warranty of
###  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
###  GNU General Public License for more details.
###
###  You should have received a copy of the GNU General Public License
###  along with this program; if not, write to the Free Software
###  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
###  USA.
###

evalPol <- function(x, beta){
  res <- 0
  for(bi in rev(beta))
    res <- res*x + bi
  res
}

curvPol <- function(x, beta){
  b1der <- beta[-1]
  b1der <- b1der*(1:length(b1der))
  b2der <- b1der[-1]
  b2der <- b2der*(1:length(b2der))
  res1 <- 0
  for(bi in rev(b1der))
    res1 <- res1*x + bi
  res2 <- 0
  for(bi in rev(b2der))
    res2 <- res2*x + bi
  res2/(1+res1^2)^(1.5)
}

coef.monpol <- function(object,
                        scale=c("original", "fitted"),
                        type=c("beta", "monpar"), ...){

  scale <- match.arg(scale)
  type <- match.arg(type)

  if(scale == "original"){
    if(type == "beta"){
      object$beta.raw
    }else{
      stop("scale='original' and type='monpar' not yet implemented.")
    }
  }else{
    if(type == "beta"){
      object$beta
    }else{
      if(object$algorithm %in% c("Hawkins", "NoIntercept"))
        stop("Algorithms 'Hawkins' and 'NoIntercept' do not provide monotone parameterisation.")
      object$par
    }
  }
}

fitted.monpol <- function(object,
                          scale=c("original", "fitted"), ...){

  scale <- match.arg(scale)
  val <- napredict(object$na.action, object$fitted.value)
  if(scale == "original")
    val <- object$scly*(val+1)/2+object$miny
  val
}
  
residuals.monpol <- function(object,
                             scale=c("original", "fitted"), ...){

  scale <- match.arg(scale)
  val <- naresid(object$na.action, object$residuals)
  if(scale == "original")
    val <- object$scly*val/2
  val
}
  
model.matrix.monpol <- function(object,
                                scale=c("original", "fitted"), ...){

  scale <- match.arg(scale)

  if (n_match <- match("x", names(object), 0L)){
    val <- object[[n_match]]
  }else{ 
    data <- model.frame(object, ...)
    val <- model.matrix.default(object, data = data)
  }
  if(scale=="fitted")
    val <- 2*(val-object$minx)/object$sclx - 1
  val
}

formula.monpol <- function(x, ...){
  formula(x$terms)
}

print.monpol <- function(x, digits = max(3, getOption("digits") - 3), ...){
  cat("\nMonotone polynomial model\n")
  cat("Call:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  cff <- coef(x, ...)
  if (length(cff)) {
    cat("Coefficients:\n")
    print.default(format(cff, digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  else cat("No coefficients\n")
  cat("\n")
  invisible(x)
}

predict.monpol <- function(object, newdata,
                           scale=c("original", "fitted"), ...){

  scale <- match.arg(scale)
  if (missing(newdata)) 
    return(fitted(object, scale))

  newdata <- as.data.frame(newdata)
  Terms <- delete.response(object$terms)
  m <- model.frame(Terms, newdata, na.action = na.omit) 
  if (!is.null(cl <- attr(Terms, "dataClasses"))) 
    .checkMFClasses(cl, m)
  x <- model.matrix(Terms, m, contrasts.arg = object$contrasts)

  if(scale=="original"){
    evalPol(x, object$beta.raw)
  }else{
    evalPol(x, object$beta)
  }
}

##
## pretty much based on logLik.nls
##
logLik.monpol <- function (object, REML = FALSE, ...){
    if (REML) 
        stop("cannot calculate REML log-likelihood for \"monpol\" objects")
    n <- length(object$residuals)
    if(is.null(w <- object$weights))
      w <- rep.int(1, n)
    zw <- w==0
    rss.unscaled <- object$scly^2*object$RSS/4
    val<- -n*(log(2*pi) + 1 - log(n) - sum(log(w + zw)) + log(rss.unscaled))/2
    attr(val, "df") <- 1L + ifelse(object$algorithm=="NoIntercept", 
                                   length(coef(object)) - 1, 
                                   length(coef(object)))
    attr(val, "nobs") <- attr(val, "nall") <- sum(!zw)
    class(val) <- "logLik"
    val
}

nobs.monpol <- function(object, ...){
  if (is.null(w <- object$weights)) length(object$residuals) else sum(w != 0) 
}
