##### Internal Functions  - Not for user interaction #####

# scale data
##' @importFrom stats model.matrix
model_frame_scale <- function(formula, data, scale = TRUE){
  x <- model.matrix(formula, data)[ ,-1, drop = FALSE]
  if(scale) x  <- scale(x)
  y <- data[deparse(formula[[2]])]
  if(scale) y <- scale(y)
  return(data.frame(y, x))
}

# Generates Prior Variance given favorites
.priorVar <- function(prior_R2, ols, favorites = NULL, R2_favorites = NULL){
  k <- ols$rank - 1
  priors <- rep(prior_R2/(k), k)
  
  if(!is.null(favorites)){
    fav <- which(names(ols$coefficients[-1] )%in% favorites)
    if(length(fav) == 0) stop(paste("No variable matches", favorites))
    if(length(fav) == length(names(ols$coefficients[-1]))) stop("You can't choose all variables as favorites")
    priors[fav]  <-  (R2_favorites)/(length(fav))
    priors[-fav] <-  (prior_R2 - R2_favorites)/(length(priors[-fav]))
  }
  
  V <- diag(k)*(priors)
  
  return(V)
}


# calculates s-values given a vector
.s_value <- function(vector){
  min <- min(vector)
  max <- max(vector)
  (max +  min) / (max -  min)
}



##' @import caTools
r2_combs <- function(R2_bounds){
  r2_combs  <- as.data.frame(t(combs(R2_bounds, 2)))
  ord       <- order(t(r2_combs)[,1, drop = FALSE], -order(t(r2_combs)[,2, drop = FALSE]))
  r2_combs  <- r2_combs[,ord, drop = FALSE]
  names(r2_combs) <- paste("R2", lapply(r2_combs, paste, collapse = "_"), sep = "_")
  r2_combs
}


# Calculates the bayesian estimates given a prior R2
##' @importFrom stats vcov
.bayes <- function(ols, prior_R2, favorites = NULL, R2_favorites = NULL){
  
  if(class(ols)!= "lm") stop("ols must be an object of class 'lm'")
  
  #--- prior variance
  V <- .priorVar(prior_R2, ols, favorites, R2_favorites)
  #--- ols coefficients
  
  b <- ols$coefficients[-1]
  
  #--- ols precision
  H <- solve(vcov(ols)); H <- H[-1, , drop = FALSE]; H <- H[,-1, drop = FALSE]
  
  #-- bayes estimate
  bayes <- solve((H+solve(V)))%*%H%*%b
  
  colnames(bayes) <- paste("b_bayes", prior_R2, sep="_")
  
  #-- bayes var-covar
  bayes_cov <- solve((H+solve(V)))
  
  #-- bayes t
  t <- bayes/sqrt(diag(bayes_cov))
  colnames(t) <- paste("t_bayes", prior_R2, sep="_")

  #-- formatting data before returning
  coefs <- data.frame(coef=bayes)
  names(coefs) <- paste(c("b"), prior_R2, sep="_")
  
  ret <- list(coefficients=bayes, vcov=bayes_cov,  t=t)
  return(ret)
}




# Extreme bounds given R2 bounds
##' @importFrom stats vcov
.bounds <- function(ols, coef_indices, R2_bounds, favorites = NULL, R2_favorites = NULL){
  
  #-- linear combination --#
  k <- ols$rank - 1
  phi <- rep(0, k)
  phi[coef_indices] <- 1
  
  #-- prior variances --#
  V1 <- .priorVar(min(R2_bounds), ols, favorites, min(R2_favorites))
  V2 <- .priorVar(max(R2_bounds), ols, favorites, max(R2_favorites))
  
  #--- Extreme bounds --#
  b <- ols$coefficients[-1]
  H <- solve(vcov(ols)); H <- H[-1, , drop = FALSE]; H <- H[,-1, drop = FALSE]
  G <- (H + solve(V2))%*%solve((solve(V1)-solve(V2)))%*%(H+solve(V2)) + (H+solve(V2))
  f <- solve(H+solve(V1))%*%(H%*%b + (solve(V1)-solve(V2))%*%solve(H+solve(V2))%*%(H%*%b)/2)
  c <- t(b)%*%H%*%solve(H + solve(V2))%*%(solve(V1)-solve(V2))%*%solve(H+solve(V1))%*%(H%*%b)/4
  up <- t(phi)%*%f + sqrt(t(phi)%*%solve(G)%*%phi)%*%sqrt(c)
  low <- t(phi)%*%f - sqrt(t(phi)%*%solve(G)%*%phi)%*%sqrt(c)
  data.frame(low, up)
}


# helper function to melt coefficients for plotting
.melt.coef <- function(x, type){
  coef <- coef(x, type)
  coef$coefficients <- rownames(coef)
  melted <- suppressMessages(melt(coef, stringsAsFactors=FALSE, value.name = type))
  return(melted)
}


# breaks text for pretty printing
text_break <- function(text, print.length = 5){
  if(print.length < 1) print.length <- 1
  if(length(text) > print.length){
    more <- paste("and", length(text) - print.length, "more.")
    text <- c(text[1:print.length], more)
  }
  text
}