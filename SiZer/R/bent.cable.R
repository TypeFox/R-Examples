
# this implements the bent-cable model
# and it finds the MLE of the threshold and
# bend width by exhaustively searching the
# parameter space.  This is rather time 
# intensive.  
#
# What we should probably do is create a wrapper
# function that first tries a search algorithm
# and if it fails then resorts to the exhaustive
# search.   

bent.cable <- function(x,y, grid.size=100){
  r <- range(x);
  alpha.min <- r[1];
  alpha.max <- r[2];
  gamma.max <- (r[2] - r[1])/2;
  gamma.min <- gamma.max / grid.size;

  A <- seq(alpha.min, alpha.max, length.out=grid.size);
  G <- seq(gamma.min, gamma.max, length.out=grid.size);

  L <- matrix(nrow=grid.size, ncol=grid.size);
  SSE <- matrix(nrow=grid.size, ncol=grid.size);

  max.llikelihood <- -Inf;
  hat.alpha <- 0;
  hat.gamma <- 0;

  count.a <- 1;
  for(a in A){
    count.g <- 1;
    for(g in G){
        I.middle <- 1* (abs(x-a) < g);
        I.left   <- 1* ( x >= a+g );
        q <- (x-a+g)^2 / (4*g) * I.middle +
             (x-a)             * I.left;
        fit <- lm(y ~ x + q);
        Beta <- fit$coefficients;
        sigma <- sqrt( sum(fit$residuals^2)/length(x) );
        mu <- Beta[1] + Beta[2]*x + Beta[3]*q;

        llikelihood <- sum( log( dnorm(y, mean=mu, sd=sigma) ) );
        L[count.a, count.g] <- llikelihood;
		SSE[count.a, count.g] <- sum(fit$residuals^2);
        if(llikelihood > max.llikelihood){
  	      max.llikelihood <- llikelihood;
          hat.alpha <- a;
          hat.gamma <- g;
        }
      count.g <- count.g + 1;	 
    }
    count.a <- count.a + 1;
  }
  out <- NULL;
  out$log.likelihood <- L;
  out$SSE <- SSE;
  out$alphas <- A;
  out$gammas <- G;
  out$alpha <- hat.alpha; 
  out$gamma <- hat.gamma; 

  I.middle <- 1* (abs(x-hat.alpha) < hat.gamma);
  I.left   <- 1* ( x >= hat.alpha + hat.gamma );
  q <- (x-hat.alpha+hat.gamma)^2 / (4*hat.gamma) * I.middle +
       (x-hat.alpha)             * I.left;
  fit <- lm(y ~ x + q);
  
  out$model <- fit;
  class(out) <- 'bent_cable';
  return(out);
}
  
predict.bent_cable <- function(object, x, ...){
   alpha <- object$alpha;
   gamma <- object$gamma;
   beta <- object$model$coefficients;

  I.middle <- 1* (abs(x-alpha) < gamma);
  I.left   <- 1* ( x >= alpha + gamma );
  q <- (x - alpha + gamma)^2 / (4*gamma) * I.middle +
       (x - alpha) * I.left;

  yhat <- beta[1] + beta[2]*x + beta[3]*q;
  return(yhat);
}

logLik.bent_cable <- function(object, ...){
	out <- logLik(object$model)
	attr(out,'df') <- 6  # 3 betas, 1 gamma, 1 alpha, 1 sigma
	return(out)
}



