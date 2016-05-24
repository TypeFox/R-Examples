#' d-Jackknife estimation.
#' @param x vector with simulations or matrix with variables in columns, simulations in rows.
#' @param fun function to apply
#' @param d how many elements to leave out. For "leave one out",  d= 1 (default)
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @keywords internal


int.jackknife <- function(x, fun, d = 1, alpha = 0.05) {  
                                     
  x <- data.frame(x)
  nloci <- nrow(x)                           #INDIVIDUALS OR POPULATIONS IN COLUMNS, 
  dcomp = nloci - d  
  N <- choose(nloci, dcomp)                      # TRAIT TO LEAVE OUT IN ROWS
  crit <- abs(qt(1-(alpha/2), N  -1))
  combinations <- combn(nloci, dcomp)
  
    obs <- apply(x, 2, fun)
    obs2 <- matrix(obs, nrow= N, ncol= length(obs), byrow=T)
    jack <- x - x
    for(i in 1:N) {
      icomb <- combinations[, i]
      jack[i, ] <- apply(x[icomb, , drop = FALSE], 2, fun)
    }

    pseudo <- (N * obs2) - ((N - 1) * jack)
    theta <- apply(pseudo, 2, mean)
    bias <- (N - 1) * (obs - theta)
    pseudo.variance <- apply(pseudo, 2, var)		
    sd.jack <- sqrt(pseudo.variance / N)
    interval <- sd.jack * crit
    CI<- rbind(theta - interval, theta + interval)
    colnames(CI) <- colnames(x)
    rownames(CI) <- c("lwr", "uppr")
  
  result <- list(obs = obs, 
                 theta = theta,
                 sd = sd.jack,
                 bias = bias,
                 CI = CI,
                 alpha = alpha,
                 t.crit = crit,
                 d = d,
                 comb = N)
  
  result
}
