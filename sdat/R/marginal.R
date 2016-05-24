marginal.test <- function(x, y, n_sim = 5000L)
{
  time0 <- proc.time()[3L]
  
  if (!is.matrix(x)) x <- cbind(x)
  n <- nrow(x)
  p <- ncol(x)
  
  alpha <- c(1 : p)
  n_alpha <- length(alpha)
  
  result <- .C("R_detect_effect", 
    x = as.double(x), y = as.double(y),
    alpha = as.double(alpha), extreme = double(n_alpha + 1L), 
    n = as.integer(n), p = as.integer(p), 
    n_alpha = as.integer(n_alpha), n_sim = as.integer(n_sim))
  
  time0 <- as.double(proc.time()[3L] - time0)
  c("Max" = result$extreme[1L], 
    "Sum" = result$extreme[n_alpha], 
    "Adaptive" = result$extreme[n_alpha + 1L],
    "Time" = time0)
}

