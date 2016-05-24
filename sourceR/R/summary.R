# Works out a 95% HPD region based on Chen-Shao 1999

summary.source_attribution <- function(object, alpha = 0.05, burn_in = 0, thin = 1, interval_type = "chen-shao", ...) {
  
#   chen_shao_hpd <- function(x, alpha) {
#     
#     # x is the MCMC output
#     
#     n <- length(x)
#     sorted <- sort(x)
#     
#     upper_pos <- round(n * (1 - alpha))
#     
#     ci.lower <- 0
#     ci.upper <- 0
#     region_size <- Inf
#     
#     for (i in 1 : (n - upper_pos)) {
#       test_interval <- sorted[upper_pos + i] - sorted[i]
#       if (test_interval < region_size) {
#         region_size <- test_interval
#         ci.lower <- sorted[i]
#         ci.upper <- sorted[upper_pos + i]
#       }
#     }
#     
#     return(c(median = median(sorted), lower = ci.lower, upper = ci.upper))
#   }
  
#   SPIn_hpd <- function(x, alpha) {
#     region <- SPIn(x, conf = 1 - alpha)
#     return(c(median = median(x), lower = region$spin[1], upper = region$spin[2]))
#   }

  create_calc_interval <- function(interval_type) { 
    if (interval_type == "chen-shao") {
    return(function(x, alpha) {
      # x is the MCMC output
      
      n <- length(x)
      sorted <- sort(x)
      
      upper_pos <- round(n * (1 - alpha))
 
      if (upper_pos == n) {
        stop("The number of iterations must be larger, or alpha must be larger.")
      }
      
      ci.lower <- 0
      ci.upper <- 0
      region_size <- Inf
      
      for (i in 1 : (n - upper_pos)) {
        test_interval <- sorted[upper_pos + i] - sorted[i]
        if (test_interval < region_size) {
          region_size <- test_interval
          ci.lower <- sorted[i]
          ci.upper <- sorted[upper_pos + i]
        }
      }
      
      return(c(median = median(sorted), lower = ci.lower, upper = ci.upper))
    })
    } else if (interval_type == "spin") {
      return(function(x, alpha) {
        region <- tryCatch(
          {
            SPIn(x, conf = 1 - alpha)$spin
          },
          error = function(cond) {
            print("Error calculating SPIn interval.")
            return(c(NA, NA))}
          )
        return(c(median = median(x), lower = region[1], upper = region[2]))
      })
    } else {
      stop("interval_type must be \"spin\" or \"chen-shao\"")
    }
  }
  calc_interval <- create_calc_interval(interval_type)
  
  summary_res <- list()
  n <- dim(object$posterior$a[[1]][[1]])[1]

  sub <- seq(burn_in + 1, n, by = thin)
  
  sources <- dim(object$posterior$a[[1]][[1]])[2]
  times <- length(object$posterior$a)
  locations <- length(object$posterior$a[[1]])
  types <- length(object$data_nested$humans[[1]][[1]])
  
  if ("q" %in% names(object$posterior)) {
    summary_res$q <- matrix(NA, ncol = 3, nrow = dim(object$posterior$q)[2])
    colnames(summary_res$q) <- c("median", "lower", "upper")
    rownames(summary_res$q) <- colnames(object$posterior$q)
    for (i in 1 : dim(object$posterior$q)[2]) {
      summary_res$q[i, ] <- calc_interval(object$posterior$q[sub, i], alpha)
    }
  }

  if ("alpha" %in% names(object$posterior)) {
    summary_res$alpha <- calc_interval(object$posterior$alpha[sub], alpha)
    names(summary_res$alpha) <- c("median", "lower", "upper")
  }
  
  if ("d" %in% names(object$posterior)) {
    summary_res$d <- calc_interval(object$posterior$d[sub], alpha)
    names(summary_res$d) <- c("median", "lower", "upper")
  }
  
  summary_res$a <- list()
  for (t in 1 : times) {
    summary_res$a[[t]] <- list()
    names(summary_res$a)[t] <- names(object$posterior$a)[t]
    for (l in 1 : locations) {
      summary_res$a[[t]][[l]] <- matrix(NA, ncol = 3, nrow = sources)
      names(summary_res$a[[t]])[l] <- names(object$posterior$a[[t]])[l]
      colnames(summary_res$a[[t]][[l]]) <- c("median", "lower", "upper")
      rownames(summary_res$a[[t]][[l]]) <- colnames(object$posterior$a[[t]][[l]])
      for (j in 1 : sources) {
        summary_res$a[[t]][[l]][j, ] <- calc_interval(object$posterior$a[[t]][[l]][sub, j], alpha)
      }
    }
  }
  
  if ("lj" %in% names(object$posterior)) {
    summary_res$lj <- list()
    for (t in 1 : times) {
      summary_res$lj[[t]] <- list()
      names(summary_res$lj)[t] <- names(object$posterior$lj)[t]
      for (l in 1 : locations) {
        summary_res$lj[[t]][[l]] <- matrix(NA, ncol = 3, nrow = sources)
        names(summary_res$lj[[t]])[l] <- names(object$posterior$lj[[t]])[l]
        colnames(summary_res$lj[[t]][[l]]) <- c("median", "lower", "upper")
        rownames(summary_res$lj[[t]][[l]]) <- colnames(object$posterior$lj[[t]][[l]])
        for (j in 1 : sources) {
          summary_res$lj[[t]][[l]][j, ] <- calc_interval(object$posterior$lj[[t]][[l]][sub, j], alpha)
        }
      }
    }
    
    lj_proportion <- list()
    summary_res$lj_proportion <- list()
    for (t in 1 : times) {
      summary_res$lj_proportion[[t]] <- list()
      lj_proportion[[t]] <- list()
      names(summary_res$lj_proportion)[t] <- names(object$posterior$lj)[t]
      for (l in 1 : locations) {
        lj_proportion[[t]][[l]] <- t(apply(X=object$posterior$lj[[t]][[l]], 1, function(x) x / sum(x)))
        summary_res$lj_proportion[[t]][[l]] <- matrix(NA, ncol = 3, nrow = sources)
        names(summary_res$lj_proportion[[t]])[l] <- names(object$posterior$lj[[t]])[l]
        colnames(summary_res$lj_proportion[[t]][[l]]) <- c("median", "lower", "upper")
        rownames(summary_res$lj_proportion[[t]][[l]]) <- colnames(object$posterior$lj[[t]][[l]])
        for (j in 1 : sources) {
          summary_res$lj_proportion[[t]][[l]][j, ] <- calc_interval(lj_proportion[[t]][[l]][sub, j], alpha)
        }
      }
    }
  }

  if ("li" %in% names(object$posterior)) {
    summary_res$li <- list()
    for (t in 1 : times) {
      summary_res$li[[t]] <- list()
      names(summary_res$li)[t] <- names(object$posterior$li)[t]
      for (l in 1 : locations) {
        summary_res$li[[t]][[l]] <- matrix(NA, ncol = 3, nrow = types)
        names(summary_res$li[[t]])[l] <- names(object$posterior$li[[t]])[l]
        colnames(summary_res$li[[t]][[l]]) <- c("median", "lower", "upper")
        rownames(summary_res$li[[t]][[l]]) <- colnames(object$posterior$li[[t]][[l]])
        for (i in 1 : types) {
          summary_res$li[[t]][[l]][i, ] <- calc_interval(object$posterior$li[[t]][[l]][sub, i], alpha)
        }
      }
    }
  }
  
  if ("r" %in% names(object$posterior)) {
    summary_res$r <- list()
    for (t in 1 : times) {
      summary_res$r[[t]] <- array(NA,  dim = (c(types, 3, sources)), dimnames = list(rownames(object$posterior$r[[t]]), c("median", "lower", "upper"), colnames(object$posterior$r[[t]])))
      names(summary_res$r)[t] <- names(object$posterior$r)[t]
      for (j in 1 : sources) {
        for (i in 1 : types) {
          summary_res$r[[t]][i, , j] <- calc_interval(object$posterior$r[[t]][i, j, sub], alpha)
        }
      }
    }
  }

  return(summary_res)
}
