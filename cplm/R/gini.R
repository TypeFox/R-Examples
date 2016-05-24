#######################################################
##    The Gini index from an ordered Lorenz curve    ##
## Author: Wayne Zhang, actuary_zhang@hotmail.com    ##
#######################################################

# class defining slots common to all derived classes 
setClass("gini", 
  representation(
    call = "call",
    gini = "matrix",
    sd = "matrix",
    lorenz = "list")
)

# compute gini indices etc for model selections
gini <- function(loss, score, base = NULL, data, ...){
  call <- match.call()  
  
  if (missing(loss) || missing(score))
    stop("loss and score must be specified")
  if (!is.character(loss) || !is.character(score) ||
     (!is.null(base) && !is.character(base)))
    stop("Arguments 'loss', 'score' and 'base' must be characters!")
  # check base 
  if (length(base) > 1) {
    base <- base[1]
    warning("Only the first base premium is used!")    
  }    
  # check score 
  if (length(base) && (pos <- match(base, score, nomatch = 0))){
    score <- score[-pos]
    warning(paste(base, "is removed from score!"))
  }
  
  # get loss and standardize so that mean is 1
  y <- data[, loss] / mean(data[, loss])
  sc <- data[, score]
  
  # when the base premium is specified
  if (!is.null(base)) {
    P <- data[, base] / mean(data[, base])
    ans <- do_gini(y, P, sc)
    dimnames(ans$gini) <- dimnames(ans$sd) <- list(base, score)
    ans$lorenz <- list(ans$lorenz)
    names(ans$lorenz) <- base
  } else { # model comparison using each score as a base
    p <- length(score)
    gi <- sd <- matrix(0, p, p)
    lrz <- vector("list", p)    
    for (i in 1:p){
      P <- sc[, i] / mean(sc[, i])
      tmp <- do_gini(y, P, sc[, -i, drop = FALSE])
      gi[i, -i] <- tmp$gini
      sd[i, -i] <- tmp$sd
      lrz[[i]] <- tmp$lorenz      
    }    
    dimnames(gi) <- dimnames(sd) <- list(score, score)
    names(lrz) <- score
    ans <- list(gini = gi, sd = sd, lorenz = lrz)
  }
  out <- new("gini", call = call, gini = ans$gini, 
             sd = ans$sd, lorenz = ans$lorenz)    
  return (out)
}

# compute the gini index, std errs, and the graph of the lorenz curve
do_gini <- function(y, P, S){    
  p <- ncol(S) 
  gi <- sd <- matrix(NA, 1, p)
  ni <- 50
  lrz <- matrix(0, ni, p)
  n <- length(y)
  for (i in 1:p){
    # order by relativity
    yP <- cbind(y, P)[order(S[, i] / P), ]
    # distribution function    
    DF <- apply(yP / colSums(yP), 2, cumsum)
    DF <- rbind(c(0, 0), DF)
    # compute the gini index
    gi[1, i] <- 1 - sum(diff(DF[, 2]) * (DF[2:(n + 1), 1] + DF[1:n, 1]))    
    
    # compute the standard errors
    h1 <- 0.5 * (yP[, 2] * DF[-1, 1] + yP[, 1] * (1 - DF[-1, 2]))    
    mh1 <- mean(h1)
    mh <- 0.5 * (1 - gi[1, i])
    vh <- var(h1)  
    vhy <- cov(h1, yP[, 1])
    vhP <- cov(h1, yP[, 2])
    vy <- var(y)
    vP <- var(P)
    vyP <- cov(yP[, 1], yP[, 2])
    v <-  4 * (4 * vh + mh^2 * (vy + vP) - 4 * mh * (vhy  + vhP)  + 2 * mh^2 * vyP)      
    sd[1, i] <- sqrt(v / n) 
    
    # compute the interpolations used for plots
    lrz[, i] <- approx(DF[, 2:1], n = ni)[[2]]
  }
  #add cdf for premiums (this is the same across models)
  lrz <- cbind(seq(0, 1, length = ni), lrz) 
  dimnames(lrz)[[2]] <- c(".P.", dimnames(S)[[2]]) 
  # scale by 100
  return(list(gini = gi * 100, sd = sd *100, lorenz = lrz * 100))
}


setMethod("show", 
  signature(object = "gini"),
  function(object){
    digits <- max(3, getOption("digits") - 3)
    gi <- format(object@gini, digits = digits)
    sd <- format(object@sd, digits = digits)
    cat("\nCall:\n", paste(deparse(object@call), sep = "\n", collapse = "\n"), 
          "\n\n", sep = "")
    cat("Gini indices:\n")
    print.default(gi, print.gap = 2,  quote = FALSE)
    cat("\n")
    cat("Standard errors:\n")
    print.default(sd, print.gap = 2,  quote = FALSE) 
    cat("\n")
    if (nrow(object@gini) > 1) {
      bm <- names(which.min(apply(object@gini, 1, max)))
    } else {
      bm <- colnames(object@gini)[which.max(object@gini[1, ])]
    }
    cat(paste("The selected score is ", bm, ".\n", sep = ""))  
  }  
)

setMethod("plot", 
  signature(x = "gini", y = "missing"),
  function(x, y, overlay = TRUE, ...) {
    # -Wall: no binding or global variables in R check 
    Base <- Premium <- Score <- Loss <- NULL
    lrz <- lapply(x@lorenz, as.data.frame)
    pd <- lapply(1:length(lrz), function(t){
            lrz[[t]]$Base <- rep(names(lrz)[t], nrow(lrz[[t]]))
            melt(lrz[[t]], c("Base", ".P.")) 
            })
    pd <- do.call("rbind", pd)
    names(pd) <- c("Base", "Premium", "Score", "Loss")
    pd$Score <- factor(pd$Score, levels = eval(x@call$score))
    pp <- ggplot(pd, aes(Premium, Loss))
    pp <- pp + geom_line(aes(color = Score, linetype = Score)) + 
           geom_line(aes(Premium, Premium), color = gray(0.35))
    if (overlay)
      pp <- pp + facet_grid(Base~.) else 
      pp <- pp + facet_grid(Base~Score) 
    print(pp)
  }
)
