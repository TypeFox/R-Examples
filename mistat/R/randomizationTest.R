randomizationTest <- function(list, R = 500, 
                              calc, fun=NA, seed=NA, 
                              printSummary=TRUE){
  
  list <- as.list(list)
  list <- lapply(list, na.omit)
  calc <- as.function(calc)
  
  if(!is.function(fun) && !is.na(fun))
    fun <- as.function(fun)
  if(length(list) < 1)
    stop("argument list should be of lenght >= 1")
  
  len <- sapply(X=list, FUN=length)
  dat <- stack(list)
  
  calcs <- function(x, len, calc){
    x[,2] <- as.factor(rep(1:length(len), len))
    by(data=x[,1], INDICES=x[,2], FUN=calc, simplify=T)
  }
  stats <- function(x, i, len, fun, calc){
    x <- calcs(x[i,], len, calc)
    if(is.function(fun))
      return(fun(x))
    else
      return(x)
  }
  
  if(!is.na(seed))
    set.seed(seed)
  
  b <- boot::boot(data=dat, statistic=stats, R=R, 
            len=len, fun=fun, calc=calc)
  
  res <- data.frame(D = c(b$t, b$t0), 
                    O=c(rep(FALSE, R), TRUE), 
                    I = -1)
  res[order(res$D), "I"] <- 1:nrow(res)
  
  if(printSummary){
    cat(paste("Original stat is at quantile", 
              res[res$O == TRUE, "I"], 
              "over", R+1, 
              " (", 
              round((res[res$O == TRUE, "I"]/(R+1))*100, 2), 
              "%) ", "\n"))
    cat(paste("Original stat is", 
              round(res[res$O == TRUE, "D"], 6), "\n"))
  }
  invisible(b)
}