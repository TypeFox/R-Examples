makedata <-
function(seed=1024, intercept=1, beta=1, n=4, g=20, areaID=NULL, ve=1, ve.contam=41, ve.epsilon=0, vu=1, vu.contam=41, vu.epsilon=0){
   # prepare the area-size-specific issues
   if (is.null(areaID)){
      N <- n*g 
      areaID <- rep(1:g, each=n)
      n <- rep(n, g)
   } else{
      n <- as.vector(table(areaID))
      g <- length(n)
      N <- length(areaID)
   }
   # generate the fixed effects
   #number of regressonrs, excl. intercept
   p <- length(beta)
   if (is.null(intercept)){
      x <- matrix(NA, N, p)
      beta.names <- paste(rep("x", p), 1:p, sep="")
      hasintercept <- 0
      for (i in 1:(p)){
	 x[, i] <- rnorm(N)
      }
   }else{
      beta <- c(intercept, beta)
      beta.names <- c("(Intercept)", paste(rep("x", p), 1:p, sep=""))
      x <- matrix(NA, N, (p+1))
      x[, 1] <- rep(1, N)
      hasintercept <- 1
      for (i in 2:(p+1)){
	 x[, i] <- rnorm(N)
      }
   }
   # get the final p
   p <- dim(x)[2]
   # give the columns of x names
   colnames(x) <- beta.names
   # generate y
   y <- as.vector(as.matrix(x) %*% beta)
   # add a random error for the model error from the Huber-Tukey mixture        
   if (ve.epsilon == 0){
      y <- y + rnorm(N, 0, sqrt(ve))
   }else{
      # contaminated observations
      outliers <- sample(1:N, floor(ve.epsilon*N))
      y[outliers] <- y[outliers] + rnorm(length(outliers), 0, sqrt(ve.contam))
      # un-contaminated observations
      non.outliers <- setdiff(1:N, outliers)
      y[non.outliers] <- y[non.outliers] + rnorm(length(non.outliers), 0, sqrt(ve))
   }
   # y as list
   y.list <- split(y, areaID)
   # add area-specific variation
   addraneff <- function(u, s){
      raneff <- rnorm(1, 0, s)
      u <- u + rep(raneff, length(u))
   }
   outlyingAreas <- sample(1:g, floor(g*vu.epsilon))
   if (length(outlyingAreas) == 0){
      r <- lapply(y.list, addraneff, s=sqrt(vu))
   }else{
      non.outlyingAreas <- setdiff(1:g, outlyingAreas)     
      r <- as.list(1:g)
      r[outlyingAreas] <- lapply(y.list[outlyingAreas], addraneff, s=sqrt(vu.contam))
      r[non.outlyingAreas] <- lapply(y.list[non.outlyingAreas], addraneff, s=sqrt(vu))
   } 
   # y as vector (not list) again
   y <- unsplit(r, areaID) 
   # prepare return value
   res <- list(y=y, X=x, areaID=areaID, nsize=n, g=g, p=p, n=N, intercept=hasintercept) 
   attr(res, "areaNames") <- paste(rep("A", g), 1:g, sep="")
   attr(res, "areadef") <- paste("area-specific ranef") 
   attr(res, "yname") <- "y"
   attr(res, "xnames")<- beta.names 
   attr(res, "call") <- match.call() 
   class(res) <- "saemodel"
   # additional attributes (do not belong to the "saemodel" class)
   skeleton <- list(intercept = intercept,
     beta = beta, 
     ve = ve,
     n = n,
     g = g,
     vu = vu,
     ve.contam = ve.contam,
     ve.epsilon = ve.epsilon,
     vu.contam = vu.contam,
     vu.epsilon = vu.epsilon) 
   attr(res, "contam") <- list(skeleton=skeleton, outlyingAreas=outlyingAreas)
   return(res)
}

