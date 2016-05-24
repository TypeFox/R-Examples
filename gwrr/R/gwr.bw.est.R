gwr.bw.est <-
function(form, locs, data, kernel="exp", cv.tol){
   # Parse variables in formula to pass to function
   lhs <- as.character(form)[2]
   rhs <- as.character(form)[3]
   rhs.v <- strsplit(rhs, " + ", fixed=TRUE)   # Returns a list with 1 first element, unknown 2 elements
   n.l <- length(rhs.v[[1]])   # get number of x variables

   # Create y vector and design matrix
   db <- data
   y <- db[,lhs]
   N <- dim(db)[1]
   X <- rep(1,N)   # Assume intercept for now
   for(i in 1:n.l) X <- cbind(X, db[,rhs.v[[1]][i]])
   
   # Calculate pairwise distances
   library(fields)   
   S <- rdist(locs)   # Assume Euclidean distance is appropriate for now
   
   # Set boundaries and tolerances for CV
   band.ub <- ceiling(max(S))
   band.lb <- min(S) + 0.01 * band.ub   # Add a small amount to min(S) to have non-zero value; ad hoc
   if(missing(cv.tol)){
      lm1 <- lm(form, data=db)
      lm.rmse <- gwr.rmse(y, lm1$fitted.values)
      cv.tol <- lm.rmse * 0.05    # Set CV tolerance as small % of RMSE from linear model; ad hoc
   }
   a <- band.lb
   b <- band.ub
   c <- (a+b)/2
   diff <- b - a
   N <- dim(X)[1]
   
   while (diff > cv.tol){
      a.c <- (a+c)/2
      c.b <- (c+b)/2
      RMSE.a.c <- gwr.cv.err(a.c, X, y, S, N, kernel)
      RMSE.c.b <- gwr.cv.err(c.b, X, y, S, N, kernel)

      if (RMSE.a.c < RMSE.c.b){
         b <- c.b
         RMSE.b <- RMSE.c.b
         print(paste("Bandwidth: ", format(b,digits=4), " RMSPE :", format(RMSE.b,digits=4)))
      }

      if (RMSE.a.c > RMSE.c.b){
         a <- a.c
         RMSE.a <- RMSE.a.c
         print(paste("Bandwidth: ", format(a,digits=4), " RMSPE :", format(RMSE.a,digits=4)))
      }

      c <- (a+b)/2
      diff <- abs(b - a)
   }

   RMSE.lb <- gwr.cv.err(band.lb, X, y, S, N, kernel)
   RMSE.ub <- gwr.cv.err(band.ub, X, y, S, N, kernel)
   RMSE.c <- gwr.cv.err(c, X, y, S, N, kernel)

   # Check bounds
   if (RMSE.lb < RMSE.c){
      c <- band.lb
      RMSE.c <- RMSE.lb
   }
   if (RMSE.ub < RMSE.c){
      c <- band.ub
      RMSE.c <- RMSE.ub
   }
   print(paste("Bandwidth: ", format(c,digits=4), " RMSPE :", format(RMSE.c,digits=4)))
   params <- list(c, RMSE.c, RMSE.c * N)
   names(params) <- c("phi", "RMSPE", "cv.score")
   params
}

