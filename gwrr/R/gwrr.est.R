gwrr.est <-
function(form, locs, data, kernel="exp", bw=TRUE, rd=TRUE, cv.tol){   # User-called function
   # Parse variables in formula to pass to function
   lhs <- as.character(form)[2]
   rhs <- as.character(form)[3]
   rhs.v <- strsplit(rhs, " + ", fixed=TRUE)   # Returns a list with 1 first element, unknown 2 elements
   n.l <- length(rhs.v[[1]])   # get number of x variables

   # Create y vector and design matrix
   db <- data
   y <- db[,lhs]
   N <- dim(db)[1]
   X <- rep(1,N)   # Assume an intercept
   for(i in 1:n.l) X <- cbind(X, db[,rhs.v[[1]][i]])

   # Calculate pairwise distances
   library(fields)   
   S <- rdist(locs)   # Assume Euclidean distance is appropriate for now

   rmspe <- NA   # RMSPE for CV, return NA if bandwidth and ridge are input
      
   # Check what should be estimated
   if (bw==TRUE){   # Estimate kernel bandwidth
      band.ub <- ceiling(max(S))
      band.lb <- min(S) + 0.01 * band.ub   # Add a small amount to min(S) to have non-zero value; ad hoc
      if(missing(cv.tol)){
         lm1 <- lm(form, data=db)
         lm.rmse <- gwr.rmse(y, lm1$fitted.values)
         cv.tol <- lm.rmse * 0.05    # Set CV tolerance as small % of RMSE from linear model; ad hoc
      }
      if (rd==TRUE){   # Estimate ridge parameter
         lam.lb <- 0
         lam.ub <- 2   # Set reasonable range for shrinkage parameter; ad hoc
         lam.tol <- 0.01 * (lam.ub-lam.lb)      
         g.bw <- gwrr.bw.rd.cv(band.lb, band.ub, cv.tol, lam.lb, lam.ub, lam.tol, X, y, S, kernel)
         bw <- g.bw$phi
         rd <- g.bw$lambda
         rmspe <- g.bw$RMSPE
      }
   }

   # Call estimation functions
   g.beta <- gwrr.beta(bw, rd, X, y, S, N, kernel)
   g.yhat <- gwr.yhat(g.beta, X)
   g.rmse <- gwr.rmse(y, g.yhat)
   g.rsquare <- gwr.rsquare(y, g.yhat)
   
   # Return estimates
   params <- list(bw, rd, rmspe, g.beta, g.yhat, g.rmse, g.rsquare)
   names(params) <- c("phi", "lambda", "RMSPE", "beta", "yhat", "RMSE", "rsquare")
   params
}

