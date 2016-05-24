SSlbc <- 
selfStart( ~ 10^(b0 + b1*theta),
   initial = function(mCall, LHS, data)
   {
      xy <- sortedXyData(mCall[["theta"]], LHS, data)
      if(nrow(xy) < 3) {
         stop("Too few distinct theta values to fit a lbc model!")
      }
      z <- xy[["y"]]
      z[z == 0] <- 0.001
      xy[["z"]] <- log10(z)
      aux <- coef(lm(z ~ x, data = xy))
      pars <- coef(nls(y ~ 10^(b0 + b1*x),
         data = xy, start = list(b0 = aux[1], b1 = aux[2])))
      names(pars) <- mCall[c("b0", "b1")]
      return(pars)
   },
   parameters = c("b0", "b1") )