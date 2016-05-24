qqPlotDemo <-
function(n = 25, distribution = "normal", mu = 0,sigma = 1, df = 10, 
            lambda = 10, numdf = 10, dendf = 16, shape1 = 40, shape2 = 5)
{      #get random sample from a distribution and plot its
        #histogram and normal quantile-quantile plot
      
        distr <- pmatch(distribution, c("normal", "t", "exponential", "chi.square", "F", "beta"), nomatch = NA)   
        if (is.na(distr))  stop("Distribution must be one of \"normal\", \"t\", \"exponential\", \"chi.square\", \"F\", or \"beta\" ")
                     
        if (sigma <=0|| df <=0 || lambda <= 0 || numdf <= 0 || dendf <= 0 || shape1 <= 0 || shape2 <= 0) stop("Parameter must be positive.")
        
        x <- switch(distr,  
        normal = rnorm(n, mu, sigma),
        t = rt(n, df),
        exponential = rexp(n, rate = lambda),
        chi.square = rchisq(n, df),
        F = rf(n, numdf, dendf),
        beta = rbeta(n, shape1, shape2))
        
        distr.expand <- char.expand(distribution, c("normal", "t", "exponential", "chi.square", "F", "beta"), nomatch = warning("No match"))

     par.old <- par(mfrow = c(2, 1), mar=c(2.1, 4.1, 2, 2), cex.main = .8, cex.axis = .8, cex.lab = .8)

      hist(x, main = "", xlab = "")
      title(paste("Sample size ", n , "; ", distr.expand, "distribution", sep=" "))
      par(mar=c(4.1,4.1,2,3))
      qqnorm(x)
      qqline(x)
      
     on.exit(par(par.old))
     invisible(x)

}
