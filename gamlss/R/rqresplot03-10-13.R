# no provision for "mixed" type of distribution is done here
# last change Tuesday, March 28, 2006 MS
 rqres.plot <- function (obj = NULL, 
                     howmany = 6, 
                        plot = c("all", "average"), 
                        type = c("wp", "QQ"), 
                        ...)
{
#-------------------------------
# local function
#-------------------------------	
  newres <- function(obj = NULL)
    {  if (!is.gamlss(obj))  stop(paste("This is not an gamlss object", "\n", ""))
      # if (obj$type!="Discrete" ) stop(paste("This is not discrete distribution ", "\n", ""))
       if (howmany>10&plot=="all")  stop(paste("You can only have 10 or less plots" , "\n", ""))
       w <- obj$weights
       if (all(trunc(w)==w)) # if frequenies as weights
           { 
            y  <- rep(obj$y, w)
            mu <- rep(fitted(obj, "mu"),w)
            if(any(obj$family%in%.gamlss.bi.list)){ bd <- rep(obj$bd,w)} # MS Wednesday, July 23, 2003 at 12:03   
            if ("sigma"%in%obj$parameters)  sigma <- rep(fitted(obj,"sigma"),w)
            if ("nu"%in%obj$parameters)        nu <- rep(fitted(obj,"nu"),w)
            if ("tau"%in%obj$parameters)      tau <- rep(fitted(obj,"tau"),w)  
           }
       else   # note that weights=1 and weights not frequencies are treated equal here and this could create problems in the future
           {
            y  <- obj$y
            mu <- fitted(obj)
            if(any(obj$family%in%.gamlss.bi.list)){ bd <- obj$bd} # MS Wednesday, July 23, 2003 at 12:03   
            if ("sigma"%in%obj$parameters)  sigma <- fitted(obj,"sigma")
            if ("nu"%in%obj$parameters)        nu <- fitted(obj,"nu")
            if ("tau"%in%obj$parameters)      tau <- fitted(obj,"tau")
           }
       res <- eval(obj$rqres)
     }
#-------------------------------
# local function finish here
#------------------------------- 
   plot <- match.arg(plot)
   type <- match.arg(type)
   var1 <- ceiling(howmany/2)
   lobj <- obj$noObs# length(fitted(obj))
     rs <- matrix(0,ncol=howmany, nrow =lobj ) 
 switch(plot,
        "all"=
        {
        op <- par(mfrow=c(var1,2), col.axis="blue4", col.main="blue4", col.lab="blue4",col="darkgreen", bg="beige")
         on.exit(par(op))
           for (i in 1:howmany)
             {
                   res <-  newres(obj)
               if (type=="QQ"){
                     rs[,i] <-    qqnorm(res, main = "Normal Q-Q Plot",
                                         xlab = "Theoretical Quantiles",
                                         ylab = "Sample Quantiles", 
                                         plot.it = TRUE, 
                                         #  frame.plot = TRUE, # MS taken out 
                                         col = "darkgreen", 
                                         #points(par(col="darkgreen"))# MS taken out
                     )$y
                     lines(res, res, col="red" , lwd=.4, cex=.4 )     
                           } 
              if (type=="wp"){
                rs[,i] <- res
                wp(resid=res, ...)
              }     
            }
        invisible(rs)     
         },
        "average"=
        {
          op <- par(col.axis="blue4", col.main="blue4", col.lab="blue4",col="darkgreen",bg="beige")
          on.exit(par(op))
        for (i in 1:howmany)
             {
                   res <-  newres(obj)
                rs[,i] <-    qqnorm(res, main = "Normal Q-Q Plot",
                          xlab = "Theoretical Quantiles",
                          ylab = "Sample Quantiles", 
                       plot.it = FALSE, 
                           #  frame.plot = TRUE, # MS taken out 
                           col = "darkgreen", 
                           #points(par(col="darkgreen"))# MS taken out
                           )$y
               }
           rmean <- apply(rs, MARGIN=1, "mean")
          if (type=="QQ"){
        	qqnorm(rmean, main = "Normal Q-Q Plot",
                xlab = "Theoretical Quantiles",
                ylab = "Sample Quantiles", 
                plot.it = TRUE, 
                frame.plot = TRUE, 
                col="darkgreen")
                lines(rmean, rmean, col="red" , lwd=.4, cex=.4 )
            }
          if (type=="wp"){
            wp(resid=rmean, ...)
          }
                invisible(rmean)     
        }) 
}
