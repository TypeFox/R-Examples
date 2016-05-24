# this function is created for comparing centiles from 
# different models. MS Tuesday, February 17, 2004 at 10:28
centiles.com  <- function( obj, 
                       ...,
                       xvar = NULL,  
                       cent = c(.4,10,50,90,99.6), 
                     legend = TRUE, 
                       ylab = "y", 
                       xlab = "x", 
                       xleg = min(xvar), 
                       yleg = max(obj$y), 
                       xlim = range(xvar), 
                       ylim = NULL, 
                    no.data = FALSE,
                      color = TRUE,
                       main = NULL, 
                       plot = TRUE
                       )
{
if (length(list(...)))   # more than one fitted model  
  {  
      object <- list(obj, ...)
        nobj <- length(object)
    isgamlss <- unlist(lapply(object, is.gamlss))
      if (!any(isgamlss)) stop("some of the objects are not gamlss")
      if(is.null(xvar)) stop(paste("The xvar argument is not specified", "\n", ""))
     #   type <- unlist(lapply(object, function(x) x$type=="Continuous"))
     # if (!any(type)) stop("The centiles are working only with continuous distributions")# ms Sunday, April 2, 2006 
       fname <- lapply(object, function(x) x$family[1])      
        qfun <- lapply(fname, function(x) paste("q",x,sep=""))
      lenpar <- lapply(object, function(x) length(x$parameters) )
       oxvar <- xvar[order(xvar)]
       oyvar <- object[[1]]$y[order(xvar)]
      if (is.null(ylim)) ylim <- range( object[[1]]$y) 
       Title <- if (is.null(main)) paste("Centile curves") else main
      if (plot)
        {
        if (no.data==FALSE) type<-"p" else type<-"n"
        plot(oxvar, oyvar, type=type, pch = 15, cex = 0.5, col = gray(0.7),
             xlab= xlab, ylab=ylab, xlim=xlim, ylim=ylim)
        title(Title)#  
        }
       ltype <- 0 
     for (iii in 1:nobj) # over models
       { 
          cat("********  Model", iii,"******** \n" )
           lpar <- lenpar[[iii]]
         if (color==TRUE)  col <- 3  else col <- 1  
          ltype <- ltype+1
             ii <- 0
             per <- rep(0,length(cent))    
         for(var in cent) 
           { 
            if(lpar==1) 
              {
               newcall <-call(qfun[[iii]],var/100,
                    mu=fitted(object[[iii]],"mu")[order(xvar)]) 
              }
            else if(lpar==2)
              {
               newcall <-call(qfun[[iii]],var/100,
                   mu=fitted(object[[iii]],"mu")[order(xvar)],
                   sigma=fitted(object[[iii]],"sigma")[order(xvar)]) 
              }
            else if(lpar==3)
              {
               newcall <-call(qfun[[iii]],var/100,
                    mu=fitted(object[[iii]],"mu")[order(xvar)],
                    sigma=fitted(object[[iii]],"sigma")[order(xvar)],
                    nu=fitted(object[[iii]],"nu")[order(xvar)])
              }
            else 
              {
               newcall <-call(qfun[[iii]],var/100,
                    mu=fitted(object[[iii]],"mu")[order(xvar)],
                    sigma=fitted(object[[iii]],"sigma")[order(xvar)],
                    nu=fitted(object[[iii]],"nu")[order(xvar)],
                    tau=fitted(object[[iii]],"tau")[order(xvar)]) 
              }
              ii <- ii+1
               ll<- eval(newcall)
           if (plot)
            { 
            lines(oxvar,ll,col=col, lty=ltype)
            if (color==TRUE)  colleg <- c(3,4,5,6,7,8,9,10) else colleg <- c(1)  
            if (legend==TRUE) legend(list(x=xleg,y=yleg), legend = cent, 
              col=colleg, lty=1, ncol=1, bg="white")# 
            } 
           if (color==TRUE)  col <- col+1  #1
           per[ii]<-(1-sum(oyvar>ll)/length(oyvar))*100
           cat("% of cases below ", var,"centile is ", per[ii], "\n" )
         }
                       
      } 
  }
else
  {  
   if (!is.gamlss(obj))  stop(paste("This is not an gamlss object", "\n", ""))
   if(is.null(xvar)) stop(paste("The xvar argument is not specified", "\n", ""))
#   if(!obj$type=="Continuous") stop(paste("The centiles are working only with continuous distributions", "\n", ""))
       fname <- obj$family[1]
        qfun <- paste("q",fname,sep="")
       Title <- paste("Centile curves using",fname, sep=" ")
       oxvar <- xvar[order(xvar)]
       oyvar <- obj$y[order(xvar)] 
        if (plot)
        {
     if (no.data==FALSE) type <- "p" else type <- "n"    
    plot(oxvar, oyvar, type = type ,  pch =  15, cex = 0.5, col =  gray(0.7),
           xlab = xlab, ylab = ylab ,xlim = xlim, ylim, ...)#
                  title(Title)#  , cex.main = 1.1
        }          
          if (color==TRUE)  col <- 3 else col <- 1 
        lpar <- length(obj$parameters)
          ii <- 0
         per <- rep(0,length(cent))
      for(var in cent) 
        { 
          if(lpar==1) 
        {
        newcall <-call(qfun,var/100,
                    mu=fitted(obj,"mu")[order(xvar)]) 
        }
      else if(lpar==2)
        {
        newcall <-call(qfun,var/100,
                    mu=fitted(obj,"mu")[order(xvar)],
                    sigma=fitted(obj,"sigma")[order(xvar)]) 
        }
      else if(lpar==3)
        {
        newcall <-call(qfun,var/100,
                    mu=fitted(obj,"mu")[order(xvar)],
                    sigma=fitted(obj,"sigma")[order(xvar)],
                    nu=fitted(obj,"nu")[order(xvar)])
        }
      else 
        {
        newcall <-call(qfun,var/100,
                    mu=fitted(obj,"mu")[order(xvar)],
                    sigma=fitted(obj,"sigma")[order(xvar)],
                    nu=fitted(obj,"nu")[order(xvar)],
                    tau=fitted(obj,"tau")[order(xvar)]) 
        }
          ii <- ii+1
          ll<- eval(newcall)
         if (plot)
          { 
          lines(oxvar,ll,col=col, lty=1)
          if (color==TRUE)  colleg <- c(3,4,5,6,7,8,9,10) else colleg <- c(1)  
          if (legend==TRUE) legend(list(x=xleg,y=yleg), legend = cent, 
             col=colleg, lty=1, ncol=1, bg="white")# ,
          } 
           if (color==TRUE)  col <- col+1  #
          per[ii]<-(1-sum(oyvar>ll)/length(oyvar))*100
          cat("% of cases below ", var,"centile is ", per[ii], "\n" )
         }
  }
}
