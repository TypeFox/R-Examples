# functions for fanchart
# last change MS Friday, August 13, 2010 at 09:40
# a new argument is added 
#------------------------------------------------------------------------------------------
centiles.fan<-function (obj, 
                    xvar = NULL, 
                    cent = c(0.4, 2, 10, 25, 50, 75, 90, 98, 99.6),
                    ylab = "y", 
                    xlab = "x", 
                    main = NULL, 
               main.gsub = "@", 
                    xleg = min(xvar), #Note main and a substitution marker (an experiment!)
                    yleg = max(obj$y), 
                    xlim = range(xvar), 
                    ylim = range(obj$y), 
                  points = FALSE,
                  median = TRUE,  
                     pch =  15,
                     cex = 0.5,
                     col = gray(0.7),  
                  colors = c("cm","gray", "rainbow", "heat", "terrain", "topo"),  
                           ...)        
{
   if (!is.gamlss(obj))  stop(paste("This is not an gamlss object", "\n", ""))
   if(is.null(xvar)) stop(paste("The xvar argument is not specified", "\n", ""))
      colors <- match.arg(colors) 
       fname <- obj$family[1]
        qfun <- paste("q",fname,sep="")
       Title <- paste("Centile fanchart using", fname, sep = " ") #Note change to default title handling
        main <- if (is.null(main)) paste("Centile curves using", fname, sep=" ") 
                else gsub(main.gsub,Title,main)    
       oxvar <- xvar[order(xvar)]
       oyvar <- obj$y[order(xvar)]  
       if (is.matrix(obj$y)) # Monday, March 26, 2007 at 14:12
          {oyvar <- obj$y[,1][order(xvar)] 
           ylim  <- range(obj$y[,1])
            yleg <- max(obj$y[,1])
          }
          plot(oxvar, oyvar, type = "n", col = col, pch = pch, 
               xlab = xlab, ylab = ylab, xlim = xlim, ylim, ...) 
          title(main)     
        lpar <- length(obj$parameters)
          ii <- 0
         #per <- rep(0,length(cent))
    LL <- matrix(0, ncol=length(cent), nrow=length(xvar))
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
     LL[,ii] <- eval(newcall)
     }
          xx <- c(oxvar, rev(oxvar)) # the x variable in the poligon
          ll <- floor(dim(LL)[2]/2)
       color <- switch(colors, 
                  "cm"=cm.colors(ll) ,
                  "gray"=rev(gray(sqrt(seq(from = 0.1, to = .85, length = ll)))), 
                  "rainbow"=rainbow(ll, start=.5, end=.9), 
                  "heat"=rev(heat.colors(ll)), 
                  "terrain"=rev(terrain.colors(ll)), 
                  "topo"=rev(topo.colors(ll)))      # getting the type of color scheam
         ii <- 0
    for (i in 1:ll)
    {
         yy <- c(LL[,i], rev(LL[,dim(LL)[2]-ii]))
         polygon(xx,yy, col=color[i], border=color[i])
         ii <- ii+1
    }
 if (points==TRUE)
    {
          points(oxvar, oyvar,  col = col , cex = cex, pch = pch,  ...) 
    } 
 if (median==TRUE)
    {
        if(lpar==1) 
        {
        newcall <-call(qfun,.5,
                    mu=fitted(obj,"mu")[order(xvar)]) 
        }
     else if(lpar==2)
        {
        newcall <-call(qfun,.5,
                    mu=fitted(obj,"mu")[order(xvar)],
                    sigma=fitted(obj,"sigma")[order(xvar)]) 
        }
     else if(lpar==3)
        {
        newcall <-call(qfun,.5,
                    mu=fitted(obj,"mu")[order(xvar)],
                    sigma=fitted(obj,"sigma")[order(xvar)],
                    nu=fitted(obj,"nu")[order(xvar)])
        }
     else 
        {
        newcall <-call(qfun,.5,
                    mu=fitted(obj,"mu")[order(xvar)],
                    sigma=fitted(obj,"sigma")[order(xvar)],
                    nu=fitted(obj,"nu")[order(xvar)],
                    tau=fitted(obj,"tau")[order(xvar)]) 
        }     
        med <- eval(newcall)
        lines(oxvar, med,  col = "black" , pch =pch,  ...)  
    }     
}
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
