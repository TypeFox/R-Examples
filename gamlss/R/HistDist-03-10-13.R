#---------------------------------------------------------------------------------------
# 2-12-2011 MS: the gamlssML() is tried first before using gamlss(). This should be faster
# gamlss() is used if gamlssML() failed.
# MS+CA Saturday, April 22, 2006 
# last 11-07-07 :... is taken out from barplot 
# last 12-2-10 : ammend to cope with the right frequency for ZIP ans ZAP type distributions
histDist <- function(y, 
                  family = NO, 
                    freq = NULL,
                 density = FALSE,
                   nbins = 10,
                    xlim = NULL, 
                    ylim = NULL,
                    main = NULL,
                    xlab = NULL, 
                    ylab = NULL, 
                    data = NULL,
                      ... )
{
             FA <- as.gamlss.family(family)    
          fname <- FA$family[1]
           dfun <- paste("d",fname,sep="")
           lpar <- length(FA$parameters)
       typeDist <- FA$type
    #if (!is.null(data)) {attach(data); on.exit(detach(data))}
         subsY <- substitute(y)
             y <- if (!is.null(data)) get(deparse(substitute(y)), envir=as.environment(data)) else y
          freq <- if (!is.null(data)&&!is.null(freq)) get(deparse(substitute(freq)), envir=as.environment(data)) else freq  
  switch(typeDist, 
     "Continuous"=
           {
         extra <- (max(y,na.rm = TRUE)-min(y, na.rm = TRUE))/5
          xmin <- if (is.null(xlim))  min(y,na.rm = TRUE)-extra else xlim[1]
          xmax <- if (is.null(xlim))  max(y,na.rm = TRUE)+extra else xlim[2]
           # for the possitive distributions start from 0.001  
           if (!FA$y.valid(xmin)) xmin <- 0.001
           if (!FA$y.valid(xmax)) xmax <-0.999 # for the beta      y
                x1 <- seq(xmin,xmax,length=101)
           },
      "Discrete" =
           {
            if (fname %in% .gamlss.bi.list)# if binomial
              {
               if (NCOL(y) == 1) 
                 {
                #  if (is.factor(y))   
                   x1 <- c(0,1)
                   bd <- rep(1, 2)
                  # y1 <- y 
                 }
                else                 
                 {
                 if (any(abs(y - round(y)) > 0.001)) {warning("non-integer counts in a binomial GAMLSS!")}
                 bd <- y[, 1] + y[, 2]
                 y1 <- y[, 1]
                 if (any(y1 < 0 | y1 > bd)) stop("y values must be 0 <= y <= N")
               xmin <- if (is.null(xlim[1]))  
                           {if (all(bd==bd[1])) min(y1, na.rm = TRUE) else 0 }
                       else xlim[1]    
               xmax <- if (is.null(xlim)) 
                           {if (all(bd==bd[1])) max(bd ,na.rm = TRUE) else 1 }
                       else xlim[2]
                 x1 <- seq(xmin,xmax,by=1)
                 }
              }
            else
              { # if other discrete
          xmin <- if (is.null(xlim))  min(y,na.rm = TRUE) else xlim[1]
          xmax <- if (is.null(xlim))  max(y,na.rm = TRUE) else xlim[2]
            x1 <- seq(xmin,xmax,by=1)     
              }
           },
      "Mixed"=
           {
              stop("Mixed distributions are not implemented yet")
           })          
           # fit the model  
            if (is.null(freq)) 
            {  # no frequencies
                 # if (is.null(data)) 
                    {
                      mod <- try(gamlssML(y,family=fname), silent=TRUE)
                            if (any(class(mod)%in%"try-error")) 
                                { # if gamlssML fails use gamlss
                                  mod <-  try( gamlss(y~1, family=fname,  ...))
                                }
                   
                    }
                  #else
                  #  {
                   #  mod <- try(gamlssML(y,family=fname, data=data), silent=TRUE)
                  #          if (any(class(mod)%in%"try-error")) 
                  #              { # if gamlssML fails use gamlss
                  #                mod <-  try( gamlss(y~1, family=fname, data=data, ...))
                  #              }
                  #  }
             }
            else
             { # with frequencies
                 # .freq <-freq
                  # if (is.null(data)) 
                     {
                       mod <- try(gamlssML(y, weights=freq, family=fname, ...), silent=TRUE)
                            if (any(class(mod)%in%"try-error")) 
                                { # if gamlssML fails use gamlss
                                  mod <-  try( gamlss(y~1, weights=freq, family=fname,  ...))
                                }
                     }
                   #else 
                  #   {
                   #    mod <- try(gamlssML(y, weights=freq, family=fname, data=data, ...), silent=TRUE)
                    #        if (any(class(mod)%in%"try-error")) 
                    #            { # if gamlssML fails use gamlss
                    #              mod <-  try( gamlss(y~1, weights=freq, family=fname, data=data,  ...))
                     #           }
                    # }
             } 
             
             
             
                         mod$call$family <- eval(as.expression(fname)) 
                   if (mod$method=="BFGS"||mod$method=="nlminb") 
                    {
                     mod$call$formula <-  subsY
                    }   
                 else  
                   {
                     mod$call$formula[[2]] <-  subsY
                   }   
       if (!is.null(data)) mod$call$data <- substitute(data)
           # get the pdf
     switch(lpar, 
          "1" =  {
             newcall <- if ((fname %in% .gamlss.bi.list)) 
                               call(dfun,x1, mu = fitted(mod)[1], bd=bd[1] )
                        else   call(dfun,x1, mu = fitted(mod)[1])
                 },
          "2" =  {
             newcall <- if ((fname %in% .gamlss.bi.list))
                         {  
  call(dfun,x1, mu = fitted(mod)[1],sigma = fitted(mod,"sigma")[1], bd=bd[1])

                         }  
                       else                            
                         {
                            call(dfun,x1, mu = fitted(mod)[1], sigma = fitted(mod,"sigma")[1])  
                         }
                 },
          "3" = {      
            newcall <- if ((fname %in% .gamlss.bi.list))
                         {  
                         call(dfun,x1, mu = fitted(mod)[1],sigma = fitted(mod,"sigma")[1],
                                        nu = fitted(mod,"nu")[1], bd=bd[1])

                         }  
                 else    {call(dfun,x1, mu = fitted(mod)[1], 
                                  sigma = fitted(mod,"sigma")[1], 
                                     nu = fitted(mod,"nu")[1])
                         }
                 },
          "4" =
                {
               newcall <- if ((fname %in% .gamlss.bi.list))
                         {  
                         call(dfun,x1, mu = fitted(mod)[1],
                                    sigma = fitted(mod,"sigma")[1],
                                       nu = fitted(mod,"nu")[1],
                                      tau = fitted(mod,"tau")[1],
                                       bd = bd[1])

                         }  
                 else    {call(dfun, x1, mu = fitted(mod)[1], 
                                      sigma = fitted(mod,"sigma")[1],
                                         nu = fitted(mod,"nu")[1],
                                        tau = fitted(mod,"tau")[1])
                         } 
                })
switch(typeDist,  
    "Continuous"=  
              { 
                y1 <- eval(newcall)
                xlim <- c(xmin,xmax)
           #===================================================================
           #MIKIS functions in version of histDist (incl. the discrete case) 
           #truehist(y, xlim = mm, main=main, ...)   #MIKIS OLD function
           #lines(x1,y1, col="blue3") #MIKIS OLD function
           #===================================================================
           #Poppy's Function -  Date: Friday 21, 2006 at 13:46
           #main1 <- paste("Histogram of y-variable and ", sep = "") #main <- paste("Histogram of y and a fitted (GAMLSS family), ", family," distribution to y", sep = "")
           #main2 <- paste("the fitted ", FA$family[2]," distribution to y-variable", sep = "") #main <- paste("Histogram of y and a fitted (GAMLSS family), ", family," distribution to y", sep = "")
           # main <- c(main1, main2) 
              main <-  if (is.null(main))  
                      paste("The ", deparse(subsY), " and the fitted ", FA$family[1]," distribution", sep = "")
                      else main
              if (!is.null(freq)) y <- rep(y,freq) # in case frequencies are used
              xlab <- if (is.null(xlab)) deparse(subsY) else xlab
              ylab <- if (is.null(ylab)) paste("f()") else ylab
            truehist(y, xlim = xlim, ylim = ylim, xlab=xlab, ylab=ylab, nbins=nbins, main = main, col="gray", lty=3, 
                     border="blue", fg=rainbow(12)[9], col.main = "blue4", col.lab = "blue4", col.axis = "blue")
                     #Poppy's Function
            lines(x1, y1, col = "red", col.axis = "blue", col.main = "blue4", col.lab = "blue4", 
                       lwd=2,fg = gray(0.7))
            if (density==TRUE) 
                 { dens<-density(y)
                  lines(dens$x, dens$y, col="blue") 
                 }
                   
           #===================================================================
              },
      "Discrete"=
              {  
               if (fname %in% .gamlss.bi.list) # if binomial 
                  {                  
                      if (all(bd==bd[1]))     # binary + equal bd
                        { 
                        y1 <- eval(newcall) # find the fitted values 
                        tabley <- if (NCOL(y) == 1) table(y) else  table(y[,1])
                        dft <- data.frame(tabley)# get the table  
                        r<-barplot(tabley/sum(tabley), fg = "blue", col="gray", axis.lty=1, 
                                     border="blue", col.axis = "blue4", col.main = "blue4", 
                                     col.lab = "blue4", ylim=ylim, xlim=xlim, ylab=ylab, 
                                     xlab=xlab, main=main)
                        yy1 <- y1[ x1%in%dft[[1]]] # this is to make sure that we have as many points as in the bar plot
                        lines(r, yy1, type='h', col = "red",  lwd=2, col.axis = "blue", 
                                    col.main = "blue4", col.lab = "blue4", fg = gray(0.7))  
                        points(r, yy1,  col="red", pch=21, lwd=2,col.axis = "blue") # extra points  
                        }
                      else  
                        {
                        xlim <- c(xmin, xmax)
                        main <- paste("proportions of ", deparse(subsY), sepp = "")
                        
                        main <-  if (is.null(main))  paste("proportions of ", deparse(subsY), sepp = "")
                                 else main
                        xlab <- if (is.null(xlab)) deparse(subsY) else xlab
                        ylab <- if (is.null(ylab)) paste("f()") else ylab
                        truehist(y1/bd, xlim = xlim, xlab="proportions", ylab=ylab, ylim=ylim, 
                                nbins=nbins, col="gray", 
                                lty=3, border="blue", col.axis = "blue",col.main = "blue4", 
                                col.lab = "blue4", main=main, fg=rainbow(12)[9])#P
                        #plot(y1/bd~I(y1/bd), type="h", xlim=mm, main=main)
                        #points(y1/bd, y1/bd,  cex= bd/mean(bd))
                        lines( fitted(mod), rep(1,length(fitted(mod))), type="h", col="red" )
                        points(fitted(mod), rep(0, length(fitted(mod))), col = "red", pch=25)
                        points(fitted(mod),rep(1,length(fitted(mod))),   col="red", pch=24)
                        }

                }  
             else
                {
                 y1 <- eval(newcall) # find the fitted values 
                 fy <- if (is.null(freq)) factor(y,levels=x1) else factor(rep(y,freq),levels=x1)
            notresy <- if (is.null(freq)) factor(y) else factor(rep(y,freq))# RR DS Friday, February 12, 2010 at 17:41
                dft <-data.frame(tabley<- xtabs(~fy))# get the table    
               #===================================================================
               #MIKIS functions in version of histDist (incl. the discrete case)
               #r <- barplot(tabley/sum(tabley), col="lightblue",  ...)## rainbow(20))#rainbow(20
               #yy1<- y1[ x1%in%dft[[1]]] # this is to make sure that we have as many points as in the bar plot
               #lines(r,  yy1, type='h', col="red", lwd=3) # fitted lines 
               #points(r, yy1,  col="red", pch=21, lwd=3) # extra points
               #===================================================================
               #Poppy's Function -  Date: Thursday 20, 2006 at 13:46
            # main1 <- paste("Histogram of y-variable and ", sep = "") #main <- paste("Histogram of y and a fitted (GAMLSS family), ", family," distribution to y", sep = "")
            # main2 <- paste("the fitted ", FA$family[2]," distribution to y-variable", sep = "") #main <- paste("Histogram of y and a fitted (GAMLSS family), ", family," distribution to y", sep = "")
            #  main <- c(main1, main2)
               main <-  if (is.null(main))  
                      paste("Barplot of the ", deparse(subsY), " and the fitted ", FA$family[2]," distribution", sep = "")
                      else main 
                 r <-  barplot(tabley/ sum(xtabs(~notresy)),  # Friday, February 12, 2010 at 17:41
                                    fg = "blue", col="gray", axis.lty=1, 
                                    border="blue", col.axis = "blue4",col.main = "blue4", 
                                    col.lab = "blue4", main=main, ylim=ylim, ylab=ylab, xlab=xlab ) 
                      #if (is.null(freq)) 
                      #    { barplot(tabley/sum(tabley), 
                      #              fg = "blue", col="gray", axis.lty=1, 
                      #              border="blue", col.axis = "blue4",col.main = "blue4", 
                      #              col.lab = "blue4", main=main, ...)
                      #    }
                      #  else
                      #    {
                      #      barplot(freq/sum(freq),  names.arg = y, 
                      #              fg = "blue", col="gray", axis.lty=1, 
                      #              border="blue", col.axis = "blue4",col.main = "blue4", 
                      #              col.lab = "blue4", main=main, ...)
                      #    }        
                            yy1<- y1[ x1%in%dft[[1]]] # this is to make sure that we have as many points as in the bar plot
                            lines(r, yy1, type='h', col = "red",  lwd=2, col.axis = "blue", 
                                    col.main = "blue4", col.lab = "blue4", fg = gray(0.7))  
                            points(r, yy1,  col="red", pch=21, lwd=2,col.axis = "blue") # extra points  
                 }       
              },
        "Mixed"=
              {
              stop("Mixed distributions are not implemented yet")
              }) 
         #on.exit(if (!is.null(freq)) rm(.freq,pos=1))          
            mod
}
                 
