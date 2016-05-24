# functions for centile plots for time series data
# last change MS Saturday, jan 24 2012
# needs zoo check why?
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
centiles.ts<-function (obj, 
                    xvar = NULL, 
                    cent = c(0.5, 2.5,  50, 95.5, 99.5 ),
                  legend = TRUE, 
                    ylab = "y", 
                    xlab = "x", 
                    main = NULL, 
               main.gsub = "@", 
                    xleg = min(xvar), #Note main and a substitution marker (an experiment!)
                    yleg = max(obj$y), 
                    xlim = range(xvar), 
                    ylim = range(obj$y), 
                    save = FALSE, 
                    plot = TRUE,
                    type = "l",
                  points = TRUE,  # this is new Saturday, August 14, 2010 MS
                     pch = "+", 
                     col = "blue", 
            col.centiles = 1:length(cent)+2, 
            lty.centiles = 1, 
            lwd.centiles = 1,  #Handling for line appearance 
                           ...)        
{
#	require(zoo)
  # if (!is.gamlss(obj))  stop(paste("This is not an gamlss object", "\n", ""))
   if(is.null(xvar)) stop(paste("The xvar argument is not specified", "\n", ""))
  # if(!obj$type=="Continuous") stop(paste("The centiles are working only with continuous distributions", "\n", ""))
  # if(any(family$family%in%.gamlss.bi.list))  we need to deal with binary response
  # the above is not needed MS Sunday, April 2, 2006 at 16:23
       fname <- obj$family[1]
        qfun <- paste("q",fname,sep="")
       Title <- paste("Centile curves using", fname, sep = " ") #Note change to default title handling
        main <- if (is.null(main)) paste("Centile curves using", fname, sep=" ") 
                else gsub(main.gsub,Title,main)    
       oxvar <- xvar[order(xvar)]
       oyvar <- obj$y[order(xvar)]  
       if (is.matrix(obj$y)) # Monday, March 26, 2007 at 14:12
          {oyvar <-  obj$y[,1][order(xvar)] 
           ylim  <-  range(obj$y[,1])
           yleg = max(obj$y[,1])
          }
       if (plot) 
          {
          lty.centiles <- rep(lty.centiles,length(cent))
          lwd.centiles <- rep(lwd.centiles,length(cent))
          col.centiles <- rep(col.centiles,length(cent))
          if (points==TRUE)
          {
          	y <- if (!is(obj$y, "zoo")) zoo(obj$y, order.by =xvar) else obj$y
         	plot(y, col="lightslategrey") 
         # plot(oxvar, oyvar, type = type, col = col, pch = pch, 
         #      xlab = xlab, ylab = ylab, xlim = xlim, ylim,  ...) 
          }
          else
          {
          plot(oxvar, oyvar, type = "n", col = col, pch = pch, 
               xlab = xlab, ylab = ylab, xlim = xlim, ylim, ...) 
          }
          title(main)
         }      
         col <- 3 # set this to 1 if you do not want colour 
        lpar <- length(obj$parameters)
          ii <- 0
         per <- rep(0,length(cent))
      # per <- cent/100 
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
     ll <- eval(newcall)
     if (plot)
        { 
         lines(oxvar, ll, col = col.centiles[ii], lty = lty.centiles[ii],lwd=lwd.centiles[ii], ...)
       # Modified by  Steve to include line types, colours etc
       # legend is moved below
       # if (legend==TRUE) legend(list(x=xleg,y=yleg), legend = cent, 
       #    col=c(3,4,5,6,7,8,9,10), lty=1, ncol=1, bg="white")# ,merge=TRUE)#, trace=TRUE)
        } 
      per[ii] <- (1-sum(oyvar>ll)/length(oyvar))*100
      if (!save) cat("% of cases below ", var,"centile is ", per[ii], "\n" )
     }
if (plot) 
     {         #Legend moved outside plotting loop
     if (legend == TRUE) 
         legend(list(x = xleg, y = yleg), legend = cent, 
                col = col.centiles, lty = lty.centiles, lwd=lwd.centiles, 
                ncol = 1, ...)
     }     
if (save) { return(cbind(cent,per))}
}
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
# more than one centiles plot 
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
# centiles.split will not work with less than 4 parameter distribution ??????

# centiles.split<- function(obj, 
#                         xvar = NULL, 
#                  xcut.points = NULL, 
#                      n.inter = 4 , 
#                         cent = c(.4,2,10,25,50,75,90,98,99.6), 
#                       legend = FALSE, 
#                         main = NULL,
#                    main.gsub = "@",   
#                         ylab = "y",
#                         xlab = "x",
#                         ylim = NULL, # ms Friday, February 6, 2004 at 19:40
#                      overlap = 0, 
#                         save = TRUE,
#                         plot = TRUE,
#                         ...
#                         )
# {
# # function 1 -------------------------
#   check.overlap <- function(interval)
#   {
#     if (!is.matrix(interval) ) {stop(paste("The interval specified is not a matrix."))}
#     if (dim(interval)[2] !=2) 
#      {
# stop(paste("The interval specified is not a valid matrix.\nThe number of columns should be equal to 2."))
#      }
#     crows = dim(interval)[1]
#     for (i in 1:(crows-1))
#     {
#         #if (interval[i,2] != interval[i+1,1]) {interval[i+1,1]=interval[i,2]}
#         if (!(abs(interval[i,2]-interval[i+1,1])<0.0001)) {interval[i+1,1]=interval[i,2]}
#     }
#     return(interval)
#   }
# #------------------------------------------
# # function 2 
#   get.intervals <- function (xvar, xcut.points ) 
# {
#     if (!is.vector(xcut.points))  {stop(paste("The interval is not a vector."))}
#     if ( any((xcut.points < min(xvar)) | any(xcut.points > max(xvar))))
#         {stop(paste("The specified `xcut.points' are not within the range of the x: (", min(xvar),
#         " , ", max(xvar), ")"))}
#     int <- c(min(xvar), xcut.points, max(xvar))
#      ii <- 1:(length(int)-1)
#       r <- 2:length(int)
#      x1 <- int[ii]
#      xr <- int[r]
#      if (any(x1>xr)) {stop(paste("The interval is are not in a increasing order."))}
#     cbind(x1,xr)
# }
# #------------------------------------------
# # the main function starts here
# #------------------------------------------
#     if (!is.gamlss(obj))  stop(paste("This is not an gamlss object", "\n", ""))
#     if(is.null(xvar)) 
#       { xvar <- seq(1,length(obj$y),1)
#         warning(paste("The xvar has been replace by an index 1:n", "\n", ""))  
#       }    # do we need this for centiles()            
#     if(is.null(xcut.points)) 
#       { # getting the intervals automatic
#       g.in <- co.intervals(xvar, number=n.inter, overlap=overlap )
#       if (overlap==0) g.in <- check.overlap(g.in) 
#       }                 
#     else
#       { # if xcut.points is set
#       g.in <- get.intervals(xvar, xcut.points ) 
#       }
#       howmany <- dim(g.in)[1]
#        fname <- obj$family[1]
#         qfun <- paste("q",fname,sep="")
#        Title <-vector(length = howmany)
#        for (i in 1:howmany)  Title[i] <- paste("(",paste(g.in[i,1],g.in[i,2], sep=", "), ")" )  
#        if (is.null(main)) main <- Title
#        else for (i in 1:howmany) main[i] <-gsub(main.gsub,Title[i],main[i])         
#              #  if (is.null(main)) 
#              #      {#rep(paste("Centile curves using", fname, sep=" "),howmany )
#              #       main<-vector(length = howmany)
#              #       for (i in 1:howmany)  main[i] <- paste("(",paste(g.in[i,1],g.in[i,2], sep=", "), ")" )
#              #      }               
#              #  else main #gsub(main.gsub,Title,main)    
#         if (length(main)!=howmany) rep(main, lenght=howmany) 
#        ncent <- length(cent)
#            X <- matrix(0, nrow = ncent, ncol = howmany, 
#                        dimnames=list(as.character(seq(1,ncent)),
#                        as.character(seq(1,howmany))))
#      dimnames(X)[[1]] <- as.character(cent)
#      dimnames(X)[[2]] <-  paste(substr(as.character(g.in[1:howmany,1]),1,7),"to",
#                                 substr(as.character(g.in[1:howmany,2]),1,7))
#         var1 <- ceiling(howmany/2)
#          
#        #   op <- par(mfrow=c(var1,2),  col.axis="blue4", col.main="blue4", 
#        #              col.lab="blue4",  col="darkgreen", bg="beige", ... )#
#         if (howmany==2)
#            layout(matrix(c(0,1,1,1,1,1,1,0,0,2,2,2,2,2,2,0), 8, 2))
#         else if (howmany==3)
#            layout(matrix(c(1,1,1,1,0,0,0,0,1,1,1,1,3,3,3,3,2,2,2,2,3,3,3,3,2,2,2,2,0,0,0,0), 8, 4))
#         else if (howmany==4)
#         layout(matrix(c(1,1,1,1,3,3,3,3,1,1,1,1,3,3,3,3,2,2,2,2,4,4,4,4,2,2,2,2,4,4,4,4), 8, 4))
#         else {
#            op <- par(mfrow=c(var1,2))#
#              }   
#        oxvar <- xvar[order(xvar)]
#        oyvar <- obj$y[order(xvar)] 
#    for (i in 1:howmany)
#        {
#           if(i==howmany) 
#             {   ##### Include points at the end of the last interval HB
#               yvar1 <- subset(oyvar, oxvar>=g.in[i,1]&oxvar<=g.in[i,2])
#               xvar1 <- subset(oxvar, oxvar>=g.in[i,1]&oxvar<=g.in[i,2])
#                 mu1 <- subset(obj$mu.fv[order(xvar)],oxvar>=g.in[i,1]&oxvar<=g.in[i,2])
#              sigma1 <- subset(obj$sigma.fv[order(xvar)],oxvar>=g.in[i,1]&oxvar<=g.in[i,2])
#                 nu1 <- subset(obj$nu.fv[order(xvar)], oxvar>=g.in[i,1]&oxvar<=g.in[i,2])
#                tau1 <- subset(obj$tau.fv[order(xvar)],oxvar>=g.in[i,1]&oxvar<=g.in[i,2])
#             }
#           else
#             {
#               yvar1 <- subset(oyvar, oxvar>=g.in[i,1]&oxvar<g.in[i,2])
#               xvar1 <- subset(oxvar, oxvar>=g.in[i,1]&oxvar<g.in[i,2]) 
#                 mu1 <- subset(obj$mu.fv[order(xvar)],oxvar>=g.in[i,1]&oxvar<g.in[i,2])
#              sigma1 <- subset(obj$sigma.fv[order(xvar)],oxvar>=g.in[i,1]&oxvar<g.in[i,2])
#                 nu1 <- subset(obj$nu.fv[order(xvar)], oxvar>=g.in[i,1]&oxvar<g.in[i,2])
#                tau1 <- subset(obj$tau.fv[order(xvar)],oxvar>=g.in[i,1]&oxvar<g.in[i,2])
#           }   
#               xlim1 <- c(g.in[i,1], g.in[i,2])
#               ylim1 <- if (is.null(ylim)) range(yvar1)
#                        else ylim      
#               xleg  <- min(xvar1) 
#               yleg  <- max(yvar1) 
#                obj1 <- obj
#              obj1$y <- yvar1
#          obj1$mu.fv <- mu1
#       obj1$sigma.fv <- sigma1
#          obj1$nu.fv <- nu1
#         obj1$tau.fv <- tau1
#               main1 <- main[i]
#           sampleval <- centiles(obj1, xvar1, cent, legend=legend, main=main1, main.gsub = main.gsub, 
#                            ylab, xlab ,xleg ,yleg, xlim=xlim1, ylim=ylim1, 
#                            save=TRUE, plot=plot, ...)
#               X[,i] <- sampleval[,2] 
#          }
#       if (howmany>4) par(op)
#       layout(mat=1)             
# if (save) { return(X)}
# }      
