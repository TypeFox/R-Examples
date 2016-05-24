################################################################################
#                                                                              #
# utility and value function package                                           #
# ==================================                                           #
#                                                                              #
# version 1.3                                        Peter Reichert 05.10.2014 #
#                                                                              #
################################################################################


# ==============================================================================
# plotting functions (see also the object-specific plotting functions)
# ==============================================================================


utility.endnode.plot1d <- 
                   function(node,
                            col       = utility.calc.colors(),
                            gridlines = c(0.2,0.4,0.6,0.8),
                            main      = "",
                            cex.main  = 1,
                            ...)
{
   length <- 101
   x <- seq(node$range[1],node$range[2],length=length)
   u <- evaluate(node,attrib=x)
   title <- main; if ( nchar(title) == 0 ) title <- node$name
   funtype <- "utility"; if ( !node$utility ) funtype <- "value"
   plot(numeric(0),numeric(0),type="l",
        xlim=node$range,ylim=c(0,1),
        xlab=node$attrib,ylab=funtype,main=title,
        xaxs="i",yaxs="i",xaxt="n",yaxt="n",cex.main=cex.main,...)

   # colored bar along y axis:
   
   if ( length(col)>1 & !node$utility )
   {
      num.grid = 100
     
      # y-axix:
      endpoints <- seq(0,1,length.out=num.grid+1)+1/(2*num.grid)
      midpoints <- 0.5*(endpoints[-1]+endpoints[-length(endpoints)])
      cols <- utility.get.colors(midpoints,col)
      for ( i in 1:(num.grid-1) )
      {
         lines((node$range[1]-0.01*(node$range[2]-node$range[1]))*c(1,1),
               endpoints[c(i,i+1)],
               col=cols[i],lwd=3,lend=2,xpd=TRUE)
      }
     
      # x-axis:
      midpoints <- 0.5*(u[-1]+u[-length(u)])
      cols <- utility.get.colors(u,col)
      for ( i in 1:length(midpoints) )
      {
         lines(c(x[i],x[i+1]),
               -0.01*c(1,1),
               col=cols[i],lwd=3,lend=2,xpd=TRUE)
      }
   }
   
   # axes (should overly colored bar):
   
   axis(side=1)
   axis(side=2)
   
   # plot gridlines:
   
   if ( !node$utility )
   {
      if ( ! is.na(gridlines[1]) )
      {
         for ( level in gridlines )
         {
            abline(h=level,lty="dashed")
            for ( i in 1:(length-1) )
            {
               if ( !is.na(u[i]) & !is.na(u[i+1]) )
               {
                  if ( (u[i] <= level & u[i+1] > level) |
                       (u[i] > level & u[i+1] <= level) )
                  {
                     x.level <- x[i] + (level-u[i])/(u[i+1]-u[i])*(x[i+1]-x[i])
                     lines(c(x.level,x.level),c(0,level),lty="dashed")
                  }
               }
            }
         }
      }
   }
   
   # plot value/utility function:
   
   color <- "black"
   if ( length(col) == 1 ) color <- col
   lines(x,u,lwd=2,col=color)
}


utility.endnode.plot2d <- function(node,
                                   col       = utility.calc.colors(),
                                   gridlines = c(0.2,0.4,0.6,0.8),
                                   main      = "",
                                   cex.main  = 1,
                                   ...)
{
   num.grid <- 100
   x <- node$ranges[[1]][1] + 
        ((1:num.grid)-0.5)/num.grid*(node$ranges[[1]][2]-node$ranges[[1]][1])
   y <- node$ranges[[2]][1] + 
        ((1:num.grid)-0.5)/num.grid*(node$ranges[[2]][2]-node$ranges[[2]][1])
   
   array.x <- sort(rep(x,num.grid))
   array.y <- rep(y,num.grid)
   array.xy <- cbind(array.x,array.y)
   colnames(array.xy) <- node$attrib
   
   u <- evaluate(node,as.data.frame(array.xy))
   u <- t(matrix(u,ncol=num.grid,byrow=FALSE))
   
   title <- main; if ( nchar(title) == 0 ) title <- node$name
   image(x=x,y=y,z=u,xlim=node$ranges[[1]],ylim=node$ranges[[2]],zlim=c(0,1),
         col=col,xlab=node$attrib[1],ylab=node$attrib[2],main=title,
         cex.main=cex.main)
}


utility.conversion.plot <- function(node,
                                    col       = "black",
                                    gridlines = NA,
                                    cex.main  = 1,
                                    ...)
{
   length <- 101
   x <- ((1:length)-1)/(length-1)
   u <- evaluate.cond(node,x)
   plot(numeric(0),numeric(0),type="l",
        xlim=c(0,1),ylim=c(0,1),
        xlab=paste("value(",node$nodes[[1]]$name,")",sep=""),ylab="utility",
        main=node$name,xaxs="i",yaxs="i",cex.main=cex.main)
   color <- "black"; if ( length(col) == 1 ) color <- col
   lines(x,u,lwd=2,col=color)
   lines(c(0,1),c(0,1))
   if ( length(node$x) > 0 & length(node$u) > 0 )
   {
      if ( length(node$x) == length(node$u) )
      {
         points(node$x,node$u,cex=1.5,xpd=TRUE)
      }
   }
}


utility.aggregation.plot <- function(node           = node,
                                     col            = col,
                                     gridlines      = gridlines,
                                     cex.main       = 1,
                                     cex.attrib     = 1,
                                     cex.nodes      = 1,
                                     ...)
{
  nodes.names <- rep(NA,length(node$nodes))
  for ( i in 1:length(node$nodes) ) nodes.names[i] <- node$nodes[[i]]$name
  if ( length(node$nodes) == 2 )
  {
    num.grid <- 100
    x <- ((1:num.grid)-0.5)/num.grid
    y <- ((1:num.grid)-0.5)/num.grid
    
    array.x <- sort(rep(x,num.grid))
    array.y <- rep(y,num.grid)
    array.xy <- cbind(array.x,array.y)
    
    v <- apply(array.xy,1,node$name.fun,node$par)
    v <- t(matrix(v,ncol=num.grid,byrow=FALSE))
    
    if ( node$utility )
    {
      contour(x=x,y=y,z=v,levels=gridlines,xlim=c(0,1),ylim=c(0,1),zlim=c(0,1),
              axes=FALSE,add=FALSE,lty="solid",lwd=2,
              xlab=node$nodes[[1]]$name,ylab=node$nodes[[2]]$name,
              main=node$name,...)
    }
    else
    {
      # area coloring:
      
      image(x=x,y=y,z=v,xlim=c(0,1),ylim=c(0,1),zlim=c(0,1),
            col=col,xaxt="n",yaxt="n",
            xlab=node$nodes[[1]]$name,ylab=node$nodes[[2]]$name,
            main=node$name,...)

      # colored bar along axes:
      
      endpoints <- seq(0,1,length.out=num.grid+1)+1/(2*num.grid)
      midpoints <- 0.5*(endpoints[-1]+endpoints[-length(endpoints)])
      cols <- utility.get.colors(midpoints,col)
      for ( i in 1:(num.grid-1) )
      {
        lines(-0.01*c(1,1),endpoints[c(i,i+1)],col=cols[i],lwd=3,lend=2,xpd=TRUE)
        lines(endpoints[c(i,i+1)],-0.01*c(1,1),col=cols[i],lwd=3,lend=2,xpd=TRUE)
      }
      
      # axes (should overly colored bar):
      
      axis(1)
      axis(2)
      lines(c(1,1,0),c(0,1,1))
      
      # contour lines:
      
      contour(x=x,y=y,z=v,levels=gridlines,xlim=c(0,1),ylim=c(0,1),zlim=c(0,1),
              axes=FALSE,add=TRUE,lty="solid",lwd=2,...)
    }
  }
  else
  {
    if ( node$name.fun == "utility.aggregate.add" |
         node$name.fun == "utility.aggregate.geo" |  
         node$name.fun == "utility.aggregate.cobbdouglas" |  
         node$name.fun == "utility.aggregate.harmo")
    {
      type <- "Additive"
      if ( node$name.fun == "utility.aggregate.geo" |  
             node$name.fun == "utility.aggregate.cobbdouglas" ) type = "Geometric"
      if ( node$name.fun == "utility.aggregate.harmo" ) type = "Harmonic"  
      w <- node$par/sum(node$par)
      w.max <- max(w)
      if ( length(w) != length(nodes.names) )
      {
        warning("Node \"",node$name,"\": ",
                "length of sub-nodes and weights not equal: ",
                length(nodes.names)," ",length(w),sep="")
      }
      else
      {
        barplot(w,names.arg=nodes.names,ylim=c(0,1.2*w.max),
                ylab="weight",main=node$name,cex.main=cex.main,cex.names=cex.nodes)
        text(0.5*1.3*length(w),1.1*w.max,paste(type,"aggregation with weights:"))
      }
    }
    else
    {
      if ( node$name.fun == "utility.aggregate.mult" )
      {
        w <- node$par
        w.max <- max(w)
        if ( length(w) != length(nodes.names) )
        {
          warning("Node \"",node$name,"\": ",
                  "length of sub-nodes and weights not equal: ",
                  length(nodes.names)," ",length(w),sep="")
        }
        else
        {
          barplot(w,names.arg=nodes.names,ylim=c(0,1.2*w.max),
                  ylab="weight",main=node$name,cex.main=cex.main,cex.names=cex.nodes)
          text(0.5*1.3*length(w),1.1*w.max,
               "Multiplicative aggregation with weights:")
        }
      }
      else
      {
        if ( node$name.fun == "utility.aggregate.min" |
             node$name.fun == "utility.aggregate.max" )
        {
          type <- "Minimum (worst-case)"
          if ( node$name.fun == "utility.aggregate.max" ) type <- "Maximum"
          plot(numeric(0),numeric(0),xlim=c(0,1),ylim=c(0,1),
               xaxt="n",yaxt="n",main=node$name,xlab="",ylab="",
               cex.main=cex.main)
          text(0.5,0.9,paste(type,"aggregation of nodes:"))
          for ( i in 1:length(nodes.names) )
          {
            text(0.5,0.7*i/length(nodes.names),nodes.names[i])
          }
        }
        else
        {
          plot(numeric(0),numeric(0),xlim=c(0,1),ylim=c(0,1),
               xaxt="n",yaxt="n",main=node$name,xlab="",ylab="",
               cex.main=cex.main)
          text(0.5,0.9,paste("aggregation with function \"",
                             node$name.fun,"\" of nodes:",sep=""))
          for ( i in 1:length(nodes.names) )
          {
            text(0.5,0.7*i/length(nodes.names),nodes.names[i])
          }
        }
      }
    }
  }
}


utility.plotcolbox <- function(x,y,col,val=NA,plot.val=FALSE)
{
  # check for availability of data:
  
  if ( length(val) == 0 ) return()
  if ( is.na(val[1]) & length(col)>1 ) return()
  
  # plot colored box (without border):
  
  color <- col
  if ( length(col) > 1 ) color <- utility.get.colors(val[1],col)
  polygon(x      = c(x[1],x[2],x[2],x[1],x[1]),
          y      = c(y[1],y[1],y[2],y[2],y[1]),
          col    = color,
          border = NA)
  
  # optionally plot value line:
  
  if ( plot.val & !is.na(val[1]) )
  {
    lines((x[1]+val[1]*(x[2]-x[1]))*c(1,1),y,lwd=1.0)
  }
}


utility.plotquantbox <- function(x,y,col,val,num.stripes=500)
{
  min.halfwidth <- 0.02
  
  # check for availability of data:
  
  if ( length(val) == 0 ) return()
  if ( sum(is.na(val)) == length(val) ) return()
  
  # get quantiles:
  
  quant <- quantile(val[!is.na(val)],probs=c(0.05,0.5,0.95))
  if ( quant[3]-quant[1] < 2*min.halfwidth )
  {
    quant[1] <- max(0,quant[1]-min.halfwidth)
    quant[3] <- min(1,quant[3]+min.halfwidth)
  }
  
  # plot colored quantile box:
  for ( j in floor(num.stripes*quant[1]):ceiling(num.stripes*quant[3]) )
  {
    lines((x[1]+j/num.stripes*(x[2]-x[1]))*c(1,1),y,
          col=utility.get.colors(j/num.stripes,col))
  }
  
  # plot median line:
  
  lines((x[1]+quant[2]*(x[2]-x[1]))*c(1,1),y,lwd=1.5)
  
  # return:
  
  return()
}


utility.plothierarchy <- 
   function(node,
            u           = NA,
            uref        = NA,
            col         = utility.calc.colors(),
            main        = "",
            cex.main    = 1,
            cex.nodes   = 1,
            cex.attrib  = 1,
            with.attrib = TRUE,
            levels      = NA,
            plot.val    = TRUE,
            ...)
{
   # call multiple times if u and possibly uref are lists:
     
   if ( is.list(u) & !is.data.frame(u) )
   {
      if ( is.list(uref) & !is.data.frame(uref) )
      {
         if ( length(u) == length(uref) )
         {
            for ( i in 1:length(u) )
            {
               utility.plothierarchy(node        = node,
                                     u           = u[[i]],
                                     uref        = uref[[i]],
                                     col         = col,
                                     main        = main,
                                     cex.main    = cex.main,
                                     cex.nodes   = cex.nodes,
                                     cex.attrib  = cex.attrib,
                                     with.attrib = with.attrib,
                                     levels      = levels,
                                     plot.val    = plot.val,
                                     ...)
            }
         }
         else
         {
            warning("if u and uref are lists, their lengths must be equal")
         }
      }
      else
      {
        utility.plothierarchy(node        = node,
                              u           = u[[i]],
                              uref        = uref,
                              col         = col,
                              main        = main,
                              cex.main    = cex.main,
                              cex.nodes   = cex.nodes,
                              cex.attrib  = cex.attrib,
                              with.attrib = with.attrib,
                              levels      = levels,
                              plot.val    = plot.val,
                              ...)
      }
      return()
   }
     
   # global parameters:

   delta.x        <- 0.1
   delta.y        <- 0.1
   dh.rel.utility <- 0.1

   # get hierarchy structure and define positions of boxes:
         
   str <- utility.structure(node)
   if ( ! is.data.frame(str) )
   {
      warning("unable to identify structure of objectives hierarchy")
      return()
   }
   if ( !is.na(levels) ) str <- utility.prune(str,levels)
   w <- 1/max(str$level)
   if ( with.attrib ) w <- 1/(max(str$level)+1)
   h <- 1/str$endnodes[1]
   str$x <- (str$level-0.5)*w
   str$y <- 1-(str$offset+0.5*str$endnodes)*h
   x.attrib <- max(str$level)*w + delta.y*w 

   # convert u and uref to data frames:

   u.local <- u
   if ( is.vector(u.local) ) u.local <- t(u.local)         
   u.local <- as.data.frame(u.local)
   uref.local <- uref
   if ( is.vector(uref.local) ) uref.local <- t(uref.local)         
   uref.local <- as.data.frame(uref.local)
   
   # plot indvidual plots per row if the same number of titles is provided;
   # plot quantile summary if not the same number of titles is provided and 
   # if the number of rows is > 1
   
   quant.summary <- length(main) != nrow(u.local) & nrow(u.local) > 1
   
   # find out if u and uref are available (otherwise plot required/not required shading)

   u.available <- FALSE
   if ( nrow(u.local)>1 | ncol(u.local)>1 | !is.na(u.local[1,1]) )
   {
      u.available <- TRUE
   }
   uref.available <- FALSE 
   ind.uref.local <- rep(1,nrow(u.local))
   if ( nrow(uref.local)>1 | ncol(uref.local)>1 | !is.na(uref.local[1,1]) )
   {
      uref.available <- TRUE
      if ( !quant.summary ) # number of rows must be unity or equal to nrow(u)
      {
         if ( nrow(uref.local) == nrow(u.local) )
         {
            ind.uref.local <- 1:nrow(u.local)
         }
         else
         {
            if ( nrow(uref.local) != 1 ) uref.available <- FALSE
         }
      }
   }
   
   # loop over rows of utilities/values:

   num.plots <- nrow(u.local)
   if ( !u.available | quant.summary ) num.plots <- 1
   for ( k in 1:num.plots )
   {
      # set-up plot frame:
         
      par.def <- par(no.readonly=TRUE)
      par(mar=c(0,0,0,0))
      plot(numeric(0),numeric(0),xlim=c(0,1),ylim=c(0,1),
           xaxt="n",yaxt="n",xlab="",ylab="",cex.main=cex.main)
           
      # write title
      
      title <- main[1]
      if ( length(main) == nrow(u.local) ) title <- main[k]
      text(0,1-0.5*h,title,adj=c(0,0.5),cex=cex.main,...)
      
      # draw color code legend:
      
      if ( u.available )
      {
         x.l <- delta.x*w
         x.r <- (1-delta.x)*w
         y   <- 0.8*h
         num.col <- 100
         v <- (1:num.col - 0.5)/num.col
         colors <- utility.get.colors(v,col)
         for ( i in 1:num.col ) 
         {
            lines(x.l+(x.r-x.l)/num.col*c(i-1,i),c(y,y),col=colors[i],lwd=3,lend=2)
         }
         text(x.l,y,"0",pos=1,cex=cex.nodes)
         text(x.r,y,"1",pos=1,cex=cex.nodes)
      }
      
      # loop over all boxes in the hierarchy:
      
      for ( i in 1:nrow(str) )
      {
         # calculate box edge coordinates:
            
         x  <- str$x[i] + (0.5-delta.x)*w*c(-1,1)
         y  <- str$y[i] + (0.5-delta.y)*h*c(-1,1)
         y1 <- c(0.5*(y[1]+y[2]),y[2])   # upper part, uref
         y2 <- c(y[1],0.5*(y[1]+y[2]))   # lower part, u
         
         # plot background color or quantile boxes:
            
         if ( !u.available ) # plot required/not required nodes in differnt grey
         {
            if ( str$required[i] ) color <- grey(0.7)
            else                   color <- grey(0.9)
            utility.plotcolbox(x,y,color)
         }
         else
         {
            if ( !quant.summary ) # plot hierarchy for each row of u
            {
               # plot background color and vertical line:
              
               val <- u.local[k,rownames(str)[i]]
               color <- col
               if ( str$utility[i] ) color <- "white"
               if ( !uref.available )
               {
                 utility.plotcolbox(x,y,color,val,plot.val)
               }
               else
               {
                 valref <- uref.local[ind.uref.local[k],rownames(str)[i]]
                 utility.plotcolbox(x,y1,color,valref,plot.val)
                 utility.plotcolbox(x,y2,color,val,plot.val)                 
               }
            }
            else # plot quantile summary of v or expected u
            {
               if ( !str$utility[i] ) # plot quantile summary
               {
                  val <- u.local[,rownames(str)[i]]
                  if ( !uref.available )
                  {
                    utility.plotquantbox(x,y,col,val)
                  }
                  else
                  {
                    valref <- uref.local[,rownames(str)[i]]
                    utility.plotquantbox(x,y1,col,valref)
                    utility.plotquantbox(x,y2,col,val)
                  }                 
               }
               else   # plot expected utility
               {
                  u.exp <- NA
                  column <- match(rownames(str)[i],colnames(u.local))
                  if ( !is.na(column) )
                  {
                    u.exp <- mean(u.local[,column],na.rm=TRUE)
                  }
                  if ( !uref.available )
                  {
                    utility.plotcolbox(x,y,"white",u.exp)
                  }
                  else
                  {
                    uref.exp <- NA
                    column <- match(rownames(str)[i],colnames(uref.local))
                    if ( !is.na(column) )
                    {
                      uref.exp <- mean(uref.local[,column],na.rm=TRUE)
                    }
                    col1 <- "lightgreen"
                    col2 <- "tomato"
                    if ( u.exp > uref.exp )
                    {
                      col1 <- "tomato"
                      col2 <- "lightgreen"
                    }
                    utility.plotcolbox(x,y1,col1,uref.exp)
                    utility.plotcolbox(x,y2,col2,u.exp)                    
                  }
               }
            }
         }
                        
         # plot bounding box:

         lines(x   = c(x[1],x[2],x[2],x[1],x[1]),
               y   = c(y[1],y[1],y[2],y[2],y[1]),
               col = as.character(str$color[i]))
         if ( str$utility[i] )
         {
            dh <- dh.rel.utility*(y[2]-y[1])
            lines(x,(y[1]+dh)*c(1,1))
            lines(x,(y[2]-dh)*c(1,1))
         }
                  
         # write text into box:
            
         text(str$x[i],str$y[i],rownames(str)[i],cex=cex.nodes,...)

         # plot connecting lines:
                           
         upper <- str$upper[i]
         if ( ! is.na(upper) )
         {
            x.line.l <- str[upper,"x"] + (0.5-delta.x)*w
            x.line.r <- str$x[i] - (0.5-delta.x)*w
            x.line.v <- str[upper,"x"] + 0.5*w
            y.line.l <- str[upper,"y"]
            y.line.r <- str$y[i]
            lines(x = c(x.line.l,x.line.v,x.line.v,x.line.r), 
                  y = c(y.line.l,y.line.l,y.line.r,y.line.r))
         }
            
         # write attribute names:
                 
         if ( with.attrib )
         {
            if ( str$endnode[i] )
            {
               attributes <- strsplit(str$attributes[i],split=";")[[1]]
               n <- length(attributes)
               for ( j in 1:n )
               {
                  y.attrib <- str$y[i] +  (0.5 - (j-0.5)/n)*(1-delta.y)*h
                  text(x.attrib,y.attrib,attributes[j],pos=4,cex=cex.attrib,...)
                  lines(c(x[2],x.attrib),c(y.attrib,y.attrib),lty="dotted")
               }
            }
         }
      } # end for i
      par(par.def)
   } # end for k
}


utility.plottable <- 
   function(node,
            u,
            uref       = NA,
            nodes      = NA,
            col        = utility.calc.colors(),
            main       = "",
            cex.main   = 1,
            cex.nodes  = 1,
            f.reaches  = 0.2,
            f.nodes    = 0.2,
            levels     = NA,
            plot.val   = FALSE,
            print.val  = TRUE,
            ...)
{
   # global parameters:

   delta.x        <- 0.2
   delta.y        <- 0.2
   delta.main     <- 0.05
   dh.rel.utility <- 0.1
   
   # initializations:
   
   if ( !is.list(u) )
   {
      warning("unable to interpret u")
      return()
   }
   if ( length(nodes)==1 & is.na(nodes[1]) ) nodes <- character(0)
   str <- utility.structure(node)
   if ( !is.na(levels) )
   {
     if ( is.data.frame(str) )
     {
       str1 <- utility.prune(str,levels) 
       ind <- order(str1$level)
       nodes <- unique(c(nodes,rownames(str1)[ind][str1$level[ind]<=levels]))
     }
   }
   uref.available <- FALSE
   ind.uref <- NA
   uref.local <- uref
   if ( is.data.frame(u) | is.matrix(u) )
   {
     if ( length(nodes)==0 ) nodes <- colnames(u)
     reaches <- rownames(u)
     if ( is.data.frame(uref) | is.matrix(uref) )
     {
       if ( nrow(u) == nrow(uref) )
       {
         uref.available <- TRUE
         ind.uref <- 1:nrow(uref)
       }
       else
       {
         if ( nrow(uref) == 1 )
         {
           uref.available <- TRUE
           ind.uref <- rep(1,nrow(u))
         }
       }
     }
   }
   else
   {
     if( length(nodes)==0 ) nodes <- colnames(u[[1]])
     reaches <- names(u)
     if ( is.list(uref) | is.matrix(uref) )
     {
       if ( !is.data.frame(uref) & !is.matrix(uref) )
       {
         if ( length(uref) == length(u) )
         {
           ind.uref <- 1:length(u)
           uref.available <- TRUE
         }
         else
         {
           if ( length(uref) == 1 )
           {
             ind.uref <- rep(1,length(u))
             uref.available <- TRUE
           }
         }
       }
       else
       {
         uref.local <- list()
         uref.local[[1]] <- uref
         ind.uref <- rep(1,length(u))
         uref.available <- TRUE
       }
     }
   }
   
   # set-up plotting parameters and plot frame:

   dx <- (1-f.reaches)/length(nodes)
   dy <- (1-f.nodes)/length(reaches)
   x <- f.reaches+(1:length(nodes)-0.5)*dx
   y <- 1-f.nodes-(1:length(reaches)-0.5)*dy
   if ( nchar(main[1]) > 0 )
   {
      y  <- (1-delta.main)*y
      dy <- (1-delta.main)*dy
   }
   par.def <- par(no.readonly=TRUE)
   par(mar=c(0,0,0,0))
   plot(numeric(0),numeric(0),xlim=c(0,1),ylim=c(0,1),
        xaxt="n",yaxt="n",xlab="",ylab="")
   
   # write and color values:
   
   for ( i in 1:length(reaches) )
   {
      for ( j in 1:length(nodes) )
      {
         xbox <- x[j]+0.5*(1-delta.x)*dx*c(-1,1)
         ybox <- y[i]+0.5*(1-delta.y)*dy*c(-1,1)
         if ( is.data.frame(u) | is.matrix(u) )
         {
           if ( !is.na(match(reaches[i],rownames(u))) &
                !is.na(match(nodes[j]  ,colnames(u))) )
           {
             yb <- ybox; if ( uref.available ) yb[2] <- 0.5*(ybox[1]+ybox[2])
             yt <- y[i]; if ( uref.available ) yt <- y[i] - 0.25*(ybox[2]-ybox[1])
             val <- u[reaches[i],nodes[j]]
             color <- col
             if ( !is.na(match(nodes[j],rownames(str))) )
             {
               if ( str[nodes[j],"utility"] ) color <- "white"
             }
             utility.plotcolbox(xbox,yb,color,val=val,plot.val=plot.val)
             if ( !is.na(val) & print.val )
             {
               val.str <- paste(round(val,2))
               if ( nchar(val.str) > 1 & substring(val.str,1,1) == "0" )
               {
                 val.str <- substring(val.str,2)
                 if ( nchar(val.str) == 2 ) val.str <- paste(val.str,"0",sep="")
               }
               text(x=x[j],y=yt,val.str,cex=cex.nodes) 
             }
           }
           if ( uref.available )
           {
             if ( !is.na(match(nodes[j],colnames(uref))) )
             {
               yb <- ybox; if ( uref.available ) yb[1] <- 0.5*(ybox[1]+ybox[2])
               yt <- y[i]; if ( uref.available ) yt <- y[i] + 0.25*(ybox[2]-ybox[1])
               val <- uref[ind.uref[i],nodes[j]]
               color <- col
               if ( !is.na(match(nodes[j],rownames(str))) )
               {
                 if ( str[nodes[j],"utility"] ) color <- "white"
               }
               utility.plotcolbox(xbox,yb,color,val=val,plot.val=plot.val)
               if ( !is.na(val) & print.val )
               {
                 val.str <- paste(round(val,2))
                 if ( nchar(val.str) > 1 & substring(val.str,1,1) == "0" )
                 {
                   val.str <- substring(val.str,2)
                   if ( nchar(val.str) == 2 ) val.str <- paste(val.str,"0",sep="")
                 }
                 text(x=x[j],y=yt,val.str,cex=cex.nodes) 
               }
             }
           }
         }
         else
         {
           yb <- ybox; if ( uref.available ) yb[2] <- 0.5*(ybox[1]+ybox[2])
           if ( !is.na(match(reaches[i],names(u))) &
                !is.na(match(nodes[j],colnames(u[[reaches[i]]]))) )
           {
             val <- u[[reaches[i]]][,nodes[j]]
             utility.plotquantbox(xbox,yb,col,val,num.stripes=500)
           }
           if ( uref.available )
           {
             yb <- ybox; yb[1] <- 0.5*(ybox[1]+ybox[2])
             val <- uref.local[[ind.uref[i]]][,nodes[j]]
             if ( length(val) > 1 )
             {
               utility.plotquantbox(xbox,yb,col,val,num.stripes=500)
             }
           }
         }
         
         # plot bounding box:
         
         lines(x   = c(xbox[1],xbox[2],xbox[2],xbox[1],xbox[1]),
               y   = c(ybox[1],ybox[1],ybox[2],ybox[2],ybox[1]),
               col = as.character(str$color[j]))
         if ( !is.na(match(nodes[j],rownames(str))) )
         {
           if ( str[nodes[j],"utility"] )
           {
             dh <- dh.rel.utility*(ybox[2]-ybox[1])
             lines(xbox,(ybox[1]+dh)*c(1,1))
             lines(xbox,(ybox[2]-dh)*c(1,1))
           }
         }
      }
   }

   # write title and names of nodes and reaches:
   
   if ( nchar(main[1]) > 0 ) text(x=0.5,y=1-0.5*delta.main,label=main[1],cex=cex.main)
   for ( i in 1:length(reaches) ) 
   {
     text(x=0,y=y[i],label=reaches[i],adj=c(0,0.5),cex=cex.nodes)
   }
   
   par(srt=90)
   for ( j in 1:length(nodes) ) 
   {
     text(x=x[j],y=1-f.nodes,label=nodes[j],adj=c(0,0.5),cex=cex.nodes)
   }
   par(srt=0)
   
   # reset plotting parameters:
   
   par(par.def)
}


utility.plot <- function(node,
                         u           = NA,
                         uref        = NA,
                         type        = c("hierarchy","table","node","nodes"),
                         nodes       = NA,
                         col         = utility.calc.colors(),
                         gridlines   = c(0.2,0.4,0.6,0.8),
                         main        = "",
                         cex.main    = 1,
                         cex.nodes   = 1,
                         cex.attrib  = 1,
                         f.reaches   = 0.2,
                         f.nodes     = 0.2,
                         with.attrib = TRUE,
                         levels      = NA,
                         plot.val    = TRUE,
                         print.val   = TRUE,
                         ...)
{
   if ( type[1] == "nodes" | type[1] == "node" )
   {
      # plot current node:
      
      if ( is.na(nodes[1]) | ! is.na(match(node$name,nodes)) )
      {
         if ( substring(class(node),1,18) == "utility.conversion" )
         {
            utility.conversion.plot(node       = node,
                                    col        = col,
                                    gridlines  = gridlines,
                                    cex.main   = cex.main,
                                    cex.nodes  = cex.nodes,
                                    cex.attrib = cex.attrib,
                                    ...)
         }
         else
         {
            if ( substring(class(node),1,19) == "utility.aggregation" )
            {
               utility.aggregation.plot(node       = node,
                                        col        = col,
                                        gridlines  = gridlines,
                                        cex.main   = cex.main,
                                        cex.nodes  = cex.nodes,
                                        cex.attrib = cex.attrib,
                                        ...)
            }
            else
            {
               if ( node$type == "endnode" )
               {
                  if ( class(node) == "utility.endnode.cond" )
                  {
                    plot(node$nodes[[i]],
                          par       = NA,
                          col       = col,
                          gridlines = gridlines,
                          cex.main  = cex.main,
                          nodes     = nodes,
                          ...)
                  }
                  else
                  {
                     plot(node$nodes[[i]],
                          par       = NA,
                          col       = col,
                          gridlines = gridlines,
                          cex.main  = cex.main,
                          ...)
                  }
               }
               else
               {
                  # unknown node type; not plotted
               }
            }
         }
      }
      
      # plot other nodes:
      
      if ( type == "nodes" )
      {
         if ( length(node$nodes) > 0 )
         {
         for ( i in 1:length(node$nodes) )
         {
            # initiate plot of subnodes:

            if ( node$nodes[[i]]$type == "endnode" )
            {
              if ( class(node$nodes[[i]]) == "utility.endnode.cond" )
              {
                plot(node$nodes[[i]],
                     par       = NA,
                     col       = col,
                     gridlines = gridlines,
                     cex.main  = cex.main,
                     nodes     = nodes,
                     ...)
              }
              else
              {
                 if ( is.na(nodes[1]) | ! is.na(match(node$nodes[[i]]$name,nodes)) )
                 {
                    plot(node$nodes[[i]],
                         par       = NA,
                         col       = col,
                         gridlines = gridlines,
                         cex.main  = cex.main,
                         ...)
                  }
               }
            }
            else
            {
               plot(node$nodes[[i]],
                    u          = u,
                    par        = NA,
                    type       = type,
                    nodes      = nodes,
                    col        = col,
                    gridlines  = gridlines,
                    cex.main   = cex.main,
                    cex.nodes  = cex.nodes,
                    cex.attrib = cex.attrib,
                    ...)
            }
         }
         }
      }
   }
   else
   {
      if ( type[1] == "hierarchy" )
      {
         if ( is.na(nodes[1]) | ! is.na(match(node$name,nodes)) )
         {
            utility.plothierarchy(node        = node,
                                  u           = u,
                                  uref        = uref,
                                  col         = col,
                                  main        = main,
                                  cex.main    = cex.main,
                                  cex.nodes   = cex.nodes,
                                  cex.attrib  = cex.attrib,
                                  with.attrib = with.attrib,
                                  levels      = levels,
                                  plot.val    = plot.val,
                                  ...)
         }
         if ( ! is.na(nodes[1]) )
         {
            if ( node$type != "endnode" )
            {
               for ( i in 1:length(node$nodes) )
               {
                  utility.plot(node$nodes[[i]],
                               u           = u,
                               uref        = uref,
                               type        = type,
                               nodes       = nodes,
                               col         = col,
                               gridlines   = gridlines,
                               main        = main,
                               cex.main    = cex.main,
                               cex.nodes   = cex.nodes,
                               cex.attrib  = cex.attrib,
                               with.attrib = with.attrib,
                               ...)
               } 
            }
         }
      }
      else
      {
         if ( type[1] == "table" )
         {
           utility.plottable(node       = node,
                             u          = u,
                             uref       = uref,
                             nodes      = nodes,
                             col        = col,
                             main       = main,
                             cex.main   = cex.main,
                             cex.nodes  = cex.nodes,
                             f.reaches  = f.reaches,
                             f.nodes    = f.nodes,
                             levels     = levels,
                             plot.val   = plot.val,
                             print.val  = print.val,
                             ...)
         }
         else
         {
            cat("unknown plot type:",type[1],"\n")
         }
      }
   }
}



# ==============================================================================





