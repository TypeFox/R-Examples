plot.PCAmix <- function(x,axes = c(1, 2), choice = "ind",label=TRUE,
                        coloring.ind=NULL,col.ind=NULL, coloring.var=NULL,
                        lim.cos2.plot=0,lim.contrib.plot=0, posleg="topleft",
                        xlim=NULL,ylim=NULL, cex=1,leg=TRUE,main=NULL,cex.leg=1, ...)
{
  cl<-match.call()
  if (!inherits(x, "PCAmix")) 
    stop("use only with \"PCAmix\" objects")
  
  res.pca <-x
  p1 <- res.pca$rec$p1
  p <- res.pca$rec$p
  p2<-res.pca$rec$p2
  m<-nrow(res.pca$levels$coord)
  quanti.coord <- res.pca$quanti$coord
  n<-nrow(res.pca$ind$coord)
  
  eig.axes<-res.pca$eig[axes,1]
  
  if (max(axes) > res.pca$ndim) 
    stop(paste("axes must be between 1 and ", res.pca$ndim, sep = ""))
  
  if (!(choice %in% c("ind", "sqload", "levels", "cor"))) 
    stop("\"choice\" must be either \"ind\",\"sqload\",\"cor\" or \"levels\"")  
  
  if (lim.cos2.plot != 0 & lim.contrib.plot!=0)
    stop("use either \"lim.cos2.plot\" OR \"lim.contrib.plot\"")
  
  if (!is.null(coloring.ind))
  {
    if (choice!="ind")
      warning("use \"coloring.ind\" only if choice=\"ind\"")
  }
  
  if (!is.null(coloring.ind))
  {
    if(!is.factor(coloring.ind) | length(coloring.ind)!=n)
      warning("\"coloring.ind\" must be either NULL or a qualitative variable of length equal to the number of individuals")
  }
  
  if (!is.null(coloring.var))
  {
    if (choice=="ind")
      warning("\"coloring.var\" is not used if choice=\"ind\"")
    if (coloring.var=="type")
    {
      if (choice=="cor" | choice=="levels")
        warning("\"coloring.var\" is not used if choice=\"cor\" or choice=\"levels\"")
    }
  }
  
  if (!is.null(coloring.var))
  {
    if(coloring.var!="type")
      warning("\"coloring.var\" must be either \"NULL\" or \"type\"")
  }
  

  
 
  dim1 <- axes[1]
  dim2 <- axes[2]
  
  lab.x <- paste("Dim ", axes[1], " (", signif(res.pca$eig[axes[1], 
                                                           2], 4), " %)", sep = "")
  lab.y <- paste("Dim ", axes[2], " (", signif(res.pca$eig[axes[2], 
                                                           2], 4), " %)", sep = "")
  
  
  if (choice == "ind") {
    if (is.null(main)) 
      main <- "Individuals component map"
    
    coord.ind<-res.pca$ind$coord
    
    if (is.null(xlim)) 
    {
      xmin <- min(coord.ind[, dim1])
      xmax <- max(coord.ind[, dim1])
      xlim <- c(xmin, xmax) * 1.2
    }
    if (is.null(ylim)) 
    {
      ymin <- min(coord.ind[, dim2])
      ymax <- max(coord.ind[, dim2])
      ylim <- c(ymin, ymax) * 1.2
    }
    
    if(is.null(col.ind) | is.null(coloring.ind))
    {
      col.plot.ind<-rep("black",nrow(coord.ind))
    }
    
    if (is.factor(coloring.ind))
    { 
      quali<-coloring.ind
      if (!is.null(col.ind))
      { 
        levels(quali)<-col.ind
        col.plot.ind<-quali
      }
      if(is.null(col.ind))
        col.plot.ind<-as.numeric(quali)
    }
    
    col.plot.ind.total<-col.plot.ind
    
    if(lim.cos2.plot == 0 & lim.contrib.plot==0)
    {
      lim.plot<-0
      select.ind<-1:nrow(coord.ind)
    }
    
    if(lim.cos2.plot != 0 & lim.contrib.plot==0)
    {
      lim.plot <- lim.cos2.plot
      base.lim <- res.pca$ind$cos2[,axes]
      select.ind <- which(apply(base.lim[,],1,sum)>=lim.plot)    
    }
    
    if(lim.cos2.plot == 0 & lim.contrib.plot!=0)
    {
      lim.plot<-lim.contrib.plot
      base.lim<-res.pca$ind$contrib[,axes]
      base.lim<-100*(base.lim/sum(eig.axes))
      select.ind<-which(apply(base.lim[,],1,sum)>=lim.plot)
    }
    
    if(length(select.ind)==0)
      warning("\"lim.cos.plot\" (or \"lim.contrib.plot\") is too large. No individuals can be plotted")
    
    coord.ind<-coord.ind[select.ind, , drop=FALSE]
    col.plot.ind<-col.plot.ind[select.ind]
    
    
    plot(coord.ind[, axes], xlim = xlim, ylim = ylim, xlab = lab.x, 
         ylab = lab.y, pch = 20, col = as.character(col.plot.ind), 
         cex = cex, main=main, ...)
    abline(h = 0, lty = 2, cex = cex)
    abline(v = 0, lty = 2, cex = cex)
    
    if(length(select.ind)!=0)
    {
      if(leg==T & is.factor(coloring.ind))
        legend(posleg, legend =paste(cl["coloring.ind"],levels(coloring.ind),sep="="), text.col = levels(as.factor(col.plot.ind.total)), 
               cex =cex.leg)
      
      if (label) 
        text(coord.ind[, axes], labels = rownames(coord.ind), 
             pos = 3, col = as.character(col.plot.ind), cex = cex, 
             ...)
    }
  }
  if (choice == "sqload") {
    if (is.null(main))
      main<-"Squared loadings"
    if (is.null(xlim)) 
    {
      xmax <- max(res.pca$sqload[, dim1])
      xlim <- c(-0.1, xmax * 1.2)
    }
    if (is.null(ylim)) 
    {
      ymax <- max(res.pca$sqload[, dim2])
      ylim <- c(-0.1, ymax * 1.2)
    }
    plot(0, 0, type = "n", xlab = lab.x, ylab = lab.y, xlim = xlim, 
         ylim = ylim, cex = cex,main=main, ...)
    abline(v = 0, lty = 2, cex = cex)
    abline(h = 0, lty = 2, cex = cex)
    
    
    if(is.null(coloring.var))
    {
      for (j in 1:nrow(res.pca$sqload)) 
      {
        arrows(0, 0, res.pca$sqload[j, dim1], res.pca$sqload[j, dim2], 
               length = 0.1, angle = 15, code = 2, cex = cex,...)
        if (label) 
        {
          if (res.pca$sqload[j, dim1] > res.pca$sqload[j, dim2]) 
          {
            pos <- 4
          }
          else pos <- 3
          text(res.pca$sqload[j, dim1], res.pca$sqload[j, dim2], labels = rownames(res.pca$sqload)[j], 
               pos = pos, cex = cex, ...)
        }
      }
      
    }
    if (!is.null(coloring.var))
    {
      if (coloring.var=="type")
      {
        for (j in 1:nrow(res.pca$sqload)) 
        {
          col.sq<-rep(c("blue","red"),c(p1,p2))
          arrows(0, 0, res.pca$sqload[j, dim1], res.pca$sqload[j, dim2], 
                 length = 0.1, angle = 15, code = 2, cex = cex, col=col.sq[j], 
                 ...)
          if (label) 
          {
            if (res.pca$sqload[j, dim1] > res.pca$sqload[j, dim2]) {
              pos <- 4
            }
            else pos <- 3
            text(res.pca$sqload[j, dim1], res.pca$sqload[j, dim2], labels = rownames(res.pca$sqload)[j], 
                 pos = pos, cex = cex, col=col.sq[j], ...)
          }
        }
        if (leg==TRUE)
          legend(posleg, legend = c("numerical","categorical"), text.col = c("blue","red"), 
                 cex = cex.leg)
      }
    }
    
  }
  if (choice == "levels") {
    if (is.null(main)) 
      main <- "Levels component map"
    
    if (lim.cos2.plot == 0 & lim.contrib.plot==0)
    {
      lim.plot<-0
      base.lim<-res.pca$levels$cos2[,axes]
    }
    
    if (lim.cos2.plot != 0 & lim.contrib.plot==0)
    {
      lim.plot<-lim.cos2.plot
      base.lim<-res.pca$levels$cos2[,axes]
    }
    
    if (lim.cos2.plot == 0 & lim.contrib.plot!=0)
    {
      lim.plot<-lim.contrib.plot
      base.lim<-res.pca$levels$contrib[,axes]
      base.lim<-100*(base.lim/sum(eig.axes))    
    }
    
    color<-rep(1,m)
    
    if (is.null(xlim)) 
    {
      xmin <- min(res.pca$levels$coord[, dim1])
      xmax <- max(res.pca$levels$coord[, dim1])
      xlim <- c(xmin, xmax) * 1.2
    }
    if (is.null(ylim)) 
    {
      ymin <- min(res.pca$levels$coord[, dim2])
      ymax <- max(res.pca$levels$coord[, dim2])
      ylim <- c(ymin, ymax) * 1.2
    }
    
    
    plot(0,0, xlim = xlim, ylim = ylim,
         xlab = lab.x, ylab = lab.y, type="n", cex = cex,main=main, ...)
    abline(h = 0, lty = 2, cex = cex)
    abline(v = 0, lty = 2, cex = cex)
    nrow.coord.lev <- 0
    if (!is.null(res.pca$levels$coord) ) 
    {
      coord.lev <- res.pca$levels$coord[, axes, drop = FALSE]
      nrow.coord.lev <- nrow(coord.lev)
      
      test.empty.plot<-c()
      for (v in 1:nrow(coord.lev)) 
      {
        if (sum(base.lim[v, ], na.rm = TRUE) >= lim.plot && !is.na(sum(base.lim[v, ], na.rm = TRUE))) {
          test.empty.plot<-c(test.empty.plot,1)
          points(coord.lev[v, 1], coord.lev[v,2], col = color[v],pch=20,cex = cex,...)
          
          if (label) 
          {
            if (abs(coord.lev[v, 1]) > abs(coord.lev[v,2])) 
            {
              if (coord.lev[v, 1] >= 0) 
                pos <- 4
              else pos <- 2
            }
            else {
              if (coord.lev[v, 2] >= 0) 
                pos <- 3
              else pos <- 1
            }
            text(coord.lev[v, 1], y = coord.lev[v, 2], 
                 labels = rownames(coord.lev)[v], pos = pos, 
                 col = color[v], cex = cex)
          }
        }
      }
      if(is.null(test.empty.plot)){
        warning("\"lim.cos.plot\" (or \"lim.contrib.plot\") is too large. No level can be plotted")
        return()
    }
    }
  }
  if (choice == "cor") {
    if (is.null(main)) 
      main <- "Correlation circle"
    if(lim.cos2.plot == 0 & lim.contrib.plot==0)
    {
      lim.plot<-0
      base.lim<-res.pca$quanti$cos2[,axes]
    }
    
    if(lim.cos2.plot != 0 & lim.contrib.plot==0)
    {
      lim.plot<-lim.cos2.plot
      base.lim<-res.pca$quanti$cos2[,axes]
    }
    
    if(lim.cos2.plot == 0 & lim.contrib.plot!=0)
    {
      lim.plot<-lim.contrib.plot
      base.lim<-res.pca$quanti$contrib[,axes]
      base.lim<-100*(base.lim/sum(eig.axes))     
    }
    
    
    if(is.null(xlim))
    {
      xlim = c(-1.1, 1.1)
    }
    if(is.null(ylim))
    {
      ylim = c(-1.1, 1.1)
    }
    
    col<-rep(1,p1)
    
    plot(0, 0, main = main, xlab = lab.x, ylab = lab.y, 
         xlim = xlim, ylim = ylim, col = "white", 
         asp = 1, cex = cex,...)
    x.cercle <- seq(-1, 1, by = 0.01)
    y.cercle <- sqrt(1 - x.cercle^2)
    lines(x.cercle, y = y.cercle)
    lines(x.cercle, y = -y.cercle)
    abline(v = 0, lty = 2, cex = cex)
    abline(h = 0, lty = 2, cex = cex)
    
    nrow.coord.var <- 0
    if (!is.null(res.pca["quanti"]$quanti$coord) ) 
    {
      coord.var <- res.pca$quanti$coord[, axes, drop = FALSE]
      nrow.coord.var <- nrow(coord.var)
      
      test.empty.plot<-c()      
      for (v in 1:nrow(coord.var)) 
      {
        if (sum(base.lim[v, ] , na.rm = TRUE) >= lim.plot && !is.na(sum(base.lim[v, ], na.rm = TRUE))) {
          test.empty.plot<-c(test.empty.plot,1)
          arrows(0, 0, coord.var[v, 1], coord.var[v,2], length = 0.1, angle = 15, code = 2, col = col[v],cex = cex)
          
          if (label) 
          {
            if (abs(coord.var[v, 1]) > abs(coord.var[v, 
                                                     2])) 
            {
              if (coord.var[v, 1] >= 0) 
                pos <- 4
              else pos <- 2
            }
            else {
              if (coord.var[v, 2] >= 0) 
                pos <- 3
              else pos <- 1
            }
            text(coord.var[v, 1], y = coord.var[v, 2], 
                 labels = rownames(coord.var)[v], pos = pos, 
                 col = col[v], cex = cex)
          }
        }
      }
      if(is.null(test.empty.plot)){
        warning("\"lim.cos.plot\" (or \"lim.contrib.plot\") is too large. No variable can be plotted")
        return()
     }
   }
 }
}


