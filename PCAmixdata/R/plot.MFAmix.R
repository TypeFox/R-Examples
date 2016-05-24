plot.MFAmix<-function (x, axes = c(1, 2), choice = "axes", label=TRUE, coloring.var = NULL, coloring.ind=NULL,
                       col.ind=NULL, col.groups=NULL,  partial = NULL,chrono=FALSE, lim.cos2.plot = 0, lim.contrib.plot=0, xlim = NULL,  ylim = NULL,
                       cex = 1, main = NULL, new.plot = FALSE,leg=TRUE,posleg="topleft",cex.leg=0.8, ...) 
{
  
  cl<-match.call()
  if (!inherits(x, "MFAmix")) 
    stop("use only with \"MFAmix\" objects")
  
  main<-main
  res.mfa <- x
  n<-nrow(res.mfa$ind$coord)
  eig.axes<-res.mfa$eig[axes,1]
  
  
  if(posleg=="topleft") posleg2<-"topright"
  if(posleg=="topright") posleg2<-"topleft"
  if(posleg=="bottomleft") posleg2<-"bottomright"
  if(posleg=="bottomright") posleg2<-"bottomleft"
  
  if (!(choice %in% c("ind", "sqload", "levels", "cor", "axes", "groups"))) 
    stop("\"choice\" must be either \"ind\",\"sqload\",\"cor\", \"levels\",\"axes\" or \"groups\"")
  
  if (lim.cos2.plot != 0 & lim.contrib.plot!=0)
    stop("use either \"lim.cos2.plot\" OR \"lim.contrib.plot\"")
  
  if(!is.null(partial)){
    if (partial!="all" && length(which(rownames(res.mfa$ind$coord) %in% partial))==0)
      stop("invalid value for \"partial\"")
  }
  
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
  
  if (!is.null(coloring.var)){
    if(coloring.var!="groups" &  coloring.var!="type")
      stop("'coloring.var' must be one of the following: 'NULL', 'type', 'groups'")
  }
  
  if (!is.null(coloring.var))
  {
    if(coloring.var!="type" & coloring.var!="groups")
      warning("\"coloring.var\" must be either \"NULL\", \"type\" or \"groups\"")
  }
  
  if (!is.null(coloring.var))
  {   
    if (coloring.var=="type")
    {
      if (choice=="ind" | choice=="cor" | choice=="levels"| choice=="axes")
      {
        warning("coloring.var=\"type\" is not used if choice=\"ind\", \"cor\",\"levels\" or \"axes\"")
        coloring.var <- NULL
      }
    } else     
      if (choice=="ind") {
        warning("\"coloring.var\" is not used if choice=\"ind\"")
        coloring.var <- NULL
      }
    
  }
  
  lab.x <- paste("Dim ", axes[1], " (", signif(res.mfa$eig[axes[1], 
                                                           2], 4), " %)", sep = "")
  lab.y <- paste("Dim ", axes[2], " (", signif(res.mfa$eig[axes[2], 
                                                           2], 4), " %)", sep = "")
  
  eig.axes<-res.mfa$eig[axes,1]
  
  group <- res.mfa$lst.groups
  nbre.grpe <- length(unique(group))
  
  if (!is.null(col.groups)){
    if (length(col.groups)!=nbre.grpe)
      warning("invalid value of \"col.groups\"")
  }
  
  p1 <- res.mfa$global.pca$rec$p1
  p <- res.mfa$global.pca$rec$p
  p2<-res.mfa$global.pca$rec$p2
  dim1 <- axes[1]
  dim2 <- axes[2]
  
  if (is.null(col.groups) | length(col.groups) != nbre.grpe) {
    col.groups <- 2:(nbre.grpe+1)
  }
  
  
  if (choice == "axes") {
    if (new.plot) 
      dev.new()
    if(is.null(xlim)){
      xlim = c(-1.1, 1.1)
    }
    if(is.null(ylim)){
      ylim = c(-1.1, 1.1)
    }
    if (is.null(main)) 
      main <- "Partial axes"
    coord.axes <- res.mfa$partial.axes$coord[, axes, drop = FALSE]
    plot(0, 0, xlab = lab.x, ylab = lab.y, xlim = xlim, ylim = ylim, col = "white", asp = 1, 
         cex = cex, main = main,...)
    x.cercle <- seq(-1, 1, by = 0.01)
    y.cercle <- sqrt(1 - x.cercle^2)
    lines(x.cercle, y = y.cercle)
    lines(x.cercle, y = -y.cercle)
    abline(v = 0, lty = 2, cex = cex)
    abline(h = 0, lty = 2, cex = cex)
    
    if(!is.null(coloring.var)){
      if (coloring.var == "groups") {
        col.var <- col.groups
        
        i = 1
        couleur.axes <- col.var[i]
        auxil = strsplit(rownames(res.mfa$partial.axes$coord)[1], 
                         ".", fixed = TRUE)[[1]]
        auxil2 = auxil[length(auxil)]
        for (j in 2:nrow(res.mfa$partial.axes$coord)) {
          auxil = strsplit(rownames(res.mfa$partial.axes$coord)[j], 
                           ".", fixed = TRUE)[[1]]
          if (auxil2 != auxil[length(auxil)]) {
            i = i + 1
            auxil2 = auxil[length(auxil)]
          }
          couleur.axes <- c(couleur.axes, col.var[i])
        }
      } 
    }
    
    
    if (is.null(coloring.var )) {
      couleur.axes <- NULL
      for (i in 1:nbre.grpe) couleur.axes <- c(couleur.axes, 
                                               rep("black", ncol(res.mfa$partial.axes$coord)))
    }
    
    for (v in 1:nrow(coord.axes)) {
      arrows(0, 0, coord.axes[v, 1], coord.axes[v, 2], 
             length = 0.1, angle = 15, code = 2, col = couleur.axes[v], 
             cex = cex)
      if (abs(coord.axes[v, 1]) > abs(coord.axes[v, 2])) {
        if (coord.axes[v, 1] >= 0) 
          pos <- 4
        else pos <- 2
      }
      else {
        if (coord.axes[v, 2] >= 0) 
          pos <- 3
        else pos <- 1
      }
      text(coord.axes[v, 1], y = coord.axes[v, 2], labels = rownames(coord.axes)[v], 
           pos = pos, col = couleur.axes[v], cex = cex)
    }
    
    if(!is.null(coloring.var)){
      if ((coloring.var == "groups") & (leg==T)) {
        legend(posleg, legend = rownames(res.mfa$groups$Lg), 
               text.col = unique(couleur.axes), cex = cex.leg)
      }
    }  
  }
  
  if (choice == "groups") {
    if (new.plot) 
      dev.new()
    if (is.null(main)) 
      main <- "Groups contributions"
    coord.actif <- res.mfa$groups$coord[, axes, drop = FALSE]
    if (is.null(xlim)){
      xlim <- c(0,1)}
    if (is.null(ylim)){
      ylim <- c(0,1)}  
    if (is.null(coloring.var)) {
      col.var = rep("darkred", nrow(coord.actif))
    }   
    if (!is.null(coloring.var)) {
      if (coloring.var == "groups") 
        col.var<-col.groups
    }
    
    plot(coord.actif, xlab = lab.x, ylab = lab.y, xlim = xlim, ylim = ylim, pch = 17, col = col.var[1:nrow(coord.actif)], 
         cex = cex, main = main, cex.main = cex * 1.2, asp = 1,...)
    if (label) 
      text(coord.actif[, 1], y = coord.actif[, 2], labels = rownames(coord.actif), 
           pos = 3, col = col.var[1:nrow(coord.actif)], 
           cex = cex)
    
  }
  
  
  if (choice=="sqload") {
    if (is.null(main)) 
      main <- "Squared loadings"
    
    if (is.null(xlim)){
      xmax <- max(res.mfa$global.pca$sqload[,axes[1]])
      xlim <- c(-0.1, xmax*1.2)}
    if (is.null(ylim)){
      ymax <- max(res.mfa$global.pca$sqload[,axes[2]])
      ylim <- c(-0.1, ymax*1.2)}
    
    if (!is.null(coloring.var)){
      if (coloring.var == "groups") {
        col.var<-col.groups
        color<-NULL
        for (i in 1:nbre.grpe){  
          if (!is.null(nrow(res.mfa$separate.analyses[[i]]$sqload))){
            color <- c(color,rep(col.var[i], nrow(res.mfa$separate.analyses[[i]]$sqload)))}
        }
        col.var<-color
      }
    }
    
    
    if (is.null(coloring.var)) {    
      color <- rep(1,nrow(res.mfa$sqload))
    }
    
    if (!is.null(coloring.var)){
      if(coloring.var=="type"){
        names.ind<-rownames(res.mfa$sqload)
        names.ind[which(names.ind %in% rownames(res.mfa$quali$contrib))]<-"red"
        names.ind[which(names.ind!="red")]<-"blue"
        color<-names.ind
      }
    }
    
    
    plot(0, 0, type="n",xlab = lab.x, ylab = lab.y, xlim = xlim, ylim = ylim,cex=cex,main=main,...)
    abline(v = 0, lty = 2,cex=cex)
    abline(h = 0, lty = 2,cex=cex)
    
    for (j in 1:nrow(res.mfa$global.pca$sqload)) {
      arrows(0, 0, res.mfa$global.pca$sqload[j,axes[1]], res.mfa$global.pca$sqload[j,axes[2]], length = 0.1, angle = 15, code = 2,cex=cex,col=color[j])
      
      if (res.mfa$global.pca$sqload[j,axes[1]] > res.mfa$global.pca$sqload[j,axes[2]]) { pos <- 4
      } else  pos <- 3
      if(label)
        text(res.mfa$global.pca$sqload[j,axes[1]],res.mfa$global.pca$sqload[j,axes[2]],labels = rownames(res.mfa$global.pca$sqload)[j], pos = pos,cex=cex,col=color[j])  
    }
    
    if (!is.null(coloring.var)){
      if ((coloring.var == "groups")& (leg==T)) {
        legend(posleg, legend = rownames(res.mfa$groups$Lg), 
               text.col = unique(color), cex = cex.leg)
      }
      if (coloring.var=="type" & (leg==TRUE)) {
        legend(posleg, legend = c("numerical","categorical"), text.col = c("blue","red"), cex=cex.leg)
      }
    }    
  }
  
  
  if (choice == "cor") {
    if (is.null(main)) 
      main <- "Correlation circle"
    if(lim.cos2.plot == 0 & lim.contrib.plot==0){
      lim.plot<-0
      base.lim<-res.mfa$quanti$cos2[,axes]
    }
    
    if(lim.cos2.plot != 0 & lim.contrib.plot==0){
      lim.plot<-lim.cos2.plot
      base.lim<-res.mfa$quanti$cos2[,axes]
    }
    
    if(lim.cos2.plot == 0 & lim.contrib.plot!=0){
      lim.plot<-lim.contrib.plot
      base.lim<-res.mfa$quanti$contrib[,axes]
      base.lim<-100*(base.lim/sum(eig.axes))  
      
    }
    
    
    if(is.null(xlim)){
      xlim = c(-1.1, 1.1)
    }
    if(is.null(ylim)){
      ylim = c(-1.1, 1.1)
    }
    if (new.plot) 
      dev.new()
    
    if (is.null(coloring.var)) {
      col.var <- rep(1,nrow(res.mfa$quanti$coord))
      col<-col.var
    }
    
    if (!is.null(coloring.var)){
      if (coloring.var == "groups") {
        col.var <- col.groups
        col<-NULL
        for (i in 1:nbre.grpe){  
          if (!is.null(nrow(res.mfa$separate.analyses[[i]]$quanti$coord))){
            col <- c(col,rep(col.var[i], nrow(res.mfa$separate.analyses[[i]]$quanti$coord)))}
        }
      }
    }
    
    plot(0, 0, main = main, xlab = lab.x, ylab = lab.y, 
         xlim = xlim, ylim = ylim, col = "white", 
         asp = 1, cex = cex,...)
    x.cercle <- seq(-1, 1, by = 0.01)
    y.cercle <- sqrt(1 - x.cercle^2)
    lines(x.cercle, y = y.cercle)
    lines(x.cercle, y = -y.cercle)
    abline(v = 0, lty = 2, cex = cex)
    abline(h = 0, lty = 2, cex = cex)
    if (!is.null(coloring.var)){
      if (coloring.var == "groups" & (leg==T)) 
        legend(posleg, legend = rownames(res.mfa$groups$Lg), text.col = col.groups, cex = cex.leg)
    }
    
    nrow.coord.var <- 0
    if (!is.null(res.mfa["quanti"]$quanti$coord) ) {
      coord.var <- res.mfa$quanti$coord[, axes, drop = FALSE]
      nrow.coord.var <- nrow(coord.var)
      
      test.empty.plot<-c()      
      for (v in 1:nrow(coord.var)) {
        if (sum(base.lim[v, ] , na.rm = TRUE) >= lim.plot && !is.na(sum(base.lim[v, ], na.rm = TRUE))) {
          test.empty.plot<-c(test.empty.plot,1)
          
          arrows(0, 0, coord.var[v, 1], coord.var[v,2], length = 0.1, angle = 15, code = 2, col = col[v],cex = cex)
          
          if (label) {
            if (abs(coord.var[v, 1]) > abs(coord.var[v, 
                                                     2])) {
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
  
  
  
  if (choice == "ind") {
    if (is.null(main)) 
      main <- "Individuals component map"
    if (new.plot) 
      dev.new()
    coord.ind<-res.mfa$ind$coord[,axes]
    
    if(is.null(col.ind) | is.null(coloring.ind)){
      col.plot.ind<-rep("black",nrow(coord.ind))
    }
    
    if (is.factor(coloring.ind)){ 
      quali<-coloring.ind
      if (!is.null(col.ind)){ 
        levels(quali)<-col.ind
        col.plot.ind<-quali
      }
      else{col.plot.ind<-as.numeric(quali)}
    }
    col.plot.ind.total<-col.plot.ind
    
    
    if(lim.cos2.plot == 0 & lim.contrib.plot==0){
      lim.plot<-0
      select.ind<-1:n
    }
    
    if(lim.cos2.plot != 0 & lim.contrib.plot==0){
      lim.plot<-lim.cos2.plot
      base.lim<-res.mfa$ind$cos2[,axes]
      select.ind<-which(apply(base.lim[,],1,sum)>=lim.plot)    
    }
    
    if(lim.cos2.plot == 0 & lim.contrib.plot!=0){
      lim.plot<-lim.contrib.plot
      base.lim<-res.mfa$ind$contrib[,axes]
      base.lim<-100*(base.lim/sum(eig.axes))
      select.ind<-which(apply(base.lim[,],1,sum)>=lim.plot)
    }
    if (length(select.ind)==0){
      plot(NULL, xlab = lab.x,ylab=lab.y,main=main,xlim=c(-5,5),ylim=c(-5,5),...)
      abline(h = 0, lty = 2, cex = cex)
      abline(v = 0, lty = 2, cex = cex)
      warning("\"lim.cos.plot\" (or \"lim.contrib.plot\") is too large. No individuals can be plotted")
      return()
    }

      
    coord.ind<-coord.ind[select.ind, , drop=FALSE]
    col.plot.ind<-col.plot.ind[select.ind]
    
    
    
    if (!is.null(partial) & chrono==FALSE){
      
      if (length(partial)==1){
        if (partial=="all") ind.partial<-1:nrow(coord.ind)
        if (partial!="all") ind.partial<-which(rownames(coord.ind) %in% partial)
        if (length(ind.partial)==0){
          plot(NULL, xlab = lab.x,ylab=lab.y,main=main,xlim=c(-5,5),ylim=c(-5,5),...)
          abline(h = 0, lty = 2, cex = cex)
          abline(v = 0, lty = 2, cex = cex)
          warning("\"lim.cos.plot\" (or \"lim.contrib.plot\") is too large for partial individuals selected.\n No partial individuals can be plotted")
          return()
        }
      }
      if (length(partial)>1){
        ind.partial<-which(rownames(coord.ind) %in% partial)
        if (length(ind.partial)==0){
          plot(NULL, xlab = lab.x,ylab=lab.y,main=main,xlim=c(-5,5),ylim=c(-5,5),...)
          abline(h = 0, lty = 2, cex = cex)
          abline(v = 0, lty = 2, cex = cex)
          warning("\"lim.cos.plot\" (or \"lim.contrib.plot\") is too large for partial individuals selected.\n No partial individuals can be plotted")
          return()
        }
      }
      
      if(is.null(xlim) & is.null(ylim)){
        xlim<-c(0,0)
        ylim<-c(0,0)
        for (i in 1:nbre.grpe){
          t<-res.mfa$ind.partial[[i]][[1]][select.ind,,drop=FALSE][ind.partial,axes,drop=FALSE]
          xlim<-c(min(xlim[1],min(t[,1])), max(xlim[2],max(t[,1])))
          ylim<-c(min(ylim[1],min(t[,2])), max(ylim[2],max(t[,2])))
        }
        xlim<-xlim*1.2
        ylim<-ylim*1.2
        
      }
      
      plot(coord.ind, xlim = xlim, ylim = ylim, xlab = lab.x, 
           ylab = lab.y, pch = 20, col =as.character(col.plot.ind), cex = cex,main=main,...)
      abline(h = 0, lty = 2, cex = cex)
      abline(v = 0, lty = 2, cex = cex)
      
      for (i in 1:nbre.grpe){
        bary.coord<-coord.ind[ind.partial, ,drop=FALSE]
        t<-res.mfa$ind.partial[[i]][[1]][select.ind,,drop=FALSE][ind.partial,axes,drop=FALSE]
        points(t,col=col.groups[i],pch=20)
        for (j in 1:length(ind.partial)){
          m<-list(x=c(bary.coord[j,1],t[j,1]), y=c(bary.coord[j,2],t[j,2]))
          lines(m,col=col.groups[i]) 
        }
      }
      
      if (leg==T) {
        legend(posleg, legend = rownames(res.mfa$groups$Lg), 
               text.col = col.groups, cex = cex.leg)
        
      }
      
    }
    
    if (!is.null(partial) & chrono==TRUE){
      
      if (length(partial)==1){
        if (partial=="all") ind.partial<-1:nrow(coord.ind)
        if (partial!="all") ind.partial<-which(rownames(coord.ind) %in% partial)
        if (length(ind.partial)==0){
          plot(NULL, xlab = lab.x,ylab=lab.y,main=main,xlim=c(-5,5),ylim=c(-5,5),...)
          abline(h = 0, lty = 2, cex = cex)
          abline(v = 0, lty = 2, cex = cex)
          warning("\"lim.cos.plot\" (or \"lim.contrib.plot\") is too large for partial individuals selected.\n No partial individuals can be plotted")
          return()
        }
      }
      if (length(partial)>1){
        ind.partial<-which(rownames(coord.ind) %in% partial)
        if (length(ind.partial)==0){
          plot(NULL, xlab = lab.x,ylab=lab.y,main=main,xlim=c(-5,5),ylim=c(-5,5),...)
          abline(h = 0, lty = 2, cex = cex)
          abline(v = 0, lty = 2, cex = cex)
          warning("\"lim.cos.plot\" (or \"lim.contrib.plot\") is too large for partial individuals selected.\n No partial individuals can be plotted")
          return()
        }
      }
      
      if(is.null(xlim) & is.null(ylim)){
        xlim<-c(0,0)
        ylim<-c(0,0)
        for (i in 1:nbre.grpe){
          t<-res.mfa$ind.partial[[i]][[1]][select.ind,,drop=FALSE][ind.partial,axes,drop=FALSE]
          xlim<-c(min(xlim[1],min(t[,1])), max(xlim[2],max(t[,1])))
          ylim<-c(min(ylim[1],min(t[,2])), max(ylim[2],max(t[,2])))
        }
        xlim<-xlim*1.2
        ylim<-ylim*1.2
        
      }
      
      plot(coord.ind, xlim = xlim, ylim = ylim, xlab = lab.x, 
           ylab = lab.y, pch = 20, col = as.character(col.plot.ind), cex = cex,main=main,...)
      abline(h = 0, lty = 2, cex = cex)
      abline(v = 0, lty = 2, cex = cex)
      
      list.coord.partial<-list()
      for (i in 1:nbre.grpe){
        bary.coord<-coord.ind[ind.partial, ,drop=FALSE]
        t<-res.mfa$ind.partial[[i]][[1]][select.ind,,drop=FALSE][ind.partial,axes,drop=FALSE]
        points(t,col=col.groups[i],pch=20)
        list.coord.partial[[i]]<-t
      }
      
      
      for (i in 1:(nbre.grpe-1)){
        for (j in 1:length(ind.partial)){
          x1<-list.coord.partial[[i]][j,1]
          x2<-list.coord.partial[[i+1]][j,1]
          y1<-list.coord.partial[[i]][j,2]
          y2<-list.coord.partial[[i+1]][j,2]
          lines(x=c(x1,x2),y=c(y1,y2))
        }
      }
      
      
      if (leg==T) {
        legend(posleg, legend = rownames(res.mfa$groups$Lg), 
               text.col = col.groups, cex = cex.leg)
        
      }
      
    }
    
    if(is.null(partial)){
      if (is.null(xlim))xlim<-range(coord.ind[,1])*1.2
      if (is.null(ylim))ylim<-range(coord.ind[,2])*1.2
      
      plot(coord.ind, xlim = xlim, ylim = ylim, xlab = lab.x, 
           ylab = lab.y, pch = 20, col = as.character(col.plot.ind), cex = cex,main=main,...)
      abline(h = 0, lty = 2, cex = cex)
      abline(v = 0, lty = 2, cex = cex)
      
    }
    
    
    
    if(leg==TRUE & is.factor(coloring.ind) & is.null(partial)){
      legend(posleg, legend =paste(cl["coloring.ind"],levels(coloring.ind),sep="="), text.col = levels(as.factor(col.plot.ind.total)), cex = cex.leg)
    }
    
    if(leg==TRUE & is.factor(coloring.ind) & !is.null(partial)){
      legend(posleg2, legend =paste(cl["coloring.ind"],levels(coloring.ind),sep="="), text.col = levels(as.factor(col.plot.ind.total)), cex = cex.leg)
    }
    
    if (label) 
      text(coord.ind, labels = rownames(coord.ind), 
           pos = 3, col = as.character(col.plot.ind), cex = cex, 
           ...)
    
  }
  if (choice == "levels") {
    if (is.null(main)) 
      main <- "Levels component map"
    if (new.plot) 
      dev.new()
    
    if(lim.cos2.plot == 0 & lim.contrib.plot==0){
      lim.plot<-0
      base.lim<-res.mfa$levels$cos2[,axes]
    }
    
    if(lim.cos2.plot != 0 & lim.contrib.plot==0){
      lim.plot<-lim.cos2.plot
      base.lim<-res.mfa$levels$cos2[,axes]
    }
    
    if(lim.cos2.plot == 0 & lim.contrib.plot!=0){
      lim.plot<-lim.contrib.plot
      base.lim<-res.mfa$levels$contrib[,axes]
      base.lim<-100*(base.lim/sum(eig.axes))
      
      
    }
    if (!is.null(coloring.var)){
      if (coloring.var == "groups") {
        col.var<-col.groups
        color<-NULL
        for (i in 1:nbre.grpe){  
          if (!is.null(nrow(res.mfa$separate.analyses[[i]]$levels$coord))){
            color <- c(color,rep(col.var[i], nrow(res.mfa$separate.analyses[[i]]$levels$coord)))}
        }
      }
    }
    
    if (is.null(coloring.var)) {
      color<-rep(1,nrow(res.mfa$levels$coord))
    }
    
    if (is.null(xlim)) {
      xmin <- min(res.mfa$levels$coord[, dim1])
      xmax <- max(res.mfa$levels$coord[, dim1])
      xlim <- c(xmin, xmax) * 1.2
    }
    if (is.null(ylim)) {
      ymin <- min(res.mfa$levels$coord[, dim2])
      ymax <- max(res.mfa$levels$coord[, dim2])
      ylim <- c(ymin, ymax) * 1.2
    }
    
    
    plot(0,0, xlim = xlim, ylim = ylim,
         xlab = lab.x, ylab = lab.y, type="n", cex = cex,main=main, ...)
    abline(h = 0, lty = 2, cex = cex)
    abline(v = 0, lty = 2, cex = cex)
    nrow.coord.lev <- 0
    if (!is.null(res.mfa$levels$coord) ) {
      coord.lev <- res.mfa$levels$coord[, axes, drop = FALSE]
      nrow.coord.lev <- nrow(coord.lev)
      
      test.empty.plot<-c()
      for (v in 1:nrow(coord.lev)) {
        if (sum(base.lim[v, ], na.rm = TRUE) >= lim.plot && !is.na(sum(base.lim[v, ], na.rm = TRUE))) {
          test.empty.plot<-c(test.empty.plot,1)
          points(coord.lev[v, 1], coord.lev[v,2], col = color[v],pch=20,cex = cex,...)
          
          if (label) {
            if (abs(coord.lev[v, 1]) > abs(coord.lev[v,2])) {
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
    
    
    if (!is.null(coloring.var)){
      if ((coloring.var == "groups") & (leg==T)) {
        legend(posleg, legend = rownames(res.mfa$groups$Lg), 
               text.col = col.groups, cex = cex.leg)
      }
    }
  }
}
