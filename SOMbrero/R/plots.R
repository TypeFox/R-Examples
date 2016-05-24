### subfunctions
orderIndexes <- function(the.grid, type) {
  if(type=="3d") tmp <- 1:the.grid$dim[1]
  else tmp <- the.grid$dim[1]:1
  match(paste(rep(1:the.grid$dim[2],the.grid$dim[1]),"-",
              rep(tmp,
                  rep(the.grid$dim[2],the.grid$dim[1])),
              sep=""),
        paste(the.grid$coord[,1],"-",the.grid$coord[,2],
              sep=""))
}

averageByCluster <- function(x,clustering,the.grid) {
  mean.var <- matrix(NA,nrow=prod(the.grid$dim),ncol=ncol(x))
  ne.neurons <- which(as.character(1:prod(the.grid$dim))%in%
                        names(table(clustering)))
  mean.var[ne.neurons,] <- matrix(unlist(by(x,clustering,colMeans)),
                                  byrow=TRUE,ncol=ncol(x))
  colnames(mean.var) <- colnames(x)
  return(mean.var)
}

words2Freq <- function(words, clustering, the.grid, type) {
  if (type=="names") {
    freq.words <- matrix(0, ncol=length(unique(words)), nrow=prod(the.grid$dim))
    all.tables <- by(words, clustering, table)
    freq.words[as.numeric(sort(unique(clustering))),] <- matrix(unlist(
      all.tables), ncol=length(unique(words)), byrow=TRUE)
    colnames(freq.words) <- names(all.tables[[1]]) 
  } else if (type=="words") {
    freq.words <- matrix(0, ncol=ncol(words), nrow=prod(the.grid$dim))
    freq.words[as.numeric(sort(unique(clustering))),] <- apply(words,2,
                                                               function(word) {
      by(word,clustering,sum)
    })
    colnames(freq.words) <- colnames(words)
  }
  return(freq.words)
}

paramGraph <- function(the.grid, print.title, type) {
  if (print.title) {
    if(type%in%c("lines","pie","boxplot","names","words")) {
      the.mar <- c(0,0,1,0)
    } else the.mar <- c(2,1,1,1)
  } else {
    if(type%in%c("lines","pie","boxplot","names","words")) {
      the.mar <- c(0,0,0,0)
    } else the.mar <- c(2,1,0.5,1)
  }
  list("mfrow"=c(the.grid$dim[1], the.grid$dim[2]),"oma"=c(0,0,3,0), 
       "bty"="c", "mar"=the.mar)
}

myTitle <- function(args, what) {
  if (is.null(args$main)) {
    the.title <- switch(what,
                        "prototypes"="Prototypes overview",
                        "obs"="Observations overview",
                        "add"="Additional variable overview")
  } else the.title <- args$main
  return(the.title)
}

plotOneVariable <- function(ind, var.val, type, the.title, args) {
  args$ylab <- ""
  args$xlab <- ""
  if (type!="names" && type!="words" && !is.null(args$col) 
      && length(args$col)>1) {
    args$col <- args$col[ind]
    if ((type=="boxplot")&(!is.null(args$border)))
      args$border <- args$border[ind]
  }
  if (!is.null(the.title)) {
    args$main <- the.title
  } else args$main <- NULL
  if (all(is.na(var.val))) {
    plot(1,type="n",axes=F,xlab="",ylab="",main=args$main)
    if (type %in% c("lines","boxplot", "names", "words")) box()
  } else {
    if(type=="lines") {
      if (is.null(args$lwd)) args$lwd <- 2
      if (is.null(args$type)) args$type <- "l"
      if (is.null(args$col)) args$col <- "tomato"
      args$xaxt <- "n"
      args$yaxt <- "n"
      args$x <- var.val
      do.call("plot",args)
    } else if (type=="barplot") {
      args$height <- var.val
      args$axes <- FALSE
      if (is.null(args$col)) args$col <- "orange"
      if (is.null(args$names.arg)) args$names.arg <- ""
      if (is.null(args$border)) args$border <- FALSE
      do.call("barplot",args)
      abline(h=0, col=args$col)
    } else if (type=="boxplot") {
      args$x <- var.val
      args$axes <- FALSE
      if (is.null(args$names)) args$names <- FALSE
      if (is.null(args$col)) args$col <- "tan2"
      do.call("boxplot", args)
      box()
    } else if (type %in% c("names", "words")) {
      if (sum(var.val)>0) {
        if (is.null(args$min.freq)) args$min.freq <- 1
        args$words <- names(var.val)[var.val>0]
        args$freq <- sqrt(var.val[var.val>0])
        args$scale <- args$scale*max(var.val)
        if (is.null(args$ordered.colors)) args$ordered.colors <- TRUE
        if (is.null(args$rot.per)) args$rot.per <- 0
        do.call("wordcloud", args)
      } else plot(1,type="n",axes=F,xlab="",ylab="")
      box()
      title(the.title)
    }
  }
}

plotAllVariables <- function(what, type, values, clustering=NULL, print.title,
                             the.titles, is.scaled, the.grid, args) {
  par(paramGraph(the.grid, print.title, type))
  ordered.index <- orderIndexes(the.grid, type)
  if (!is.null(args$col) && length(args$col)>1 && 
        length(args$col)!=length(ordered.index)) {
    warning("unadequate number of colors; first color will be used for all\n", 
            call.=TRUE, immediate.=TRUE)
    args$col <- args$col[1]
  }
  if (print.title) {
    the.titles <- the.titles
  } else the.titles <- rep(NULL,length(ordered.index))
  if (type %in% c("names", "words")) {
    freq.words <- words2Freq(values, clustering, the.grid, type)
    if (is.null(args$colors)) {
      nb.breaks <- 6
      all.colors <- brewer.pal(9,"Purples")[5:9]
    } else if (length(args$colors)==1) {
      words.col <- matrix(args$colors, ncol=ncol(freq.words),
                          nrow=nrow(freq.words))
    } else {
      nb.breaks <- length(args$colors)+1
      all.colors <- args$colors
    }
    if (is.null(args$colors)|length(args$colors)>1) {
      the.breaks <- seq(min(freq.words[freq.words>0])-0.1,
                        max(freq.words)+0.1,length=nb.breaks)
      words.cut <- apply(freq.words,2,cut,breaks=the.breaks,labels=FALSE)
      words.col <- apply(words.cut, 2, function(wc) all.colors[wc])
    }
    if (is.null(args$scale)) {
      args$scale <- c(2,0.5)/max(freq.words)
    } else args$scale <- args$scale/max(freq.words)
    sapply(ordered.index, function(ind) {
      cur.args <- args
      if (sum(freq.words[ind,])>0) {
        cur.args$colors <- words.col[ind,freq.words[ind,]>0]
      }
      plotOneVariable(ind, freq.words[ind,], type, the.titles[ind], cur.args)
    })
  } else {
    is.scaled.values <- scale(values,is.scaled, is.scaled)
    args$ylim <- range(is.scaled.values)
    if (what=="prototypes"|type=="boxplot") {
      mean.val <- is.scaled.values
    } else mean.val <- averageByCluster(is.scaled.values, clustering, the.grid)
    if (type=="boxplot") {
      sapply(ordered.index, function(ind) {
        plotOneVariable(ind,matrix(as.matrix(mean.val)[which(clustering==ind),],
                                   ncol=ncol(as.matrix(values))), type, 
                        the.titles[ind], args)
      })
    } else {
      sapply(ordered.index, function(ind) {
        plotOneVariable(ind,mean.val[ind,], type, the.titles[ind], args)
      })
    }
  }
  title(main=myTitle(args, what), outer=TRUE)
  par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(5, 4, 4, 2)+0.1)
}

plotColor <- function(what, values, clustering, the.grid, my.palette, 
                      print.title, the.titles, args) {
  if (what!="prototypes") {
    x <- rep(NA, prod(the.grid$dim))
    ne.neurons <- which(as.character(1:prod(the.grid$dim))%in%
                          names(table(clustering)))
    x[ne.neurons] <- tapply(values, clustering, mean)
  } else x <- values
  if (is.null(my.palette)) {
    nb.breaks <- min(c(floor(1 + (3.3*log(prod(the.grid$dim),base=10))),
                       9))
  } else nb.breaks <- length(my.palette)
  the.breaks <- seq(min(x, na.rm=TRUE)-10^(-5),
                    max(x, na.rm=TRUE)+10^(-5),
                    length=nb.breaks+1)
  if (is.null(my.palette)) {
    my.colors <- heat.colors(nb.breaks)[nb.breaks:1]
  } else my.colors <- my.palette
  vect.color <- my.colors[cut(x, the.breaks, labels=FALSE)]
  if (what!="prototypes") vect.color[is.na(vect.color)] <- "white"
  plot.args <- c(list(x=the.grid, neuron.col=vect.color),
                 args)
  do.call("plot.myGrid",plot.args)
  if (print.title) {
    text(x=the.grid$coord[,1], y=the.grid$coord[,2],
         labels=the.titles, cex=0.7)
  }
}

plotRadar <- function(x,the.grid,what,print.title,the.titles,args) {
  args$main <- myTitle(args, what)
  if (print.title) {
      args$labels <- the.titles
  } else args$labels <- ""
  if (is.null(args$key.labels)) args$key.labels <- colnames(x)
  args$x <- x
  args$nrow <- the.grid$dim[1]
  args$ncol <- the.grid$dim[2]
  args$locations <- the.grid$coord
  if (is.null(args$draw.segments)) args$draw.segments <- TRUE
  if (is.null(args$len)) args$len <- 0.4
  do.call("stars", args)
}

plot3d <- function(x, the.grid, type, variable, args) {
  args$x <- 1:the.grid$dim[2]
  args$y <- 1:the.grid$dim[1]
  args$z <- t(matrix(data=x[orderIndexes(the.grid, "3d"),variable], 
                ncol=the.grid$dim[2], 
                nrow=the.grid$dim[1], byrow=TRUE))
  if (is.null(args$theta)) args$theta <- -20
  if (is.null(args$phi)) args$phi <- 20
  if (is.null(args$expand)) args$expand <- 0.4
  if (is.null(args$xlab)) args$xlab <- colnames(the.grid$coord)[1]
  if (is.null(args$ylab)) args$ylab <- colnames(the.grid$coord)[2]
  if (is.null(args$zlab)) args$zlab <- colnames(x)[variable]
  if (is.null(args$lwd)) args$lwd <- 1.5
  if (is.null(args$col)) args$col <- "tomato"
  do.call("persp", args)
}

plotOnePolygon <- function(ind, values, the.grid, col) {
  ## FIX IT think about a more general way to implement it
  cur.coords <- the.grid$coord[ind,]
  if (the.grid$topo=="square") {
    neighbors <- match(paste(rep((cur.coords[1]-1):(cur.coords[1]+1),c(3,3,3)),
                             rep((cur.coords[2]-1):(cur.coords[2]+1),3)), 
                       paste(the.grid$coord[,1], the.grid$coord[,2]))
    neighbors <- neighbors[-5]
    poly.coord <- matrix(nrow=8, ncol=2)
    poly.coord[is.na(neighbors),] <- tcrossprod(rep(1,sum(is.na(neighbors))),
                                                cur.coords)
    active.indexes <- setdiff(1:8,which(is.na(neighbors)))
    poly.coord[active.indexes,] <- t(sapply(seq_along(active.indexes),
                                            function(ind2) {
      n.ind <- neighbors[active.indexes[ind2]]
      n.coords <- the.grid$coord[n.ind,]
      cur.coords + (n.coords-cur.coords)*values[ind2]
    }))
    poly.coord <- poly.coord[c(1,4,6,7,8,5,3,2),]
    polygon(poly.coord, col=col)
  } else {
    stop("Sorry: 'polygon' is still to be implemented for this dist.type or/and
          this topology.", call.=TRUE)
  }
}

plotPolygon <- function(values, clustering, the.grid, my.palette, args) {
  plot.args <- c(list(x=the.grid),args)
  do.call("plot.myGrid",plot.args)
  if (is.null(my.palette)) {
    nb.breaks <- min(c(floor(1 + (3.3*log(prod(the.grid$dim), base=10))),9))
  } else nb.breaks <- length(my.palette)
  the.breaks <- seq(0,max(table(clustering)),length=nb.breaks+1)
  if (is.null(my.palette)) {
    my.colors <- heat.colors(nb.breaks)[nb.breaks:1]
  } else my.colors <- my.palette
  freq.clust <- sapply(1:prod(the.grid$dim), function(ind) sum(clustering==ind))
  vect.color <- my.colors[cut(freq.clust, the.breaks, labels=FALSE)]
  vect.color[is.na(vect.color)] <- "white"
  invisible(sapply(1:prod(the.grid$dim), function(ind) {
    plotOnePolygon(ind, values[[ind]], the.grid, vect.color[ind])
  }))
}

### SOM algorithm graphics
plotPrototypes <- function(sommap, type, variable, my.palette, print.title,
                           the.titles, is.scaled, view, args) {
  ## types : 3d, lines, barplot, radar, color, smooth.dist, poly.dist, umatrix,
  # mds
  
  # default value for type="lines"
  if (!is.element(type,c("3d","lines","barplot","radar","color", "poly.dist",
                         "umatrix", "smooth.dist", "mds", "grid.dist"))) {
    warning("incorrect type replaced by 'lines'\n", call.=TRUE, 
            immediate.=TRUE)
    type <- "lines"
  }
  # relational control
  if (sommap$parameters$type=="relational" && type %in% c("color", "3d"))
    stop("prototypes/", type, " plot is not available for 'relational'\n", 
         call.=TRUE)
  
  if (type=="lines" || type=="barplot") {
    if (sommap$parameters$type=="korresp") {
      if (view=="r")
        tmp.proto <- sommap$prototypes[,(ncol(sommap$data)+1):
                                         ncol(sommap$prototypes)]
      else
        tmp.proto <- sommap$prototypes[,1:ncol(sommap$data)]
    } else
      tmp.proto <- sommap$prototypes
    plotAllVariables("prototypes",type,tmp.proto,
                     print.title=print.title,the.titles=the.titles,
                     is.scaled=is.scaled,
                     the.grid=sommap$parameters$the.grid,args=args)
  } else if (type=="radar") {
    if (sommap$parameters$type=="korresp") {
      if (view=="r")
        tmp.proto <- sommap$prototypes[,(ncol(sommap$data)+1):
                                         ncol(sommap$prototypes)]
      else
        tmp.proto <- sommap$prototypes[,1:ncol(sommap$data)]
    } else
      tmp.proto <- sommap$prototypes
    plotRadar(tmp.proto, sommap$parameters$the.grid, "prototypes",
              print.title, the.titles, args)
  } else if (type=="color") {
    if (length(variable)>1) {
      warning("length(variable)>1, only first element will be considered\n", 
              call.=TRUE, immediate.=TRUE)
      variable <- variable[1]
    }
    if (sommap$parameters$type=="korresp" & (is.numeric(variable))) {
      if (view=="r") 
        tmp.var <- variable+ncol(sommap$data)
      else tmp.var <- variable
    } else tmp.var <- variable
    plotColor("prototypes", sommap$prototypes[,tmp.var], sommap$clustering,
              sommap$parameters$the.grid, my.palette, print.title, the.titles,
              args)
  } else if (type=="3d") {
    if (length(variable)>1) {
      warning("length(variable)>1, only first element will be considered\n", 
              call.=TRUE, immediate.=TRUE)
      variable <- variable[1]
    }
    if (sommap$parameters$type=="korresp" & (is.numeric(variable))) {
      if (view=="r") 
        tmp.var <- variable+ncol(sommap$data)
      else 
        tmp.var <- variable
    } else
      tmp.var <- variable
    plot3d(sommap$prototypes, sommap$parameters$the.grid, type, tmp.var, args)
  } else if (type=="poly.dist") {
    values <- protoDist(sommap, "neighbors")

    if (sommap$parameters$type=="relational") {
      if (sum(unlist(values)<0)>0) {
        stop("Impossible to plot 'poly.dist'!", call.=TRUE)
      } else values <- lapply(values,sqrt)
    }

    maxi <- max(unlist(values))
    values <- lapply(values, function(x) 0.429*((maxi-x)/maxi+0.05))
    plotPolygon(values, sommap$clustering, sommap$parameters$the.grid,
                my.palette, args)
    if (print.title) {
      text(x=sommap$parameters$the.grid$coord[,1]-0.1,
           y=sommap$parameters$the.grid$coord[,2]+0.1,
           labels=the.titles, cex=0.7)
    }
  } else if (type=="umatrix" || type=="smooth.dist") {
    values <- protoDist(sommap, "neighbors")

    if (sommap$parameters$type=="relational") {
      if (sum(unlist(values)<0)>0) {
        stop("Impossible to plot 'smooth.dist'!", call.=TRUE)
      } else values <- lapply(values,sqrt)
    }
    values <- unlist(lapply(values,mean))
    if (type=="umatrix") {
      plotColor("prototypes", values, sommap$clustering, 
                sommap$parameters$the.grid, my.palette, print.title, the.titles,
                args)
    } else {
      args$x <- 1:sommap$parameters$the.grid$dim[2]
      args$y <- 1:sommap$parameters$the.grid$dim[1]
      args$z <- matrix(data=values, nrow=sommap$parameters$the.grid$dim[2], 
                       ncol=sommap$parameters$the.grid$dim[1], byrow=TRUE)
      if (is.null(args$color.palette)) args$color.palette <- cm.colors
      if (is.null(args$main)) args$main <- "Distances between prototypes"
      if (is.null(args$xlab)) args$xlab <- "x"
      if (is.null(args$ylab)) args$ylab <- "y"
      do.call("filled.contour", args)
    }
  } else if (type=="mds") {
    if (sommap$parameters$type=="relational") {
      the.distances <- protoDist(sommap, "complete")

      if (sum(the.distances<0)>0) {
        stop("Impossible to plot 'MDS'!", call.=TRUE)
      } else the.distances <- sqrt(the.distances)
    }
    else the.distances <- dist(sommap$prototypes)
    proj.coord <- cmdscale(the.distances, 2)
    args$x <- proj.coord[,1]
    args$y <- proj.coord[,2]
    if (is.null(args$pch)) args$type <- "n"
    if (is.null(args$xlab)) args$xlab <- "x"
    if (is.null(args$ylab)) args$ylab <- "y"
    if (is.null(args$main)) args$main <- "Prototypes visualization by MDS"
    if (is.null(args$labels)) {
      the.labels <- as.character(1:nrow(sommap$prototypes))
    } else {
      the.labels <- args$labels
      args$labels <- NULL
    }
    if (is.null(args$col)) args$col <- "black"
    if (is.null(args$cex)) args$cex <- 1
    do.call("plot", args)
    text(proj.coord,the.labels,cex=args$cex,col=args$col)
  } else if (type=="grid.dist") {
    if (sommap$parameters$type=="relational") {
      the.distances <- protoDist(sommap, "complete")

     if (sum(the.distances<0)>0) {
      stop("Impossible to plot 'grid.dist'!", call.=TRUE)
     } else {
       the.distances <- sqrt(the.distances)
       the.distances <- the.distances[lower.tri(the.distances)]
     }
    } else the.distances <- dist(sommap$prototypes)
    args$x <- as.vector(the.distances)
    args$y <- as.vector(dist(sommap$parameters$the.grid$coord))
    if (is.null(args$pch)) args$pch <- '+'
    if (is.null(args$xlab)) args$xlab <- "prototype distances"
    if (is.null(args$ylab)) args$ylab <- "grid distances"
    do.call("plot", args)
  } else stop("Sorry: this type is still to be implemented.", call.=TRUE)
}

plotObs <- function(sommap, type, variable, my.palette, print.title, the.titles,
                    is.scaled, view, args) {
  ## types : hitmap, lines, names, color, barplot, boxplot, radar
  
  # default value is type="hitmap"
  if (!is.element(type,c("hitmap", "lines", "names", "color", "radar",
                         "barplot", "boxplot"))) {
    warning("incorrect type replaced by 'hitmap'\n", call.=TRUE, 
            immediate.=TRUE)
    type <- "hitmap"
  }

  # korresp control
  if (sommap$parameters$type=="korresp" && !(type%in%c("hitmap", "names"))) {
    warning("korresp SOM: incorrect type replaced by 'hitmap'\n", call.=TRUE, 
            immediate.=TRUE)
    type <- "hitmap"
  }
  # relational control
  if (sommap$parameters$type=="relational" && !(type%in%c("hitmap", "names"))) {
    warning("relational SOM: incorrect type replaced by 'hitmap'\n", 
            call.=TRUE, immediate.=TRUE)
    type <- "hitmap"
  }
  
  if (type=="lines" || type=="barplot") {
    plotAllVariables("obs", type, sommap$data, sommap$clustering, 
                     print.title, the.titles, is.scaled,
                     sommap$parameters$the.grid, args)
  } else if (type=="color") {
    if (length(variable)>1) {
      warning("length(variable)>1, only first element will be considered\n", 
              call.=TRUE, immediate.=TRUE)
      variable <- variable[1]
    } 
    plotColor("obs", sommap$data[,variable], sommap$clustering,
              sommap$parameters$the.grid, my.palette, print.title, the.titles,
              args)
  } else if (type=="radar") {
    mean.var <- averageByCluster(sommap$data, sommap$clustering,
                                 sommap$parameters$the.grid)
    plotRadar(mean.var, sommap$parameters$the.grid, "obs", print.title,
              the.titles, args)
  } else if (type=="hitmap") {
    freq <- sapply(1:nrow(sommap$prototypes), function(ind) {
      length(which(sommap$clustering==ind))
    })
    freq <- freq/sum(freq)
    # basesize is 0.45 for the maximum frequence
    basesize <- 0.45*sqrt(freq)/max(sqrt(freq))
    
    if (is.null(args$col)) {
      my.colors <- rep("pink", nrow(sommap$prototypes))
    } else if (length(args$col)==1) {
      my.colors <- rep(args$col, nrow(sommap$prototypes))
    } else {
      if(length(args$col)==nrow(sommap$prototypes)){
        my.colors <- args$col
      } else {
        warning("unadequate number of colors default color will be used\n", 
                immediate.=TRUE, call.=TRUE)
        my.colors <- rep("pink", nrow(sommap$prototypes))
      }
    }
    plot.args <- c(list(x=sommap$parameters$the.grid), args)
    do.call("plot.myGrid",plot.args)
    invisible(sapply(1:nrow(sommap$prototypes), function(ind){
      xleft <- (sommap$parameters$the.grid$coord[ind,1]-basesize[ind])
      xright <- (sommap$parameters$the.grid$coord[ind,1]+basesize[ind])
      ybottom <- (sommap$parameters$the.grid$coord[ind,2]-basesize[ind])
      ytop <- (sommap$parameters$the.grid$coord[ind,2]+basesize[ind])
      rect(xleft,ybottom,xright,ytop, col=my.colors[ind], border=NA)
    }))
  } else if (type=="boxplot") {
    if (length(variable)>5) {
      stop("maximum number of variables for type='boxplot' exceeded\n", 
           call.=TRUE)
    }
    plotAllVariables("obs", type, sommap$data[,variable], sommap$clustering, 
                     print.title, the.titles, is.scaled,
                     sommap$parameters$the.grid, args)
  } else if (type=="names") {
    if (sommap$parameters$type=="korresp") {
      values <- names(sommap$clustering)
    } else {
      if (!is.null(rownames(sommap$data))) {
        values <- rownames(sommap$data)
      } else values <- 1:nrow(sommap$data)
    }
    plotAllVariables("obs", type, values, sommap$clustering, 
                     print.title, the.titles, is.scaled,
                     sommap$parameters$the.grid, args)
  }
}

plotEnergy <- function(sommap, args) {
  # possible only if some intermediate backups have been done
  if (is.null(sommap$backup)) {
    stop("no intermediate backups have been registered\n", call.=TRUE)
  } else {
    if (is.null(args$main))
      args$main <- "Energy evolution"
    if (is.null(args$ylab)) args$ylab <- "Energy"
    if (is.null(args$xlab)) args$xlab <- "Steps"
    if (is.null(args$type)) args$type <- "b"
    if (is.null(args$pch)) args$pch <- "+"
    args$x <- sommap$backup$steps
    args$y <- sommap$backup$energy
    do.call("plot",args)
  }
}

projectFactor <- function(the.graph, clustering, the.factor, pie.color=NULL) {

  if (!is.factor(the.factor)) the.factor <- as.factor(the.factor)
  vertex.pie <- lapply(split(the.factor, factor(clustering)), table)
  if (is.null(pie.color)) {
    pie.color <- list(c(brewer.pal(8,"Set2"),
                        brewer.pal(12,"Set3"))[1:nlevels(the.factor)])
  }
  return(list("vertex.pie"=vertex.pie, "vertex.pie.color"=pie.color))
}

plotProjGraph <- function(proj.graph, print.title=FALSE, the.titles=NULL, 
                          s.radius=1, pie.graph=FALSE, pie.variable=NULL, ...) {
  
  args <- list(...)
  args$x <- proj.graph
  args$edge.width <- E(proj.graph)$weight/max(E(proj.graph)$weight)*10
  if (is.null(s.radius)) s.radius <- 1
  args$vertex.size <- s.radius*20*sqrt(V(proj.graph)$size)/
    max(sqrt(V(proj.graph)$size))
  if (is.null(args$vertex.label) & !print.title) args$vertex.label <- NA
  if (print.title) {
    if (is.null(args$vertex.label)) {
      if (is.null(the.titles)) {
        args$vertex.label <- V(proj.graph)$name
      } else 
        args$vertex.label <- the.titles[as.numeric(V(proj.graph)$name)]
    }
  }
    
  if (args$vertex.shape!="pie") {
    if (is.null(args$vertex.color))
      args$vertex.color <- brewer.pal(12,"Set3")[4]
    if (is.null(args$vertex.frame.color))
      args$vertex.frame.color <- brewer.pal(12,"Set3")[4]
  }

  par(bg="white")
  do.call("plot.igraph", args)
}

plotAdd <- function(sommap, type, variable, proportional, my.palette,
                    print.title, the.titles, is.scaled, s.radius, pie.graph,
                    pie.variable, args) {
  ## types : pie, color, lines, boxplot, names, words, graph, barplot, radar
  # to be implemented: graph
  
  # default value is type="pie"
  if (!is.element(type,c("pie", "color", "lines", "barplot", "words",
                         "boxplot", "names", "radar", "graph"))) {
    warning("incorrect type replaced by 'pie'\n", call.=TRUE, 
            immediate.=TRUE)
    type <- "pie"
  }
  if (is.null(variable)) {
    stop("for what='add', the argument 'variable' must be supplied\n", 
         call.=TRUE)
  }
  
  # korresp control
  if (sommap$parameters$type=="korresp") 
    stop("graphics of type 'add' do not exist for 'korresp'\n", call.=TRUE)
  
  if(type!="graph" && nrow(variable)!=nrow(sommap$data)){
    stop("length of additional variable does not fit length of the original
         data", call.=TRUE)
  }
  
  # switch between different types
  if (type=="pie") {
    if (!is.factor(variable)) variable <- as.factor(variable)
    cluster.freq <- tapply(variable,sommap$clustering,table)
    cluster.size <- table(sommap$clustering)
    par(paramGraph(sommap$parameters$the.grid, print.title, "pie"))
    ordered.index <- orderIndexes(sommap$parameters$the.grid, type)
    for (ind in ordered.index) {
      cur.cluster.vect.full <- cluster.freq[[as.character(ind)]]
      if (!is.null(cur.cluster.vect.full)) {
        cur.cluster.vect <- cur.cluster.vect.full[cur.cluster.vect.full>0]
        cur.args <- args
        cur.args$x <- cur.cluster.vect
        if (print.title) {
          cur.args$main <- the.titles[ind]
        } else cur.args$main <- NULL
        if (is.null(args$col)) {
          cur.args$col <- rainbow(nlevels(variable))[cur.cluster.vect.full>0]
        } else cur.args$col <- args$col[cur.cluster.vect.full>0]
        if (is.null(args$labels)) {
          cur.args$labels <- levels(variable)[cur.cluster.vect.full>0]
        } else cur.args$labels <- args$labels[cur.cluster.vect.full>0]
        if (proportional) {
          cur.args$radius=0.9*s.radius*sqrt(cluster.size[as.character(ind)]/
                                            max(cluster.size))
        } else {
          if (is.null(args$radius)) cur.args$radius <- 0.7
        }
        do.call("pie",cur.args)
      } else plot(1, type="n", bty="n", axes=FALSE)
    }
    if (is.null(args$main)) {
      title(main="Additional variable distribution", outer=TRUE)
    } else title(args$main, outer=TRUE)
    par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(5, 4, 4, 2)+0.1)
  } else if (type=="color") {
    if (!is.numeric(variable)) {
      stop("for type='color', argument 'variable' must be a numeric vector\n", 
           call.=TRUE)
    }
    plotColor("add", variable, sommap$clustering, sommap$parameters$the.grid, 
              my.palette, print.title, the.titles, args)
  } else if (type=="lines" || type=="barplot") {
    if (!all(apply(variable,2,is.numeric))) {
      stop("for type='lines' or 'barplot', argument 'variable' must be either a
           numeric matrix or a numeric data frame\n", call.=TRUE)
    }
    plotAllVariables("add", type, variable, sommap$clustering, 
                     print.title, the.titles, is.scaled,
                     sommap$parameters$the.grid, args)
  } else if (type=="radar") {
    if (ncol(variable)<2) {
      stop("for type='radar', argument 'variable' must have at least 2
           columns\n")
    }
    mean.var <- averageByCluster(variable, sommap$clustering,
                                 sommap$parameters$the.grid)
    plotRadar(mean.var, sommap$parameters$the.grid, "add", print.title,
              the.titles, args)
  } else if (type=="boxplot") {
    if(ncol(variable)>5) {
      stop("maximum number of variables (5) for type='boxplot' exceeded\n", 
           call.=TRUE)
    }
    plotAllVariables("add", type, variable, sommap$clustering, 
                     print.title, the.titles, is.scaled,
                     sommap$parameters$the.grid, args)
  } else if (type=="words") {
    if (is.null(colnames(variable))) {
      stop("no colnames for 'variable'", call.=TRUE)
    }
    plotAllVariables("add", type, variable, sommap$clustering, print.title,
                     the.titles, is.scaled, sommap$parameters$the.grid, args)
  } else if (type=="names") {
    if (ncol(variable) != 1) {
      stop("for type='names', argument 'variable' must be a vector", call.=TRUE)
    }
    plotAllVariables("add", type, as.character(variable), sommap$clustering, 
                     print.title, the.titles, is.scaled,
                     sommap$parameters$the.grid, args)
  } else if (type=="graph") {
    # controls
    if (!is.igraph(variable)){
      stop("for type='graph', argument 'variable' must be an igraph object\n", 
           call.=TRUE)
    }
    if (length(V(variable)) != nrow(sommap$data)){
      stop("length of additional variable does not fit length of the original
         data", call.=TRUE)
    }
    # case of pie
    if (pie.graph) {
      print("ok")
      if (is.null(pie.variable)) 
        stop("pie.graph is TRUE, you must supply argument 'pie.variable'\n", 
             call.=TRUE)
      
      if (nrow(as.matrix(pie.variable)) != nrow(sommap$data)) {
        stop("length of argument 'pie.variable' does not fit length of the 
             original data", call.=TRUE)
      }
      
      args$vertex.shape <- "pie"
      if (is.null(args$vertex.pie.color)) args$vertex.pie.color <- NULL
      proj.pie <- projectFactor(variable, sommap$clustering, pie.variable,
                                pie.color=args$vertex.pie.color)
      args$vertex.pie <- proj.pie$vertex.pie
      args$vertex.pie.color <- proj.pie$vertex.pie.color
    } else if (is.null(args$vertex.shape)) args$vertex.shape <- "circle"
    
    # create projected graph and plot
    args$proj.graph <- projectGraph(variable, sommap$clustering, 
                                    sommap$parameters$the.grid$coord)
    args$print.title <- print.title
    args$the.titles <- the.titles
    args$s.radius <- s.radius
    args$pie.graph <- pie.graph
    args$pie.variable <- pie.variable
    do.call("plotProjGraph", args)
  } else 
    stop("Sorry: this type is still to be implemented.", call.=TRUE)
}

plot.somRes <- function(x, what=c("obs", "prototypes", "energy", "add"), 
                        type=switch(what,
                                    "obs"="hitmap",
                                    "prototypes"="color",
                                    "add"="pie",
                                    "energy"=NULL),
                        variable = if (what=="add") NULL else 
                          if (type=="boxplot") 1:min(5,ncol(x$data)) else 1,
                        my.palette=NULL, 
                        is.scaled = if (x$parameters$type=="numeric") TRUE else
                          FALSE,
                        print.title=FALSE, the.titles=if (what!="energy") 
                          switch(type, 
                                 "graph"=1:prod(x$parameters$the.grid$dim),
                                 paste("Cluster",
                                       1:prod(x$parameters$the.grid$dim))),
                        proportional=TRUE, s.radius=1, pie.graph=FALSE, 
                        pie.variable=NULL,
                        view = if (x$parameters$type=="korresp") "r" else NULL,
                        ...) {
  args <- list(...)
  what <- match.arg(what)
  if ((x$parameters$type=="korresp")&&!(view%in%c("r","c")))
      stop("view must be one of 'r'/'c'",call.=TRUE)
  if (length(the.titles)!=prod(x$parameters$the.grid$dim) & what!="energy") {
    the.titles=switch(type,
                      "graph"=1:prod(x$parameters$the.grid$dim),
                      paste("Cluster",1:prod(x$parameters$the.grid$dim)))
    warning("unadequate length for 'the.titles'; replaced by default",
            call.=TRUE, immediate.=TRUE)
  }

  switch(what,
         "prototypes"=plotPrototypes(x, type, variable, my.palette, print.title,
                                     the.titles, is.scaled, view, args),
         "energy"=plotEnergy(x, args),
         "add"=plotAdd(x, type, if (type!="graph") as.matrix(variable) else 
           variable, proportional, my.palette, print.title, the.titles, 
                       is.scaled, s.radius, pie.graph, pie.variable, args),
         "obs"=plotObs(x, type, variable, my.palette, print.title, the.titles,
                       is.scaled, view, args))
}