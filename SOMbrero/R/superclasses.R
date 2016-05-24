############################### subfunctions ##################################
dendro3dProcess <- function(v.ind, ind, tree, coord, mat.moy, scatter) {
  if (tree$merge[ind,v.ind]<0) {
    res <- coord[abs(tree$merge[ind,v.ind]),]
    scatter$points3d(matrix(c(res,0,res,tree$height[ind]), 
                        ncol=3, byrow=TRUE), type="l")
  } else {
    res <- mat.moy[tree$merge[ind,v.ind],]
    scatter$points3d(matrix(c(res,tree$height[tree$merge[ind,v.ind]],
                          res,tree$height[ind]), ncol=3, 
                        byrow=TRUE), type="l")
  }
  return(res)
}

########################## super classes #####################################
superClass.somRes <- function(sommap, method="ward.D", members=NULL, k=NULL,
                              h=NULL, ...) {
  if (sommap$parameters$type=="relational") {
    the.distances <- protoDist(sommap, "complete")
    if (sum(the.distances<0)>0) {
      stop("Impossible to make super clustering!", call.=TRUE)
    } else the.distances <- as.dist(the.distances)
  }
  else the.distances <- as.dist(protoDist(sommap, "complete")^2)
  hc <- hclust(the.distances, method, members)
  if (!is.null(k) || !is.null(h)) {
    sc <- cutree(hc, k, h)
    res <- list("cluster"=as.numeric(sc), "tree"=hc, "som"=sommap)
  } else {
    res <- list("tree"=hc, "som"=sommap)
  }
  class(res) <- "somSC"
  return(res)
}

superClass <- function(sommap, method, members, k, h,...) {
  UseMethod("superClass")
}

################################ S3 functions #################################
print.somSC <- function(x, ...) {
  cat("\n   SOM Super Classes\n")
  cat("     Initial number of clusters : ", prod(x$som$parameters$the.grid$dim),
      "\n")
  if (length(x)>2) {
    cat("     Number of super clusters   : ", length(unique(x$cluster)), "\n\n")
  } else cat("     Number of super clusters not chosen yet.\n\n")
}

summary.somSC <- function(object, ...) {
  print(object)
  if (length(object)>2) {
    cat("\n  Frequency table")
    print(table(object$cluster))
    cat("\n  Clustering\n")
    output.clustering <- object$cluster
    names(output.clustering) <- seq_along(object$cluster)
    print(output.clustering)
    cat("\n")
  
    if (object$som$parameters$type=="numeric") {
      sc.clustering <- object$cluster[object$som$clustering]
      cat("\n  ANOVA\n")
      res.anova <- as.data.frame(t(sapply(1:ncol(object$som$data), function(ind) {
        res.aov <- summary(aov(object$som$data[,ind]~as.factor(sc.clustering)))
        c(round(res.aov[[1]][1,4],digits=3), round(res.aov[[1]][1,5],digits=8))
      })))
      names(res.anova) <- c("F", "pvalue")
      res.anova$significativity <- rep("",ncol(object$som$data))
      res.anova$significativity[res.anova$"pvalue"<0.05] <- "*"
      res.anova$significativity[res.anova$"pvalue"<0.01] <- "**"
      res.anova$significativity[res.anova$"pvalue"<0.001] <- "***"
      rownames(res.anova) <- colnames(object$som$data)
      
      cat("\n        Degrees of freedom : ", 
          summary(aov(object$som$data[,1]~as.factor(sc.clustering)))[[1]][1,1],
          "\n\n")
      print(res.anova)
      cat("\n")
    } else if (object$som$parameters$type=="relational") {
      if (object$som$parameters$scaling=="cosine") {
        norm.data <- preprocessData(object$som$data, object$parameters$scaling)
      } else norm.data <- object$som$data
      sse.total <- sum(norm.data)/(2*nrow(norm.data))
      
      sc.clustering <- object$cluster[object$som$clustering]
      
      sse.within <- sum(sapply(unique(sc.clustering), function(clust)
        sum(norm.data[sc.clustering==clust,sc.clustering==clust])/
          (2*sum(sc.clustering==clust))))
      
      n.clusters <- length(unique(sc.clustering))
      F.stat <- ((sse.total-sse.within)/sse.within) * 
        ((nrow(norm.data)-n.clusters)/(n.clusters-1))
      
      p.value <- 1-pf(F.stat, n.clusters-1, nrow(norm.data)-n.clusters)
      if (p.value<0.001) {
        sig <- "***"
      } else if (p.value<0.1) {
        sig <- "**"
      } else if (p.value<0.05) sig <- "*"
      
      cat("\n  ANOVA\n")
      cat("         F                       : ", F.stat, "\n")
      cat("         Degrees of freedom      : ", n.clusters-1, "\n")
      cat("         p-value                 : ", p.value, "\n")
      cat("                 significativity : ", sig, "\n")
    }
  }
}

plot.somSC <- function(x, type=c("dendrogram", "grid", "hitmap", "lines", 
                                 "barplot", "boxplot", "mds", "color", 
                                 "poly.dist", "pie", "graph", "dendro3d", 
                                 "radar", "projgraph"),
                       plot.var=TRUE, plot.legend=FALSE, add.type=FALSE, 
                       ...) {
  # TODO: add types "names" and "words"
  args <- list(...)
  type <- match.arg(type)

  if (type=="dendrogram") {
    args$x <- x$tree
    if (is.null(args$xlab)) args$xlab <- ""
    if (is.null(args$ylab)) args$ylab <- ""
    if (is.null(args$sub)) args$sub <- ""
    if (is.null(args$main)) args$main <- "Super-clusters dendrogram"
    if ((x$tree$method=="ward")&(plot.var)) {
      layout(matrix(c(2,2,1),ncol=3))
      Rsq <- cumsum(x$tree$height/sum(x$tree$height))
      plot(length(x$tree$height):1, Rsq, type="b", pch="+",
           xlab="Number of clusters", ylab="proportion of unexplained variance",
           main="Proportion of variance\n not explained by\n super-clusters")
      do.call("plot", args)
    } else do.call("plot", args)
    if (length(x)>2) {
      rect.hclust(x$tree, k=max(x$cluster))
    } else warning("Impossible to plot the rectangles: no super clusters.\n",
                   call.=TRUE, immediate.=TRUE)
    par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(5, 4, 4, 2)+0.1)
  } else  if (type=="dendro3d") {
    if (length(x)==3) {
      if ((!is.null(args$col))&(length(args$col)=max(x$cluster))) {
        clust.col.pal <- args$col
        clust.col <- args$col[x$cluster]
      } else {
        if (!is.null(args$col))
          warning("Incorrect number of colors
                  (does not fit the number of super-clusters);
                  using the default palette.\n", call.=TRUE, immediate.=TRUE)
        # create a color vector from RColorBrewer palette
        clust.col.pal <- brewer.pal(max(x$cluster), "Set2")
        clust.col <- clust.col.pal[x$cluster]
      }
    } else clust.col <- rep("black",prod(x$som$parameters$the.grid$dim))
    # FIX IT! maybe some more code improvements...  
    x.y.coord <- x$som$parameters$the.grid$coord+0.5
    if (floor(max(x$tree$height[-which.max(x$tree$height)]))==0) {
      z.max <- max(x$tree$height[-which.max(x$tree$height)])
    } else {
      z.max <- ceiling(max(x$tree$height[-which.max(x$tree$height)]))
    }
    spt <- scatterplot3d(x=x.y.coord[,1], y=x.y.coord[,2], 
                         z=rep(0,prod(x$som$parameters$the.grid$dim)), 
                         zlim=c(0, z.max), 
                         pch=19, color=clust.col, xlab="x", ylab="y",  
                         zlab="", x.ticklabs="", y.ticklabs="")
    horiz.ticks <- matrix(NA, nrow=prod(x$som$parameters$the.grid$dim)-1, ncol=2)
    for (neuron in 1:(prod(x$som$parameters$the.grid$dim)-1)) {
      vert.ticks <- sapply(1:2, dendro3dProcess, ind=neuron, tree=x$tree, 
                           coord=x.y.coord, mat.moy=horiz.ticks, scatter=spt)
      horiz.ticks[neuron,] <- rowMeans(vert.ticks)
      spt$points3d(matrix(c(vert.ticks[,1], x$tree$height[neuron],
                            vert.ticks[,2], x$tree$height[neuron]), ncol=3,
                          byrow=TRUE), type="l")
    }
  } else {
    if (length(x)<3) {
      stop("No super clusters: plot unavailable.\n")
    } else {
      if ((!is.null(args$col)) & (length(args$col)=max(x$cluster))) {
        clust.col.pal <- args$col
        clust.col <- args$col[x$cluster]
      } else {
        if (!is.null(args$col))
          warning("Incorrect number of colors
                  (does not fit the number of super-clusters);
                  using the default palette.\n", call.=TRUE, immediate.=TRUE)
        # create a color vector from RColorBrewer palette
        clust.col.pal <- brewer.pal(max(x$cluster), "Set2")
        clust.col <- clust.col.pal[x$cluster]
      }
      if (type=="grid") {
        if (plot.legend) {
          layout(matrix(c(2,2,1),ncol=3))
          plot.new()
          legend(x="center", legend=paste("Super cluster", 1:max(x$cluster)), 
                 col=clust.col.pal, pch=19)
        }
        args$x <- x$som$parameters$the.grid
        args$neuron.col <- clust.col
        do.call("plot.myGrid", args)
        par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(5, 4, 4, 2)+0.1)
      } else if (type %in% c("hitmap", "lines", "barplot", "boxplot", "mds",
                             "color", "poly.dist", "pie", "graph", "radar")) {
        if ((x$som$parameters$type=="korresp") && 
              (type %in% c("boxplot", "pie", "graph")))
            stop(type, " plot is not available for 'korresp' super classes\n", 
                 call.=TRUE)
        if ((x$som$parameters$type=="relational") && 
              (type %in% c("boxplot", "color")))
          stop(type, " plot is not available for 'relational' super classes\n", 
               call.=TRUE)
          
        if ((type%in%c("poly.dist", "radar"))&(plot.legend)) {
          plot.legend <- FALSE
          warning("Impossible to plot the legend with type '",type,"'.\n",
                  call.=TRUE, immediate.=TRUE)
        }
        if (!(type%in%c("graph","pie", "radar"))) {
          args$col <- clust.col
        } else if (type=="graph") {
          neclust <- which(!is.na(match(1:prod(x$som$parameters$the.grid$dim),
                                        unique(x$som$clustering))))
          if (is.null(args$pie.graph)) args$pie.graph <- FALSE
          if (!args$pie.graph) {
            args$vertex.color <- clust.col[neclust]
            args$vertex.frame.color <- clust.col[neclust]
          } else {
            if (plot.legend)
              warning("Impossible to plot the legend with type '",type,"'.\n",
                      call.=TRUE, immediate.=TRUE)
            plot.legend <- FALSE
            print.title <- TRUE
            args$vertex.label <- paste("SC",x$cluster[neclust])
            args$vertex.label.color <- "black"
          }
        }
        if (plot.legend) {
          layout(matrix(c(2,2,1),ncol=3))
          plot.new()
          legend(x="center", legend=paste("Super cluster", 1:max(x$cluster)), 
                 col=clust.col.pal, pch=19)
          if (type%in%c("lines","barplot","boxplot","color","pie", "poly.dist", 
                        "radar"))
            warning("Impossible to plot the legend with type '",type,"'.\n",
                    call.=TRUE, immediate.=TRUE)
        }
        args$x <- x$som
        if (!add.type) {
          if (type %in% c("hitmap", "boxplot")) {
            args$what <- "obs"
          } else if (type%in%c("graph", "pie")) {
            args$what <- "add"
          } else args$what <- "prototypes"
        } else args$what <- "add"
        args$type <- type
        if (type=="boxplot") args$border <- clust.col
        args$the.titles <- paste("SC",x$cluster)
        if (type%in%c("pie", "radar")) {
          args$print.title <- TRUE
        } else if (type%in%c("color", "poly.dist")) args$print.title <- FALSE
       do.call("plot.somRes", args)
       if (type=="color") 
         text(x=x$som$parameters$the.grid$coord[,1], 
              y=x$som$parameters$the.grid$coord[,2],
              labels=paste("SC",x$cluster))
       else if (type=="poly.dist")
         text(x=x$som$parameters$the.grid$coord[,1]-0.1, 
              y=x$som$parameters$the.grid$coord[,2]+0.1,
              labels=paste("SC",x$cluster) )
       par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(5, 4, 4, 2)+0.1)
      } else if (type=="projgraph") {
        # check arguments
        if (x$som$parameters$type=="korresp")
          stop(type, " plot is not available for 'korresp' super classes\n", 
               call.=TRUE)
        if (is.null(args$variable)) {
          stop("for type='projgraph', the argument 'variable' must be supplied (igraph object)\n", 
               call.=TRUE)
        }
        if (!is.igraph(args$variable)){
          stop("for type='projgraph', argument 'variable' must be an igraph object\n", 
               call.=TRUE)
        }
        if (length(V(args$variable)) != nrow(x$som$data)){
          stop("number of nodes in graph does not fit length of the original data", call.=TRUE)
        }
        
        # colors
        if ((!is.null(args$col)) & (length(args$col)=max(x$cluster))) {
          args$vertex.color <- args$col
          args$vertex.frame.color <- args$col
        } else {
          if (!is.null(args$col))
            warning("Incorrect number of colors 
  (does not fit the number of super-clusters);
  using the default palette.\n", call.=TRUE, immediate.=TRUE)
          # create a color vector from RColorBrewer palette
          args$vertex.color <- brewer.pal(max(x$cluster), "Set2")
          args$vertex.frame.color <- brewer.pal(max(x$cluster), "Set2")
        }
        
        if (plot.legend) {
          layout(matrix(c(2,2,1),ncol=3))
          plot.new()
          legend(x="center", legend=paste("Super cluster", 1:max(x$cluster)), 
                 col=args$vertex.color, pch=19)
        }
        
        # case of pie
        if (is.null(args$pie.graph)) args$pie.graph <- FALSE
        if (args$pie.graph) {
          if (is.null(args$pie.variable)) 
            stop("pie.graph is TRUE, you must supply argument 'pie.variable'\n", 
               call.=TRUE)
          if (nrow(as.matrix(args$pie.variable)) != nrow(x$som$data)) {
            stop("length of argument 'pie.variable' does not fit length of the 
             original data", call.=TRUE)
          }
          args$vertex.shape <- "pie"
          if (is.null(args$vertex.pie.color)) args$vertex.pie.color <- NULL
          proj.pie <- projectFactor(args$variable, x$cluster[x$som$clustering],
                                    args$pie.variable,
                                    pie.color=args$vertex.pie.color)
          args$vertex.pie <- proj.pie$vertex.pie
          args$vertex.pie.color <- proj.pie$vertex.pie.color
        } else if (is.null(args$vertex.shape)) args$vertex.shape <- "circle"
        
        # find projected graph
        proj.graph <- projectIGraph.somSC(x, args$variable)
        args$proj.graph <- proj.graph
        args$variable <- NULL
        do.call("plotProjGraph", args)
        
      } else stop("Sorry, this type is not implemented yet\n", call.=TRUE) 
    }
  }
}

projectIGraph.somSC <- function(object, init.graph, ...) {
  if (length(object) <= 2) 
    stop("The number of clusters has not been chosen yet. Cannot project the graph on super-clusters.\n", 
         call.=TRUE)
  if (object$som$parameters$type=="korresp")
    stop("projectIGraph is not available for 'korresp' super classes\n", 
         call.=TRUE)
  # project the graph into the SOM grid
  proj.graph <- projectIGraph.somRes(object$som, init.graph)
  # clustering of the non empty clusters
  induced.clustering <- object$cluster[as.numeric(V(proj.graph)$name)]
  # search for the positions (center of gravity) of the superclusters
  original.positions <- object$som$parameters$the.grid$coord
  positions <- cbind(tapply(original.positions[,1], object$cluster, mean),
                     tapply(original.positions[,2], object$cluster, mean))
  
  proj.graph.sc <- projectGraph(proj.graph, induced.clustering, positions)
  
  proj.graph.sc <- set.graph.attribute(proj.graph.sc, "layout", positions)
  return(proj.graph.sc)
}