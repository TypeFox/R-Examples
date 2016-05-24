is.simone <- function(x){
  if (class(x)=="simone") TRUE else FALSE
}

plot.simone <- function(x, output=c("BIC","AIC","ROC","PR","path.edges",
                             "path.penalty","sequence"), ref.graph=NULL,
                        ask=TRUE, ...) {
  
  if (!is.simone(x))
    stop("must be a simone object")

  if ("BIC" %in% output) {
    plot.BIC(x, ...)
    if (ask)
      readline(prompt="\nPress return for next plot...")
  }

  if ("AIC" %in% output) {
    plot.AIC(x, ...)
    if (ask)
      readline(prompt="\nPress return for next plot...")
  }

  if ("ROC" %in% output & !is.null(ref.graph)) {
    plot.ROC(x, ref.graph, ...) 
    if (ask)
      readline(prompt="\nPress return for next plot...")
  }

  if ("PR" %in% output & !is.null(ref.graph)) {
    plot.PR(x, ref.graph, ...)
    if (ask)
      readline(prompt="\nPress return for next plot...")
  }

  if ("path.edges" %in% output) {
    plot.path(x, x.axis="degrees of freedom (edges)", ...)
    if (ask)
      readline(prompt="\nPress return for next plot...")
  }
  
  if ("path.penalty" %in% output) {
    plot.path(x, x.axis="penalty level", ...)
    if (ask)
      readline(prompt="\nPress return for next plot...")
  }

  if ("sequence" %in% output & !is.list(x$networks[[1]])) {
    plot.sequence(x)
  }  
}

plot.path <- function(o, x.axis = "degrees of freedom (edges)",
                      main="regularization path", ...) {

  ## don't the warning message about pch recycling
  old.warn <- getOption("warn")
  options(warn=-1)

  directed <- o$control$edges.sym.rule == "NO"

  T <- 1
  ## Collect all the graphs in the same list
  if (is.list(o$networks[[1]])) {
    T <- length(o$networks[[1]])
  }
  
  layout(matrix(1:T,floor(sqrt(T)),ceiling(sqrt(T))))
  
  for (t in 1:T) {
    y <- c()
    for (k in 1:length(o$networks)) {
      
      if (T == 1) {
        theta  <- o$networks[[k]]
        main.t <- main
      } else {
        main.t <- paste(main," network #",t)
        theta  <- o$networks[[k]][[t]]
      }
      
      if (!directed) {
        y <- rbind(y,as.vector(theta[upper.tri(theta)]))
      } else {
        y <- rbind(y,as.vector(theta))
      }
    }
    
    x  <- switch(x.axis, "degrees of freedom (edges)" = o$n.edges[,t],
                 "penalty level" = o$penalties)
    
    matplot(x, y, xlab = x.axis, type = "l", ylab = "coefficients", ...)
    title(main.t,line=2.5)
    abline(h = 0, lty = 3)
    
    abline(v=x[which.max(o$AIC)], lty=3, col="red")
    axis(1, at = x[which.max(o$AIC)], labels = "AIC", las=3, cex.axis=0.75)
    
    abline(v=x[which.max(o$BIC)], lty=3, col="blue")
    axis(3, at = x[which.max(o$BIC)], labels = "BIC", las=3, cex.axis=0.75)

  }
  
  options(warn=old.warn)
}

plot.BIC <- function(o, main="Bayesian Information Criterion", ...) {

  plot(o$penalties,o$BIC, main=main, type="b", xlab="penalty level", ...)
  abline(v=o$penalties[which.max(o$BIC)], lty=3, col="red")  
}

plot.AIC <- function(o, main="Akaike Information Criterion", ...) {

  plot(o$penalties,o$AIC, main=main, type="b", xlab="penalty level", ...)
  abline(v=o$penalties[which.max(o$AIC)], lty=3, col="blue")
  
}

plot.ROC <- function(o, Theta.star, main="ROC Curve",type="b", ...) {

  ## Collect all the graphs in the same list
  if (is.list(o$networks[[1]])) {
    T <- length(o$networks[[1]])

    if (is.matrix(Theta.star)) {
      if (dim(Theta.star) == o$networks[[1]][[1]]) {
        Theta.star.l <- list()
        for (t in 1:T) {
          Theta.star.l[[t]] <- Theta.star
        }
        Theta.star <- Theta.star.l
      }
    }
    if (length(Theta.star) != T) {
      stop("\n Theta.star should be a list of",T,"networks")
    }
  } else {
    T == 1
  }
    
  fallout <- c()
  recall <- c()

  if (T == 1) {
    for (i in 1:length(o$networks)) {
      ana <- analysis(o$networks[[i]], Theta.star)
      recall <- c(recall,ana$recall)
      fallout <- c(fallout,ana$fallout)
    }
  } else {
    for (i in 1:length(o$networks)) {
      TP <- FP <- TN <- ne <- 0
      for (t in 1:T) {
        ana <- analysis(o$networks[[i]][[t]], Theta.star[[t]])
        TP <- TP + ana$TP
        FP <- FP + ana$FP
        TN <- TN + ana$TN
        ne <- ne + ana$n.edges
      }
      recall <- c(recall,TP/ne)
      fallout <- c(fallout,FP/(FP+TN))
    }    
  }
  
  plot(fallout,recall, main=main, type=type,xlim=0:1,ylim=0:1, ...)
  abline(v=fallout[which.max(o$BIC)], col="red",  lty=3)
  abline(v=fallout[which.max(o$AIC)], col="blue", lty=3)
  abline(h=recall[which.max(o$BIC)],  col="red",  lty=3)
  abline(h=recall[which.max(o$AIC)],  col="blue", lty=3)
  legend("bottomright",c("BIC","AIC"),lty=c(3,3),col=c("red","blue"))
  
}

plot.PR <- function(o,Theta.star,main="PR Curve", type="b", ...) {

  ## Collect all the graphs in the same list
  if (is.list(o$networks[[1]])) {
    T <- length(o$networks[[1]])

    if (is.matrix(Theta.star)) {
      if (dim(Theta.star) == o$networks[[1]][[1]]) {
        Theta.star.l <- list()
        for (t in 1:T) {
          Theta.star.l[[t]] <- Theta.star
        }
        Theta.star <- Theta.star.l
      }
    }    
    if (length(Theta.star) != T) {
      stop("\n Theta.star should be a list of",T,"networks")
    }
  } else {
    T == 1
  }
  
  precision <- c()
  recall <- c()

  if (T == 1) { 
    for (i in 1:length(o$networks)) {
      ana <- analysis(o$networks[[i]], Theta.star)
      recall <- c(recall,ana$recall)
      precision <- c(precision,ana$precision)
    }
  } else {
    
    for (i in 1:length(o$networks)) {
      TP <- FP <- ne <- 0
      for (t in 1:T) {
        ana <- analysis(o$networks[[i]][[t]], Theta.star[[t]])
        TP <- TP + ana$TP
        FP <- FP + ana$FP
        ne <- ne + ana$n.edges
      }
      pr  <- TP/(TP+FP)
      pr[is.nan(pr)] <- NA
      precision <- c(precision,pr)
      recall <- c(recall,TP/ne)
    }
  }
  
  plot(recall,precision, main=main, type=type,xlim=0:1,ylim=0:1, ...)
  abline(v=recall[which.max(o$BIC)]   , col="red",  lty=3)
  abline(v=recall[which.max(o$AIC)]   , col="blue", lty=3)
  abline(h=precision[which.max(o$BIC)], col="red",  lty=3)
  abline(h=precision[which.max(o$AIC)], col="blue", lty=3)
  legend("bottomleft",c("BIC","AIC"),lty=c(3,3),col=c("red","blue"))
}

analysis <-function(Theta, Theta.star,directed=!isSymmetric(Theta.star)) {
  
  if (!directed) {
    ## Removing diagonal terms...
    diag(Theta.star) <- 0
    diag(Theta)  <- 0
    
    ## List of infered edges
    Inferred.Edges <- which(abs(Theta[upper.tri(Theta)]) > 0)
    
    ## List of true edges
    True.Edges <- which(abs(Theta.star[upper.tri(Theta.star)]) > 0)
    n.Edges  <- length(True.Edges)
    n.Bar    <- length(Theta.star[upper.tri(Theta.star)]) - length(True.Edges)
    
  } else {
    ## List of infered edges
    Inferred.Edges <- which(abs(Theta) > 0)
    
    ## List of true edges
    True.Edges <- which(abs(Theta.star) > 0)
    n.Edges  <- length(True.Edges)
    n.Bar    <- length(Theta.star) - length(True.Edges)
  }  

  ## Computing TP, FP, ...
  TP <- sum(Inferred.Edges %in% True.Edges)
  FP <- length(Inferred.Edges) - TP
  TN <- n.Bar - FP
  FN <- n.Edges - TP

  recall  <- TP/n.Edges
  precision  <- TP/(TP+FP)
  precision[is.nan(precision)] <- NA
  fallout <- FP/(FP+TN)
  
  return(list(recall = recall, precision = precision, fallout = fallout,
              TP = TP, FP = FP, TN = TN, FN = FN, n.edges = n.Edges, n.zeros = n.Bar))
}

plot.sequence <- function(o, main="Sequencing display of the network") {

  ## The first network of the list is used as a prototype to create
  ## basical stuff and initialize
  proto    <- o$networks[[1]]
  directed <- o$control$edges.sym.rule == "NO"
  p        <-  nrow(proto)
  
  if (is.list(proto)) {
    stop("not available")
  }
  
  ## Get the clusters the reordenering indices
  class <- o$clusters
  index <- order(class)
  class <- class[index]
  class.col <- sample(light.palette(max(nlevels(class),1), black=FALSE))
  proto     <- proto[index,index]
  label     <- colnames(proto)
  
  ## Build the plot windows and display the vertices
  coord <- Gplot.default(proto, directed=directed,
                         class=class,class.col=class.col,
                         label=label, margin = 0.05, main = main )
  Env <- get("Gplot.graph", envir=Gplot.env)
  
  added.edges <- matrix(0, p, p)

  cat("\nPress \"Return\" or \"n+Return\" to go forward\nPress \"p+Return\" to go backward \nPress \"q+Return\" to quit\n\n")
  k <- 1
  forward <- TRUE
  while (k <= length(o$networks)) {
    ## get the current reordered network
    net <- o$networks[[k]][index,index]    

    if (forward) {
      new.edges <- net * ( ( (net != 0) &  (! added.edges) ) )
      added.edges <- (net != 0)
      if ( sum(new.edges != 0) > 0) {
        Gplot.add( new.edges, directed=directed,
                  class=class, class.col=class.col,
                  coord=Env$default.gr$network$coord,
                  margin = 0.05, main = main )
      }
    } else {
      #
      #  Backward
      #
      del.edges <- ( (added.edges) & (net == 0) )
      added.edges <- (net != 0)
      if (sum( del.edges !=0 ) > 0) {
        Gplot.delete(del.edges, del.col="white", directed=directed,
                     class=class, class.col=class.col,
                     coord=Env$default.gr$network$coord,
                     margin = 0.05, main = main )
      }
    }
    ctrl <-readline(paste("Penalty=",round(o$penalties[k],5),"",o$n.edges[k],"edges "))
    while(!(ctrl %in% c("q","n","p",""))) {
      ctrl <- readline(paste("quit/next/previous (q/n/p)? "))
    }
    if (ctrl == "q") {
      break;
    } 
    if (ctrl %in% c("n","")) {
      forward <- TRUE
      k <- k + 1
    } else 
    if (ctrl == "p") {
      forward <- FALSE
      k <- max(1,k-1)
    }
  }
  cat("\n")    
}

getNetwork <- function(object, selection=length(object$clusters), nodes=NULL) {
  
  if (!is.simone(object))
    stop("must be a simone object")
  
  directed <- object$control$edges.sym.rule == "NO"
  clusters <- object$clusters
  
  if (selection == "BIC") {
    select <- which.max(object$BIC)
    if (is.null(dim(object$n.edges)[2])) {
      name <- paste("BIC choice (", as.character(object$n.edges[select]), " edges)",
                    sep="")
    } else {
      name <- paste("BIC choice (")
    }    
  }

  if (selection == "AIC") {
    select <- which.max(object$AIC)
    if (is.null(dim(object$n.edges)[2])) {
      name <- paste("AIC choice (", as.character(object$n.edges[select]), " edges)",
                    sep="")
    } else {
      name <- paste("AIC choice (")
    }
  }

  if (is.numeric(selection)) {
    if (selection <= min(floor(apply(object$n.edges,1,mean))) ) {
      select <- 1
    } else {
      thres  <- min(selection,max(floor(apply(object$n.edges,1,mean))))
      select <- max(which(floor(apply(object$n.edges,1,mean))<=thres))
    }
    if (dim(object$n.edges)[2] >= 2) {
      cat("\nFound networks with",object$n.edges[select,],"edges.\n")
      name <- "Choice fixed to"
    } else {
      cat("\nFound a network with",object$n.edges[select],"edges.\n")
      name <- paste("Choice fixed to", as.character(object$n.edges[select]),"edges")
    }
  } else {
    if (!is.list(object$networks[[select]])) {
      cat("\nFound a network with",object$n.edges[select],"edges.\n")
      name <- paste(name,as.character(object$n.edges[select])," edges)", sep="")
    }
  }
  Theta <- object$networks[[select]] 

  ## select a sub network if required
  if (!is.null(nodes)) {
    if (is.list(Theta)) {
      for (t in 1:length(Theta)) {
        Theta[[t]] <- Theta[[t]][nodes,nodes]
      }
    } else {
      Theta <- Theta[nodes,nodes]
    }
    clusters <- clusters[nodes]
  }

  ## contruct the assosiated network object (a list if multitask
  if (is.list(Theta)) {
    net <- list()
    for (t in 1:length(Theta)) {
      if (is.numeric(selection)) {
        mname <- paste(name,as.character(object$n.edges[select,t]),
                       "edges for task", t)
      } else {
        mname <- paste(name, as.character(object$n.edges[select,t]),
                       " edges) for task",t, sep="")
      }
      A <- sign(abs(Theta[[t]]))
      if (!directed) {diag(A) <- 0}
      net[[t]] <- structure(list(A        = A,
                                 Theta    = Theta[[t]],
                                 directed = directed,
                                 clusters = clusters,
                                 name     = mname), class="simone.network")
    }
  }  else {
    A <- sign(abs(Theta))
    if (!directed) {diag(A) <- 0}
    net <-  structure(list(A        = A,
                           Theta    = Theta,
                           directed = directed,
                           clusters = clusters,
                           name     = name), class="simone.network")
  }
  
  return(net)
}
