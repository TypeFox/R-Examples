###{{{ plot.lvm

##' Plot path diagram
##'
##' Plot the path diagram of a SEM
##'
##'
##' @aliases plot.lvmfit
##' @param x Model object
##' @param diag Logical argument indicating whether to visualize variance
##' parameters (i.e. diagonal of variance matrix)
##' @param cor Logical argument indicating whether to visualize correlation
##' parameters
##' @param labels Logical argument indiciating whether to add labels to plot
##' (Unnamed parameters will be labeled p1,p2,...)
##' @param intercept Logical argument indiciating whether to add intercept
##' labels
##' @param addcolor Logical argument indiciating whether to add colors to plot
##' (overrides \code{nodecolor} calls)
##' @param plain if TRUE strip plot of colors and boxes
##' @param cex Fontsize of node labels
##' @param fontsize1 Fontsize of edge labels
##' @param noplot if TRUE then return \code{graphNEL} object only
##' @param graph Graph attributes (Rgraphviz)
##' @param attrs Attributes (Rgraphviz)
##' @param unexpr if TRUE remove expressions from labels
##' @param addstyle Logical argument indicating whether additional style should
##' automatically be added to the plot (e.g. dashed lines to double-headed
##' arrows)
##' @param Rgraphviz if FALSE igraph is used for graphics
##' @param init Reinitialize graph (for internal use)
##' @param layout Graph layout (see Rgraphviz or igraph manual)
##' @param edgecolor if TRUE plot style with colored edges
##' @param \dots Additional arguments to be passed to the low level functions
##' @author Klaus K. Holst
##' @keywords hplot regression
##' @examples
##'
##' if (interactive()) {
##' m <- lvm(c(y1,y2) ~ eta)
##' regression(m) <- eta ~ z+x2
##' regression(m) <- c(eta,z) ~ x1
##' latent(m) <- ~eta
##' labels(m) <- c(y1=expression(y[scriptscriptstyle(1)]),
##' y2=expression(y[scriptscriptstyle(2)]),
##' x1=expression(x[scriptscriptstyle(1)]),
##' x2=expression(x[scriptscriptstyle(2)]),
##' eta=expression(eta))
##' edgelabels(m, eta ~ z+x1+x2, cex=2, lwd=3,
##'            col=c("orange","lightblue","lightblue")) <- expression(rho,phi,psi)
##' nodecolor(m, vars(m), border="white", labcol="darkblue") <- NA
##' nodecolor(m, ~y1+y2+z, labcol=c("white","white","black")) <- NA
##' plot(m,cex=1.5)
##'
##' d <- sim(m,100)
##' e <- estimate(m,d)
##' plot(e)
##'
##' m <- lvm(c(y1,y2) ~ eta)
##' regression(m) <- eta ~ z+x2
##' regression(m) <- c(eta,z) ~ x1
##' latent(m) <- ~eta
##' plot(lava:::beautify(m,edgecol=FALSE))
##' }
##' @export
##' @method plot lvm
`plot.lvm` <-
  function(x,diag=FALSE,cor=TRUE,labels=FALSE,intercept=FALSE,addcolor=TRUE,plain=FALSE,cex,fontsize1=10,noplot=FALSE,graph=list(rankdir="BT"),
         attrs=list(graph=graph),
           unexpr=FALSE,
           addstyle=TRUE,Rgraphviz=lava.options()$Rgraphviz,init=TRUE,
           layout=c("dot","fdp","circo","twopi","neato","osage"),
           edgecolor=lava.options()$edgecolor,
           ...) {
    if (is.null(vars(x))) {
      message("Nothing to plot: model has no variables.")
      return(NULL)
    }
  index(x) <- reindex(x)
  if (length(index(x)$vars)<2) {
    message("Not available for models with fewer than two variables")
    return(NULL)
  }
  if (edgecolor) x <- beautify(x)
  myhooks <- gethook("plot.post.hooks")
  for (f in myhooks) {
    x <- do.call(f, list(x=x,...))
  }

  suppressWarnings(igraphit <- !Rgraphviz || !(requireNamespace("graph",quietly=TRUE)) || !(requireNamespace("Rgraphviz",quietly=TRUE)))
  if (igraphit) {
    if (!requireNamespace("igraph",quietly=TRUE)) {
      message("package 'Rgraphviz' or 'igraph' not available")
      return(NULL)
    }
    L <- igraph::layout.sugiyama(g <- igraph.lvm(x,...))$layout
    if (noplot) return(graph::updateGraph(g))
    dots <- list(...)
    if (is.character(layout))
      plot(g,layout=L,...)
    else plot(g,layout=layout,...)
    return(invisible(g))
  }
  if (init) {
    g <- finalize(x,diag=diag,cor=cor,addcolor=addcolor,intercept=intercept,plain=plain,cex=cex,fontsize1=fontsize1,unexpr=unexpr,addstyle=addstyle)
  } else {
    g <- Graph(x)
  }
  if  (labels) {
    AP <- matrices(x,paste0("p",seq_len(index(x)$npar)))
    mylab <- AP$P; mylab[AP$A!="0"] <- AP$A[AP$A!="0"]
    mylab[!is.na(x$par)] <- x$par[!is.na(x$par)]
    mylab[!is.na(x$covpar)] <- x$covpar[!is.na(x$covpar)]
    g <- edgelabels(g, lab=mylab)
  }
  if (lava.options()$debug) {
    plot(g)
  } else {
    ## graphRenderInfo(g)$recipEdges <- "distinct"
    .savedOpt <- options(warn=-1) ## Temporarily disable warnings as renderGraph comes with a stupid warning when labels are given as "expression"
    dots <- list(...)
    dots$attrs <- attrs
    dots$x <- g
    dots$recipEdges <- "distinct"
    if (attributes(g)$feedback) dots$recipEdges <- c("combine")
    if (is.null(dots$layoutType)) dots$layoutType <- layout[1]
    if (all(index(x)$A==0))
      dots$layoutType <- "circo"

    g <- do.call(getFromNamespace("layoutGraph","Rgraphviz"), dots)
    ## Temporary work around:
    graph::nodeRenderInfo(g)$fill <- graph::nodeRenderInfo(dots$x)$fill
    graph::nodeRenderInfo(g)$col <- graph::nodeRenderInfo(dots$x)$col
    graph::edgeRenderInfo(g)$col <- graph::edgeRenderInfo(dots$x)$col
    if (noplot)
      return(g)
    res <- tryCatch(Rgraphviz::renderGraph(g),error=function(e) NULL)
    options(.savedOpt)
  }
  ## if (!is.null(legend)) {
  ##   op <- par(xpd=TRUE)
  ##   legend(legend, c("Exogenous","Endogenous","Latent","Time to event"),
  ##          pt.cex=1.5, pch=15, lty=0, col=cols[1:4], cex=0.8)
  ##   par(op)
  ## }


  myhooks <- gethook("plot.hooks")
  for (f in myhooks) {
    do.call(f, list(x=x,...))
  }

  invisible(g)
}

###}}} plot.lvm

###{{{ plot.lvmfit

##' @export
`plot.lvmfit` <-
  function(x,diag=TRUE,cor=TRUE,type,noplot=FALSE,fontsize1=5,...) {
    .savedOpt <- options(warn=-1) ## Temporarily disable warnings as renderGraph comes with a stupid warning when labels are given as "expression"
    if (!requireNamespace("graph",quietly=TRUE)) {
      plot(Model(x),...)
      return(invisible(x))
    }
    g <- Graph(x)
    newgraph <- FALSE
    if (is.null(g)) {
      newgraph <- TRUE
      Graph(x) <- finalize(Model(x), diag=TRUE, cor=FALSE, fontsize1=fontsize1, ...)
    }
    if(noplot) return(Graph(x))
    if (newgraph) {
      if (missing(type))
        type <- "est"
      x <- edgelabels(x, type=type, diag=diag, cor=cor, fontsize1=fontsize1, ...)
    } else {
      if (!missing(type)) {
        x <- edgelabels(x, type=type, diag=diag, cor=cor, fontsize1=fontsize1, ...)
      }
    }
    g <- Graph(x)
    var <- rownames(covariance(Model(x))$rel)
    if (!cor) {
       delta <- 1
      for (r in seq_len(nrow(covariance(Model(x))$rel)-delta) ) {
        for (s in seq(r+delta,ncol(covariance(Model(x))$rel)) ) {
          if (covariance(Model(x))$rel[r,s]==1) {
            g <- graph::removeEdge(var[r],var[s], g)
            g <- graph::removeEdge(var[s],var[r], g)
          }
        }
      }
    }
    if (!diag) {
      for (r in seq_len(nrow(covariance(Model(x))$rel)) ) {
        if (graph::isAdjacent(g,var[r],var[r]))
          g <- graph::removeEdge(var[r],var[r],g)
      }
    }
    m <- Model(x); Graph(m) <- g
    g <- plot(m, diag=diag, cor=cor, fontsize1=fontsize1, init=FALSE, ...)
    options(.savedOpt)
    invisible(g)
  }

###}}} plot.lvmfit

###{{{ plot.multigroup

##' @export
plot.multigroup <- function(x,diag=TRUE,labels=TRUE,...) {
  k <- x$ngroup
  for (i in seq_len(k))
    plot(x$lvm[[i]],diag=diag,labels=labels, ...)
}

##' @export
plot.multigroupfit <- function(x,...) {
  plot(Model(x),...)
}

###}}}

###{{{ igraph.lvm

##' @export
igraph.lvm <- function(x,layout=igraph::layout.kamada.kawai,...) {
  requireNamespace("igraph",quietly=TRUE)
  oC <- covariance(x)$rel
  for (i in seq_len(nrow(oC)-1))
    for (j in seq(i+1,nrow(oC))) {
      if (oC[i,j]!=0) {
        x <- regression(x,vars(x)[i],vars(x)[j])
        x <- regression(x,vars(x)[j],vars(x)[i])
      }
    }
  g <- igraph::graph.adjacency(x$M,mode="directed")
  igraph::V(g)$color <- "lightblue"
  igraph::V(g)$label <- vars(x)
  igraph::V(g)$shape <- "rectangle"
  for (i in match(latent(x),igraph::V(g)$name)) {
      igraph::V(g)$shape[i] <- "circle"
      igraph::V(g)$color[i] <- "green"
  }
  endo <- index(x)$endogenous
  for (i in match(endo,igraph::V(g)$name)) {
      igraph::V(g)$color[i] <- "orange"
  }
  igraph::E(g)$label <- as.list(rep("",length(igraph::E(g))))
  oE <- edgelabels(x)
  for (i in seq_along(igraph::E(g))) {
    st <- as.character(oE[i])
    if (length(st)>0)
      igraph::E(g)$label[[i]] <- st
  }
  g$layout <- layout(g)
  return(g)
}

###}}} igraph.lvm


beautify <- function(x,col=c("lightblue","orange","yellowgreen"),border=rep("white",3),labcol=rep("darkblue",3),edgecol=TRUE,...) {
    nodecolor(x, exogenous(x), border=border[1], labcol=labcol[1]) <- NA
    nodecolor(x, endogenous(x), border=border[1], labcol=labcol[1]) <- NA
    nodecolor(x, latent(x), border=border[1], labcol=labcol[1]) <- NA
    trimmed <- gsub("[[:digit:]]*$","",vars(x))
    num <- c()
    for (i in seq_len(length(vars(x)))) {
        num <- c(num,gsub(trimmed[i],"",vars(x)[i]))
    }
    lab <- paste0(vars(x),"=",paste0("expression(",trimmed,"[scriptscriptstyle(",num,")])"),collapse=",")
    labels(x) <- eval(parse(text=paste("c(",lab,")")))
    if (!edgecol) return(x)
    iex <- index(x)$exo.idx
    ien <- index(x)$endo.idx
    ila <- index(x)$eta.idx
    for (i in iex) {
        for (j in which(x$M[i,]==1))
            edgelabels(x, to=vars(x)[j], from=rev(vars(x)[i]), cex=2, lwd=3,col=col[1]) <- NA
    }
    for (i in ien) {
        for (j in which(x$M[i,]==1))
            edgelabels(x, to=vars(x)[j], from=rev(vars(x)[i]), cex=2, lwd=3,col=col[2]) <- NA
    }
    for (i in ila) {
        for (j in which(x$M[i,]==1))
            edgelabels(x, to=vars(x)[j], from=rev(vars(x)[i]), cex=2, lwd=3,col=col[3]) <- NA
    }
    x
}
