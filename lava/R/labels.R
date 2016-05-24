###{{{ labels

##' Define labels of graph
##'
##' Alters labels of nodes and edges in the graph of a latent variable model
##'
##'
##' @aliases labels<- labels labels<-.default labels.lvm labels.lvmfit
##' labels.graphNEL edgelabels edgelabels<- edgelabels<-.lvm nodecolor
##' nodecolor<- nodecolor<-.default
##' @author Klaus K. Holst
##' @export
##' @keywords graphs aplot
##' @examples
##' m <- lvm(c(y,v)~x+z)
##' regression(m) <- c(v,x)~z
##' labels(m) <- c(y=expression(psi), z=expression(zeta))
##' nodecolor(m,~y+z+x,border=c("white","white","black"),
##'           labcol="white", lwd=c(1,1,5),
##'           lty=c(1,2)) <-  c("orange","indianred","lightgreen")
##' edgelabels(m,y~z+x, cex=c(2,1.5), col=c("orange","black"),labcol="darkblue",
##'            arrowhead=c("tee","dot"),
##'            lwd=c(3,1)) <- expression(phi,rho)
##' edgelabels(m,c(v,x)~z, labcol="red", cex=0.8,arrowhead="none") <- 2
##' if (interactive()) {
##'     plot(m,addstyle=FALSE)
##' }
##' 
##' m <- lvm(y~x)
##' labels(m) <- list(x="multiple\nlines")
##' if (interactive()) {
##' op <- par(mfrow=c(1,2))
##' plot(m,plain=TRUE)
##' plot(m)
##' par(op)
##' 
##' d <- sim(m,100)
##' e <- estimate(m,d)
##' plot(e,type="sd")
##' }
##' @param object \code{lvm}-object.
##' @param value node label/edge label/color
##' @param to Formula specifying outcomes and predictors defining relevant
##' edges.
##' @param \dots Additional arguments (\code{lwd}, \code{cex}, \code{col},
##' \code{labcol}), \code{border}.
##' @param var Formula or character vector specifying the nodes/variables to
##' alter.
##' @param border Colors of borders
##' @param labcol Text label colors
##' @param shape Shape of node
##' @param lwd Line width of border
##' @usage
##' \method{labels}{default}(object, ...) <- value
##' \method{edgelabels}{lvm}(object, to, ...) <- value
##' \method{nodecolor}{default}(object, var=vars(object),
##' border, labcol, shape, lwd, ...) <- value
`labels<-` <- function(object,...,value) UseMethod("labels<-")

##' @export
`labels<-.default` <- function(object,...,value) {
  labels(object,value)
}

##' @export
labels.graphNEL <- function(object,lab=NULL,...) {
  if (is.null(lab))
    return(graph::nodeRenderInfo(object)$label)
  graph::nodeRenderInfo(object) <- list(label=lab)
  names(graph::nodeRenderInfo(object)$label) <- graph::nodes(object);
  return(object)
}

##' @export
labels.lvmfit <- function(object,lab=NULL,...) {
  if (is.null(lab)) return(object$noderender$label)
  object$noderender$label <- lab
  return(object)
}

##' @export
`labels.lvm` <- function(object,lab=NULL,...) {
  if (is.null(lab))
    return(object$noderender$label)
  if (is.null(object$noderender$label))
    object$noderender$label <- lab
  else
    object$noderender$label[names(lab)] <- lab
  return(object)
}
###}}} labels

###{{{ edgelabels

##' @export
"edgelabels<-.lvmfit" <- function(object,to,from,est=TRUE,edges=NULL,cex=1,...,value) {
  if (is.null(edges))  {
    if (inherits(to,"formula")) {
      yy <- decomp.specials(getoutcome(to))
      from <- setdiff(all.vars(to),yy)
      to <- yy
    }
    edges <- paste(from,to,sep="~")
  }

  edges. <- paste0("\"", edges, "\"")
  fromto <- edge2pair(edges)
  val <- c()
  for (i in seq_along(edges)) {
    val <- c(val,
             formatC(effects(object,from=fromto[[i]][1],to=fromto[[i]][2],silent=TRUE)$directef[[1]])
             )
  }
  if (est)
    mytext <- paste("c(", paste(paste0(edges.,"=expression(",as.character(value),"==\"",val,"\")"),collapse=","),")")
  else
    mytext <- paste("c(", paste(paste0(edges.,"=expression(",as.character(value),")"),collapse=","),")")
  graph::edgeRenderInfo(Graph(object))$label <- eval(parse(text=mytext))
  graph::edgeRenderInfo(Graph(object))$cex[edges] <- cex
  return(object)
}

##' @export
edgelabels.lvmfit <- function(object,value,type,pthres,intercept=FALSE,format.fun=formatC,...) {
    if (!missing(value)) {
        edgelabels(object,...) <- value
        return(object)
    }
    if (missing(type))
        return(graph::edgeRenderInfo(Graph(object))$label)
    

    Afix <- index(object)$A ## Matrix with fixed parameters and ones where parameters are free
    Pfix <- index(object)$P ## Matrix with fixed covariance parameters and ones where param
    mfix <- index(object)$v0
    
    npar.mean <- index(object)$npar.mean
    Par <- object$coef
    mpar <- c()
    if (npar.mean>0) {
        mpar <- do.call(format.fun,list(Par[seq_len(npar.mean)]))
        Par <- Par[-seq_len(npar.mean),,drop=FALSE]
    }
    Par <-
        switch(type,
               sd = paste0(do.call(format.fun,list(Par[,1,drop=FALSE])), " (", do.call(format.fun,list(Par[,2,drop=FALSE])), ")"),
               est = do.call(format.fun,list(Par[,1,drop=FALSE])),
               pval = do.call(format.fun,list(Par[,4,drop=FALSE])),
               name = rownames(Par),
               none = ""
               )
    AP <- matrices(Model(object), Par,mpar) ## Ignore expar
    A <- AP$A; P <- AP$P
    P[exogenous(object),exogenous(object)] <- NA
    
    gr <- finalize(Model(object), ...)
    Anz <- A; Anz[Afix==0] <- NA    
    gr <- edgelabels(gr, lab=Anz)
    Pnz <- P; Pnz[Model(object)$cov==0] <- NA
    if (intercept) {
        idx <- which(!is.na(diag(Pnz)))
        diag(Pnz)[idx] <- paste(paste0("[",AP$v[idx],"]"),diag(Pnz)[idx],sep="\n")
    }
    gr <- edgelabels(gr, lab=Pnz, expr=!intercept)
    Graph(object) <- gr
    return(object)
}

##' @export
`edgelabels` <- function(object, ...) UseMethod("edgelabels")

##' @export
`edgelabels<-` <- function(object,...,value) UseMethod("edgelabels<-")

##' @export
`edgelabels<-.lvm` <- function(object,to,...,value) {
  edgelabels(object,to=to, lab=value,...)
}

##' @export
`edgelabels<-.graphNEL` <- function(object,...,value) {
  edgelabels(object,lab=value,...)
}

##' @export
`edgelabels.graphNEL` <- function(object, lab=NULL, to=NULL, from=NULL, cex=1.5, lwd=1, lty=1, col="black", labcol="black", arrowhead="closed",
                                  expr=TRUE,
                                  debug=FALSE,...) {
  if (is.null(lab)) {
    return(graph::edgeRenderInfo(object)$label)
  }
  if (inherits(to,"formula")) {
    yy <- decomp.specials(getoutcome(to))
    from <- all.vars(to[[3]])##setdiff(all.vars(to),yy)
    if (length(from)==0) from <- yy
    to <- yy
  }

  M <- as(object, Class="matrix")
  nodes <- graph::nodes(object)

  if (is.null(graph::edgeRenderInfo(object)$label))
      graph::edgeRenderInfo(object)$label <- expression()


  if (!is.null(lab)) {
    if (!is.null(from) & !is.null(to)) {
      estr <- paste0("\"",from,"~",to,"\"")
      estr2 <- paste0(from,"~",to)
      if (length(lab)!=length(estr2)) lab <- rep(lab,length(estr2))
      if (length(col)!=length(estr2)) col <- rep(col,length(estr2))
      if (length(cex)!=length(estr2)) cex <- rep(cex,length(estr2))
      if (length(lwd)!=length(estr2)) lwd <- rep(lwd,length(estr2))
      if (length(lty)!=length(estr2)) lty <- rep(lty,length(estr2))
      if (length(arrowhead)!=length(estr2))
          arrowhead <- rep(arrowhead,length(estr2))
      if (length(labcol)!=length(estr2))
          labcol <- rep(labcol,length(estr2))

      curedges <- names(graph::edgeRenderInfo(object)$label)
       Debug(estr,debug)

      estr2.idx <- which(estr2%in%curedges)
      newstr.idx <- setdiff(seq_along(estr2),estr2.idx)
      newstr <- estr2[newstr.idx]
      estr2 <- estr2[estr2.idx]
      if (length(estr2)>0) {
        if (!is.null(lab))
          graph::edgeRenderInfo(object)$label[estr2] <- lab[estr2.idx]
        if (!is.null(cex))
            graph::edgeRenderInfo(object)$cex[estr2] <- cex[estr2.idx]
        if (!is.null(col))
            graph::edgeRenderInfo(object)$col[estr2] <- col[estr2.idx]
        if (!is.null(lwd))
            graph::edgeRenderInfo(object)$lwd[estr2] <- lwd[estr2.idx]
        if (!is.null(lty))
            graph::edgeRenderInfo(object)$lty[estr2] <- lty[estr2.idx]
        if (!is.null(labcol))
            graph::edgeRenderInfo(object)$textCol[estr2] <- labcol[estr2.idx]
        if (!is.null(arrowhead))
            graph::edgeRenderInfo(object)$arrowhead[estr2] <- arrowhead[estr2.idx]
      }
      if (length(newstr)>0) {

          if (!is.null(lab))
              graph::edgeDataDefaults(object)$futureinfo$label[newstr] <-
                  lab[newstr.idx]
        if (!is.null(cex))
            graph::edgeDataDefaults(object)$futureinfo$cex[newstr] <-
                cex[newstr.idx]
        if (!is.null(col))
            graph::edgeDataDefaults(object)$futureinfo$col[newstr] <-
                col[newstr.idx]
        if (!is.null(lwd))
            graph::edgeDataDefaults(object)$futureinfo$lwd[newstr] <-
                lwd[newstr.idx]
        if (!is.null(lty))
            graph::edgeDataDefaults(object)$futureinfo$lty[newstr] <-
                lty[newstr.idx]
        if (!is.null(labcol))
            graph::edgeDataDefaults(object)$futureinfo$textCol[newstr] <-
                labcol[newstr.idx]
        if (!is.null(arrowhead))
            graph::edgeDataDefaults(object)$futureinfo$arrowhead[newstr] <-
                arrowhead[newstr.idx]
      }
      return(object)
    }

    ## Used by "edgelabels.lvmfit"
    for (r in seq_len(nrow(M)))
      for (s in seq_len(ncol(M))) {
        if (M[r,s]!=0 & !is.na(lab[r,s])) {
          estr <- paste0("\"",nodes[r],"~",nodes[s],"\"")
          estr2 <- paste0(nodes[r],"~",nodes[s])
          Debug(estr, debug)
          if (expr)
            st <- eval(parse(text=paste0("expression(",lab[r,s],")")))
          else
            st <- lab[r,s]
          graph::edgeRenderInfo(object)$label[estr2] <- st
        }
      }
  }

  return(object)
}



##' @export
`edgelabels.lvm` <- function(object, lab=NULL, to=NULL, from=NULL,
                             cex=1.5, lwd=1, lty=1, col="black",
                             labcol="black", arrowhead="closed",
                             expr=TRUE, debug=FALSE,...) {
  if (is.null(lab)) {
    return(object$edgerender$label)
  }
  if (inherits(to,"formula")) {
    yy <- decomp.specials(getoutcome(to))
    from <- all.vars(to[[3]])##setdiff(all.vars(to),yy)
    if (length(from)==0) from <- yy
    to <- yy
  }

  M <- object$M
  nodes <- colnames(M)

  if (is.null(object$edgerender$label))
    object$edgerender$label <- expression()


  if (!is.null(lab)) {
    if (!is.null(from) & !is.null(to)) {
      estr <- paste0("\"",from,"~",to,"\"")
      estr2 <- paste0(from,"~",to)
      if (length(lab)!=length(estr2)) lab <- rep(lab,length(estr2))
      if (length(col)!=length(estr2)) col <- rep(col,length(estr2))
      if (length(cex)!=length(estr2)) cex <- rep(cex,length(estr2))
      if (length(lwd)!=length(estr2)) lwd <- rep(lwd,length(estr2))
      if (length(lty)!=length(estr2)) lty <- rep(lty,length(estr2))
      if (length(labcol)!=length(estr2)) labcol <- rep(labcol,length(estr2))
      if (length(arrowhead)!=length(estr2)) arrowhead <- rep(arrowhead,length(estr2))

      curedges <- names(object$edgerender$label)
       Debug(estr,debug)

      estr2.idx <- which(estr2%in%curedges)
      newstr.idx <- setdiff(seq_along(estr2),estr2.idx)
      newstr <- estr2[newstr.idx]
      estr2 <- estr2[estr2.idx]
      if (length(estr2)>0) {
        if (!is.null(lab))
          object$edgerenderlabel[estr2] <- lab[estr2.idx]
        if (!is.null(cex))
          object$edgerender$cex[estr2] <- cex[estr2.idx]
        if (!is.null(col))
          object$edgerender$col[estr2] <- col[estr2.idx]
        if (!is.null(lwd))
          object$edgerender$lwd[estr2] <- lwd[estr2.idx]
        if (!is.null(lty))
          object$edgerender$lty[estr2] <- lty[estr2.idx]
        if (!is.null(labcol))
          object$edgerender$textCol[estr2] <- labcol[estr2.idx]
        if (!is.null(arrowhead))
          object$edgerender$arrowhead[estr2] <- arrowhead[estr2.idx]
      }
      if (length(newstr)>0) {

        if (!is.null(lab))
          object$edgerender$futureinfo$label[newstr] <-
            lab[newstr.idx]
        if (!is.null(cex))
          object$edgerender$futureinfo$cex[newstr] <-
            cex[newstr.idx]
        if (!is.null(col))
          object$edgerender$futureinfo$col[newstr] <-
            col[newstr.idx]
        if (!is.null(lwd))
          object$edgerender$futureinfo$lwd[newstr] <-
            lwd[newstr.idx]
        if (!is.null(lty))
          object$edgerender$futureinfo$lty[newstr] <-
            lty[newstr.idx]
        if (!is.null(labcol))
          object$edgerender$futureinfo$textCol[newstr] <-
            labcol[newstr.idx]
        if (!is.null(arrowhead))
          object$edgerender$futureinfo$arrowhead[newstr] <-
            arrowhead[newstr.idx]
      }
      return(object)
    }

    ## Used by "edgelabels.lvmfit"
    for (r in seq_len(nrow(M)))
      for (s in seq_len(ncol(M))) {
        if (M[r,s]!=0 & !is.na(lab[r,s])) {
          estr <- paste0("\"",nodes[r],"~",nodes[s],"\"")
          estr2 <- paste0(nodes[r],"~",nodes[s])
          Debug(estr, debug)
          if (expr)
            st <- eval(parse(text=paste0("expression(",lab[r,s],")")))
          else
            st <- lab[r,s]
          object$edgerender$label[estr2] <- st
        }
      }
  }
  return(object)
}

###}}} edgelabels
