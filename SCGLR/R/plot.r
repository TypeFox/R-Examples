if(getRversion()>="2.15.1") {
  # remove warnings due to ggplot2 syntax strangeness
  utils::globalVariables(c("comp","y","label","angle","hjust"))
}


#' SCGLR generic plot
#' @export
#' @importFrom stats cor aggregate
#' @importFrom grid circleGrob gpar
#' @method plot SCGLR
#' @description SCGLR generic plot
#' @param x an object from SCGLR class.
#' @param \dots optional arguments (see \link{customize}).
#' @param style named list of values used to customize the plot (see \link{customize})
#' @param plane a size-2 vector (or comma separated string) indicating which components are plotted (eg: c(1,2) or "1,2").
#' @return an object of class \code{\link{ggplot}}.
#' @examples \dontrun{
#' library(SCGLR)
#' 
#' # load sample data
#' data(genus)
#' 
#' # get variable names from dataset
#' n <- names(genus)
#' ny <- n[grep("^gen",n)]    # Y <- names that begins with "gen"
#' nx <- n[-grep("^gen",n)]   # X <- remaining names
#' 
#' # remove "geology" and "surface" from nx
#' # as surface is offset and we want to use geology as additional covariate
#' nx <-nx[!nx%in%c("geology","surface")]
#' 
#' # build multivariate formula
#' # we also add "lat*lon" as computed covariate
#' form <- multivariateFormula(ny,c(nx,"I(lat*lon)"),c("geology"))
#' 
#' # define family
#' fam <- rep("poisson",length(ny))
#' 
#' genus.scglr <- scglr(formula=form,data = genus,family=fam, K=4,
#'  offset=genus$surface)
#' 
#' summary(genus.scglr)
#' 
#' barplot(genus.scglr)
#' 
#' plot(genus.scglr)
#' 
#' plot(genus.scglr, predictors=TRUE, factor=TRUE)
#' 
#' pairs(genus.scglr)
#'
#' } 
plot.SCGLR <- function(x, ..., style=getOption("plot.SCGLR"), plane=c(1,2)) {
  data <- x
  if(class(data)!="SCGLR") 
    stop("This plot function need an SCGLR result")
  
  if(dim(data$compr)[2]<2)
    stop("At least two axes are needed for this kind of plot!")
  
  # customize from remaining arguments
  dots <- list(...)
  # merged with style argument
  if(!is.null(style)) {
    style[names(dots)] <- dots
    dots <- style
  }
  allowed_cust <- c(
    "title","expand","labels.offset","labels.auto","labels.size","threshold",
    "observations",
    "observations.factor","observations.size","observations.color",
    "observations.alpha",
    "predictors",
    "predictors.color","predictors.arrows","predictors.arrows.color",
    "predictors.labels","predictors.labels.color","predictors.labels.size",
    "predictors.labels.auto",
    "predictors.alpha","predictors.arrows.alpha","predictors.labels.alpha",
    "covariates",
    "covariates.color","covariates.arrows","covariates.arrows.color",
    "covariates.labels","covariates.labels.color","covariates.labels.size",
    "covariates.labels.auto",
    "covariates.alpha","covariates.arrows.alpha","covariates.labels.alpha",
    "factor",
    "factor.points","factor.points.size","factor.points.shape",
    "factor.labels","factor.labels.color","factor.labels.size"
  )
  match_cust <- function(key) {
    # build regexp  foo.bar  --> ^foo[^.]*\\.bar[^.]*$
    key2 <- paste("^",gsub("\\.","[^.]*\\\\.",key),"[^.]*$",sep="")
    # find matching cust key
    key2 <- grep(key2,allowed_cust,value=TRUE)
    # if not found (length==0) or ambiguity (length >1) revert to original
    if(length(key2)!=1)
      key2 <- key
    key2
  }
  # resolve abbreviated customization key names
  names(dots) <- sapply(names(dots),match_cust)
  # get custom value from key (if not found return default value)
  cust <- function(key, def=NULL) {
    if(is.null(key)||is.null(dots[[key]])) {
      return(def) 
    } else {
      return(dots[[key]])
    }
  }
  # check if key has a value or is not false
  has_cust <- function(key, def=NULL) {
    key <- cust(key, def)
    !is.null(key) && (!is.logical(key) || key)
  }
  
  labels.offset <- cust("labels.offset",0.01)
  labels.auto <- cust("labels.auto",TRUE)
  labels.size <- cust("labels.size",1)
    
  # process plane
  if(is.character(plane)) {
    plane <- as.integer(trim(unlist(strsplit(plane,","))))
  }
  
  # sanity checking
  if(length(plane) !=2 ) {
    stop("Plane should have two components!")
  }
  if((min(plane)<1) || (max(plane)>ncol(data$compr))) {
    stop("Invalid components for plane!")
  }
  
  # plan
  axis_names <- colnames(data$compr)[plane]
  
  # check factor
  factor <- cust("factor",NULL)
  if(!is.null(factor)&&(is.character(factor)||factor)) {
    if(is.logical(factor)) {
      if(is.null(data$xFactors))
        stop("No factor in data!")
      factor <- names(data$xFactors)[1]
      warning("No factor given, assuming first one! (",factor,")!")
    } else {
      if(!factor %in% names(data$xFactors))
        stop("Invalid factor!")
    }
  }
  
  # inertia
  inertia <- data$inertia[plane]

  # layers
  covariates <- has_cust("covariates",TRUE)
  predictors <- has_cust("predictors",FALSE)
  observations <- has_cust("observations",FALSE)
  
  # build base plot
  p <- qplot((-1:1)*cust("expand",1.0), (-1:1)*cust("expand",1.0), geom="blank")+
    coord_fixed()+
    xlab(paste(axis_names[1],"(",round(100*inertia[1],2),"%)")) + 
    ylab(paste(axis_names[2],"(",round(100*inertia[2],2),"%)")) 

  if(observations) {
    p <- p +
      geom_hline(yintercept=0,size=0.5)+
      geom_vline(xintercept=0,size=0.5)
  }
  if(covariates||predictors) {
    p <- p +
    # thicker x unit arrow
      geom_segment(aes(x=-1.1,xend=1.1,y=0,yend=0),size=1,arrow=arrow(length=unit(0.02,"npc")))+
    # thicker y unit arrow
      geom_segment(aes(y=-1.1,yend=1.1,x=0,xend=0),size=1,arrow=arrow(length=unit(0.02,"npc")))
  }
  
  # plot title
  if(has_cust("title", FALSE)) {
    p <- p + ggtitle(cust("title"))
  } else {
    if(observations) {
      if(covariates||predictors) {
        p <- p + ggtitle("Mixed individual and \ncorrelation plot\n")
      } else {
        p <- p + ggtitle("Individual plot\n")
      }
    } else {
      if(covariates||predictors) {
        p <- p + ggtitle("Correlation plot\n")      
      } else {
        p <- p + ggtitle("SCGLR Plot\n")      
      }
    }
  }

  # add unit circle
  if(covariates||predictors)
    p <- p + annotation_custom(circleGrob(r=0.5,gp=gpar(fill=NA)),-1,1,-1,1)

  # add threshold circle
  if((covariates||predictors)&&has_cust("threshold"))
    p <- p + annotation_custom(circleGrob(r=0.5*cust("threshold"),gp=gpar(fill=NA,lty="dashed")),-1,1,-1,1)
  
  # add observations
  if(observations) {
    obs <- as.data.frame(data$compr[,plane])
    names(obs) <- c("x","y")
    
    # colored according to factor ?
    if(!is.null(factor)&&(cust("observations.factor",FALSE))) {
      obs <- cbind(obs,data$xFactors[factor])
      p <- p + geom_point(
        aes_string(x="x",y="y",color=factor),
        data=obs,
        size=cust("observations.size",1),
        alpha=cust("observations.alpha",1)
      )
    } else {
      p <- p + geom_point(
        aes(x=x,y=y),
        data=obs,
        size=cust("observations.size",1),
        color=cust("observations.color","black"),
        alpha=cust("observations.alpha",1)
      )
    }
    tmp <- data.frame(x=c(0,0),y=range(obs$y))
    p <- p + geom_line(aes(x,y),data=tmp)
  }
  
  # add linear predictor arrows
  if(predictors) {
    predictors <- cust("predictors")
    co <- as.data.frame(cor(data$lin.pred, data$compr[,plane]))
    names(co) <- c("x", "y")
    co$norm <- sqrt(co$x^2+co$y^2)
    co$label <- rownames(co)
    co$color <- cust("predictors.color","red")
    co$alpha <- cust("predictors.alpha",1)
    if(cust("predictors.arrows",TRUE)) {
      co$arrows.color <- cust("predictors.arrows.color",co$color)
      co$arrows.alpha <- cust("predictors.arrows.alpha",co$alpha)
    }
    if(cust("predictors.labels",TRUE)) {
      co$labels.color <- cust("predictors.labels.color",co$color)
      co$labels.alpha <- cust("predictors.labels.alpha",co$alpha)
      co$labels.size <- cust("predictors.labels.size",4*labels.size)
    }
    
    # adjust label position
    co$angle <- atan2(co$y,co$x)*180/pi
    co$hjust <- ifelse(abs(co$angle)>90,1,0)
    if(cust("predictors.labels.auto",labels.auto)) {
      co$angle <- ifelse(abs(co$angle)>90,co$angle+180,co$angle)
    } else {
      co$angle <- 0
    }
    
    # filter according to user's will
    if(is.character(predictors)) {
      co <- co[co$label %in% predictors,]      
    }
    
    # filter according to given threshold
    if(cust("threshold",FALSE)>0)
      co <- co[co$norm>cust("threshold"),]
    
    if(nrow(co)==0) {
      warning("No predictors with correlation higher than threshold value ",cust("threshold"),"!")
    } else {      
      # draw arrows ?
      if(cust("predictors.arrows",TRUE)) {
        p <- p + geom_segment(
          aes(x=0,y=0,xend=x,yend=y),
          data=co,
          color=co$arrows.color,
          alpha=co$arrows.alpha,
          arrow=arrow(length=unit(0.02,"npc"))
        )
      } else {
        # reset label position
        co$hjust <- 0.5
        co$angle <- 0
      }
      
      # draw labels ?
      if(cust("predictors.labels",TRUE)) {
        p <- p + geom_text(
          aes(x=x*(1+labels.offset/norm),y=y*(1+labels.offset/norm),label=label,angle=angle,hjust=hjust),
          data=co,
          color=co$labels.color,
          alpha=co$labels.alpha,
          size=co$labels.size
        )
      }
    }    
  }
    
  # add co-variate arrows
  if(covariates) {
    covariates <- cust("covariates")
    co <- as.data.frame(cor(data$xNumeric, data$compr[,plane]))
    names(co) <- c("x", "y")
    co$norm <- sqrt(co$x^2+co$y^2)
    co$label <- rownames(co)
    co$color <- cust("covariates.color","black")
    co$alpha <- cust("covariates.alpha",1)
    if(cust("covariates.arrows",TRUE)) {
      co$arrows.color <- cust("covariates.arrows.color",co$color)
      co$arrows.alpha <- cust("covariates.arrows.alpha",co$alpha)
    }
    if(cust("covariates.labels",TRUE)) {
      co$labels.color <- cust("covariates.labels.color",co$color)
      co$labels.alpha <- cust("covariates.labels.alpha",co$alpha)
      co$labels.size <- cust("covariates.labels.size",4*labels.size)
    }

    # adjust label position
    co$angle <- atan2(co$y,co$x)*180/pi
    co$hjust <- ifelse(abs(co$angle)>90,1,0)
    if(cust("covariates.labels.auto",labels.auto)) {
      co$angle <- ifelse(abs(co$angle)>90,co$angle+180,co$angle)
    } else {
      co$angle <- 0
    }
    
    # filter according to users's will
    if(is.character(covariates)) {
      co <- co[co$label %in% covariates,]
    }

    # filter according to given threshold
    if(cust("threshold",FALSE)>0)
      co <- co[co$norm>cust("threshold"),]
    
    if(nrow(co)==0) {
      warning("No correlation higher than threshold value ",cust("threshold"),"!")
    } else {
      # draw arrows ?
      if(cust("covariates.arrows",TRUE)) {
        p <- p + geom_segment(
          aes(x=0,y=0,xend=x,yend=y),
          data=co,
          color=co$arrows.color,
          alpha=co$arrows.alpha,
          arrow=arrow(length=unit(0.02,"npc"))
        )
      } else {
        # reset label adjust
        co$hjust <- 0.5
        co$angle <- 0
      }
      
      # draw labels ?
      if(cust("covariates.labels",TRUE)) {
        p <- p + geom_text(
          aes(x=x*(1+labels.offset/norm),y=y*(1+labels.offset/norm),label=label,angle=angle,hjust=hjust),
          data=co,
          color=co$labels.color,
          size=co$labels.size,
          alpha=co$labels.alpha
        )
      }
    }
  }

  # add factors
  if(!is.null(factor)) {
    bary <- aggregate(data$compr[,plane],data$xFactors[factor],mean)
    names(bary) <- c(factor,"x","y")
    
    # draw points ?
    if(cust("factor.points",TRUE)) {
      p <- p + geom_point(
        aes_string(x="x",y="y",color=factor),
        data=bary,
        size=cust("factor.points.size",4),
        shape=cust("factor.points.shape",13)
        )
    }
    
    # draw labels ?
    if(cust("factor.labels",TRUE)) {
      p <- p + geom_label(
        aes_string(x="x",y="y",label=factor),
        data=bary,
        color=cust("factor.labels.color","black"),
        size=cust("factor.labels.size",4*labels.size)
        )
    }
  }
  
  # return plot
  p
}

#' @title Barplot of percent of overall X variance captured by component
#' @export
#' @method barplot SCGLR
#' @description A custom plot for SCGLR objetcs
#' @param height object of class 'SCGLR', usually a result of running \code{\link{scglr}}.
#' @param \dots optional arguments.
#' @param plane a size-2 vector (or comma separated string) indicating which components are plotted (eg: c(1,2) or "1,2").
#' @return an object of class ggplot.
#' @seealso For barplot application see examples in \code{\link{plot.SCGLR}}.
#' @importFrom scales percent_format
barplot.SCGLR <- function(height, ..., plane=NULL) {
  if(!is.null(plane)) {
    # process plane
    if(is.character(plane)) {
      plane <- as.integer(trim(unlist(strsplit(plane, ","))))
    }
    
    # sanity checking
    if(length(plane) !=2 ) {
      stop("Plane should have two components!")
    }
    if((min(plane)<1) || (max(plane)>length(height$inertia))) {
      stop("Invalid components for plane!")
    }
  }
  
  # build data frame from inertia
  inertia <- data.frame(
      inertia=height$inertia,
      comp=1:length(height$inertia),
      plane_color="black",
      stringsAsFactors = F)
  
  # recolor unused component in gray
  if(!is.null(plane))
    inertia$plane[-plane] <- "gray"
  
  ggplot(data=inertia)+geom_bar(aes(comp,inertia),fill=inertia$plane,stat="identity",width=0.5) +
    scale_x_discrete(labels=names(height$inertia),limits=1:length(inertia$comp))+
    scale_y_continuous(labels=percent_format())+
    labs(x="Components", y="Inertia", title="Inertia per component\n", ...)#+
    #geom_text(aes(comp,inertia,label=round(100*inertia,2)),vjust=0)
}

#' @title Pairwise scglr plot on components
#' @export
#' @importFrom utils combn
#' @method pairs SCGLR
#' @description Pairwise scglr plot on components
#' @param x object of class 'SCGLR', usually a result of running \code{\link{scglr}}.
#' @param \dots optionally, further arguments forwarded to \code{link{plot.SCGLR}}.
#' @param nrow number of rows of the grid layout.
#' @param ncol number of columns of the grid layout.
#' @param components vector of integers selecting components to plot (default is all components).
# @return an object of class ggplot.
#' @seealso For pairs application see examples in \code{\link{plot.SCGLR}} 
pairs.SCGLR <- function(x, ..., nrow=NULL, ncol=NULL, components=NULL) {
  prm <- list(...)

  nr <- nrow
  nc <- ncol

  prm["x"] <- list(x)
  prm[["components"]] <- NULL
  
  # pairs of components
  if(is.null(components)) {
    ncomp <- ncol(x$compr)
  } else {
    ncomp <- components
  }
  # sanity check
  if((min(ncomp)<1) || (max(ncomp)>ncol(x$compr))) {
    stop("Invalid components for plane!")
  }
  # build pairs
  cmp_pairs <- combn(ncomp, 2, simplify=FALSE)
  
  # build plot list
  one_plot <- function(cmp_pair) {
    do.call("plot.SCGLR", c(prm, plane=list(cmp_pair),title=paste(cmp_pair,collapse = "/")))
  }
  plots <- lapply(cmp_pairs, one_plot)
  
  # arrange them in a grid
  if(is.null(nc) & is.null(nr)) { 
    nr <- as.integer(sqrt(length(plots)))
  }
#  if(require("gridExtra",quietly = TRUE)) {
#    do.call("arrangeGrob", c(plots, nrow=nr, ncol=nc,main="toto"))
#  } else {
    do.call("arrange", c(plots, nrow=nr, ncol=nc))
#  }
}

## equivalent du par pour les ggplot2
#' @importFrom grid viewport grid.newpage pushViewport viewport grid.layout
vp.layout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)
arrange <- function(..., nrow=NULL, ncol=NULL, as.table=FALSE) {
  dots <- list(...)
  n <- length(dots)
  if(is.null(nrow) & is.null(ncol)) { nrow = floor(n/2) ; ncol = ceiling(n/nrow)}
  if(is.null(nrow)) { nrow = ceiling(n/ncol)}
  if(is.null(ncol)) { ncol = ceiling(n/nrow)}
  ## NOTE see n2mfrow in grDevices for possible alternative
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow,ncol) ) )
  ii.p <- 1
  for(ii.row in seq(1, nrow)){
    ii.table.row <- ii.row
    if(as.table) {ii.table.row <- nrow - ii.table.row + 1}
    for(ii.col in seq(1, ncol)){
      ii.table <- ii.p
      if(ii.p > n) break
      print(dots[[ii.table]], vp=vp.layout(ii.table.row, ii.col))
      ii.p <- ii.p + 1
    }
  }
  NULL
}
