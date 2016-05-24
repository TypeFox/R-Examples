plot.stdf <- function(x,
                      sel          = NULL,
                      plot.alpha   = TRUE,
                      mymar        = c(5,8,4,2),
                      pch.p        = c(1,5),
                      col.p        = c("black", "red"),
                      colorspace   = NULL,
                      cex.p        = 1.25,
                      line.stdf    = 1,
                      line.alpha   = 4,
                      with.legend  = TRUE,
                      legend.label = c("before", "after"),
                      legend.cex   = 1,
                      legend.xy    = NULL,
                      ...)
{

  object <- x

  ## ################
  ## check colorspace
  if (!is.null(colorspace)){
    suppressWarnings(require( "colorspace", character.only=TRUE ))
    message("Argument 'col.p' is ignored unless it has numeric values.")
    
    if (colorspace==TRUE){  
      if (any(apply(as.data.frame(col.p),1,is.numeric))==FALSE){
        col.p <- sample(rainbow_hcl(20),2)
        message("Plot colors are randomly chosen since 'col.p' are strings.")
      }else{
        col.p <- rainbow_hcl(20)[col.p]
      }
    }else{
      
      if (colorspace==FALSE){
        if (any(apply(as.data.frame(col.p),1,is.numeric))==FALSE){
          col.p <- sample(grey.colors(20),2)
          message("Plot colors are randomly chosen since 'col.p' are strings.")
        }else{
          col.p <- grey.colors(20)[col.p]
        }
      }else{
        warning("Argument 'colorspace' is ignored since it is not logical.")
      } 
    }
  }

   
  
  ## ############
  ## check object
  if(missing(object)){    
    stop("Argument 'object' is needed.")  
  }else{
    if (any(substring(class(object)[1],1,7) == "bal.str")){
      stop("Matching is not done before.")
    }else{
      if (is.null(object$bal.test$Stand.diff)){
        stop("No standardized differences are available.")
      }else{
        if (dim(object$bal.test$Stand.diff)[1] != 2)
          stop("Matching is not done before.")
      }
    }
  }
  
  
  ## ####################################
  ## find sel, only names are of interest
  bal.var <- colnames(object$bal.test$Stand.diff)

  if(is.null(sel)){
    name.sel <- bal.var
  }else{
    sel <- find.sel(data=object$data,
                    sel=sel)
    name.sel <- unique(names(sel))

    if (length(intersect(name.sel, bal.var)) == 0){
      stop("No standardized differences available for all selected variables.")
    }else{
      if (length(name.sel) > length(intersect(name.sel, bal.var))){
        cat("There are no standardized differences for selected variable(s):\n")
        print(setdiff(name.sel, bal.var))
        name.sel <- intersect(name.sel, bal.var)
      }
    }
  }

  ## ####
  ## plot
  val.b  <- object$bal.test$Stand.diff[1,][name.sel]
  val.a  <- object$bal.test$Stand.diff[2,][name.sel]
  stdf   <- object$bal.test$Stand.diff[,name.sel]

  ## recode stdf for variables with Inf
  if (any(stdf==Inf)) stdf[stdf==Inf] <- NaN  
 
  
  par(mar=mymar)
  plot(val.b,
       2:(length(name.sel)+1),
       pch=pch.p[1],
       col=col.p[1],
       cex=cex.p,
       axes=F,
       ylim=c(1, length(name.sel)+1.5),
       xlim=c(min(floor(min(round(stdf, 2),na.rm=TRUE)),0),
              ceiling(max(round(stdf, 2), na.rm=TRUE))),
       xlab="",
       ylab="",
       ...)

  points(val.a,
         2:(length(name.sel)+1),
         pch=pch.p[2],
         col=col.p[2],
         cex=cex.p)
  axis(1,
       ...)
  axis(2,
       at=2:(length(name.sel)+1),
       labels=name.sel,
       ...)
  box()

  apply(as.data.frame(cbind(val.a, val.b, c(2:(length(name.sel)+1)))), 1,
      function(x)
      lines(c(x[1],x[2]), c(x[3],x[3]), lty=line.stdf))

  if (plot.alpha)
    abline(v=object$bal.test$alpha, lty=line.alpha)

  if (with.legend){
    if (is.null(legend.xy)){
      legend(x      = ceiling(max(round(stdf, 2), na.rm=TRUE))/2,
             y      = 1.5,
             legend = legend.label,
             horiz  = TRUE,
             pch    = pch.p,
             col    = col.p,
             pt.cex = cex.p,
             cex    = legend.cex)
    }else{      
      if ( !is.numeric(legend.xy) ){
        stop("Argument 'legend.xy' must be numeric.")
      }else{
        if (length(legend.xy) != 2 ){
          stop("Argument 'legend.xy' must be a numeric vector of length 2.")
        }else{          
          legend(x      = legend.xy[1],
                 y      = legend.xy[2],
                 legend = legend.label,
                 horiz  = TRUE,
                 pch    = pch.p,
                 col    = col.p,
                 pt.cex = cex.p,
                 cex    = legend.cex)
        }
      }
    }
  }
  
}


