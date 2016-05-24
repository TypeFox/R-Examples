#' Check that a detection function is monotone
#'
#' Check that a fitted detection function is monotone non-increasing.
#'
#' Evaluates a series of points over the range of the detection function (left to right truncation) then determines:
#'
#' 1. If the detection function is always less than or equal to its value at the left truncation poin (\code{g(x)<=g(left)}, or usually \code{g(x)<=g(0)}).
#' 2. (Optionally) The detection function is always monotone decreasing (\code{g(x[i])<=g(x[i-1])}). This check is only performed when \code{strict=TRUE} (the default).
#' 3. The detection function is never less than 0 (\code{g(x)>=0}).
#' 4. The detection function is never greater than 1 (\code{g(x)<=1}).
#'
#' For models with covarates in the scale parameter of the detection function is evaluated at all observed covariate combinations.
#'
#' Currently covariates in the shape parameter are not supported.
#'
#' @param df a fitted detection function object
#' @param strict if \code{TRUE} (default) the detection function must be "strictly" monotone, that is that (\code{g(x[i])<=g(x[i-1])}) over the whole range (left to right truncation points).
#' @param n.pts number of equally-spaced points between left and right truncation at which to evaluate the detection function (default 100)
#' @param tolerance numerical tolerance for monotonicity checks (default 1e-6)
#' @param plot plot a diagnostic highlighting the non-monotonic areas (default FALSE)
#' @param max.plots when \code{plot=TRUE}, what is the maximum number of plots of non-monotone covariate combinations that should be plotted? Plotted combinations are a random sample of the non-monotonic subset of evaluations. No effect for non-covariate models.
#'
#' @return \code{TRUE} if the detection function is monotone, \code{FALSE} if it's not. \code{warning}s are issued to warn the user that the function is non-monotonic.
#'
#' @keywords utility
#' @author David L. Miller
#' @importFrom stats as.formula model.matrix
#' @importFrom graphics polygon rug
#' @importFrom grDevices rgb
#' @export
check.mono <- function(df,strict=TRUE,n.pts=100,tolerance=1e-6,plot=FALSE,max.plots=6){

  # extract the ddf object from the fitted model
  ddfobj <- df$ds$aux$ddfobj
  # extract the truncation
  right.trunc <- df$meta.data$width
  left.trunc <- df$meta.data$left

  # generate distances between truncation points
  x <- seq(left.trunc,right.trunc,len=n.pts)

  # grab the unique covariate combinations from the data
  if(ddfobj$type=="unif"){
    # no covariates or scale parameter for uniform, so
    # just return a 1x1 matrix of 1
    udat <- matrix(1,1,1)
  }else{
    udat <- unique(model.matrix(as.formula(ddfobj$scale$formula), data=df$dat))
  }

  # function to apply over the unique rows
  chpply <- function(this.udat,ddfobj,x,strict,plot=FALSE){

    # uniform doesn't need a design matrix update, as there is no formula
    if(ddfobj$type!="unif"){
      # build the design matrix for this covariate combination
      this.udat.save <- this.udat
      this.udat <- as.matrix(matrix(this.udat,nrow=1)[rep(1,length(x),by=1),])
      ddfobj$scale$dm <- this.udat

      # dummy data matrix for shape
      if(!is.null(ddfobj$shape)){
        ddfobj$shape$dm <- matrix(1,nrow=length(x),ncol=1)
      }
    }
    # make predictions over the data
    ps <- as.vector(detfct(x,ddfobj,width=right.trunc,standardize=TRUE))

    ## check for monotonicity
    # first: weak monotonicity
    #  check that all evaluations are less than that at
    #  the left truncation point
    weak.diff <- ps[2:length(ps)]-ps[1]
    ps.weak.chk <- weak.diff<=tolerance
    #  maybe some are greater than g(left), but
    #  they might be less than tolerance
    ps.weak.chk[!ps.weak.chk] <- abs(weak.diff[!ps.weak.chk])<=tolerance
    # combine these two checks and issue a warning
    if(!all(ps.weak.chk)){
      warning("Detection function is not weakly monotonic!")
    }

    # second: strict monotonicity
    if(strict){
      # check that all values are less than or equal to the previous
      strict.diff <- diff(ps)
      ps.strict.chk <- strict.diff<=tolerance
      # is the greater than difference less than tolerance?
      ps.strict.chk[!ps.strict.chk] <- abs(strict.diff[!ps.strict.chk])<=tolerance
      # combine these two checks and issue a warning
      if(!all(ps.strict.chk)){
        warning("Detection function is not strictly monotonic!")
      }
    }else{
      ps.strict.chk <- TRUE
    }

    # third: check that the detection function is always in (0,1)
    if(any(ps>1)){
      warning("Detection function is greater than 1 at some distances")
    }
    if(any(ps<0)){
      warning("Detection function is less than 0 at some distances")
    }


    ## if we were asked to provide a diagnostic plot
    if(plot){

      # min and max values of the evaluations
      gmax <- max(c(1,ps))
      gmin <- min(c(0,ps))

      # function to plot the polygons
      # definitely didn't just define this to get a cool function name
      plot.monopoly <- function(chk,gmax,gmin,col){
        start.l <- min(x[which(!chk)])
        end.l <- x[min(max(which(!chk)+2),length(x))]
        polygon(x=c(start.l,start.l,end.l,end.l,start.l),
                y=c(gmin,gmax,gmax,gmin,gmin),col=col,
                lty=0)
      }

      # if we have covariates then put that info in the plot title
      if(ncol(udat)>1){
        cov.info <- paste0("\n",paste0(
                           apply(cbind(names(this.udat.save)[-1],
                                       "=",this.udat.save[-1]),
                                 1,paste0,collapse=""),collapse=","))
      }else{
        cov.info <- ""
      }

      # need to pick out levels here for covariate models
      plot(x, ps, type="l",
           xlim=c(left.trunc,right.trunc),
           ylim=c(gmin,gmax),
           xlab="Distance",ylab="Detection function evaluation",
           main=paste0("Monotonicity diagnostic plot",cov.info))
      rug(x)

      # plot the areas of weak non-monotonicity
      if(!all(ps.weak.chk)){
        plot.monopoly(ps.weak.chk,gmax+1,gmin-1,rgb(0.5,0,0,0.5))
      }

      # plot the areas of strict non-monotonicity
      if(!all(ps.strict.chk)){
        plot.monopoly(ps.strict.chk,gmax+1,gmin-1,rgb(0.5,0,0,0.5))
      }

      # plot detection function >1
      if(any(ps>1)){
        plot.monopoly(rep(FALSE,length(ps)),gmax+1,1,rgb(0.5,0,0,0.5))
      }
      # plot detection function <0
      if(any(ps<0)){
        plot.monopoly(rep(FALSE,length(ps)),0,gmin-1,rgb(0.5,0,0,0.5))
      }
    }

    # return logical -- TRUE == montonic
    return(all(ps.weak.chk) & all(ps.strict.chk) &
                 all(ps<=1) & all(ps>=0))
  }

  # apply the check function to all the unique covariate combinations
  mono.status <- apply(udat,1,chpply,x=x,strict=strict,ddfobj=ddfobj)

  # if plotting was requested and there are non-monotonicity
  if(plot){
    # if no covariates or only 1 unique combination
    if(nrow(udat)==1){
      # re-run doing plotting but not producing the warnings a second time
      d <- suppressMessages(apply(udat,1,chpply,x=x,strict=strict,
                                  ddfobj=ddfobj,plot=TRUE))
    }else{
      if(!all(mono.status)){
        # data frame of non-monotonic covariate combinations
        plot.data <- udat[!mono.status,]
      }else{
        # all plot data
        plot.data <- udat
      }


      # take a sample
      if(is.matrix(plot.data)){
        # might be fewer combinations than max.plots
        max.plots <- min(max.plots,nrow(plot.data))
        plot.sample <- plot.data[sample(1:nrow(plot.data),max.plots),]
      }else{
        # if we have only 1 row
        plot.sample <- matrix(plot.data,nrow=1)
        max.plots <- 1
      }

      # use plot.layout to get the layout
      dd <- plot.layout(1:max.plots,pages=1)

      # make the plots
      d <- suppressMessages(apply(plot.sample,1,chpply,x=x,ddfobj=ddfobj,
                                  strict=strict,plot=TRUE))
    }
  }

  # AND together the per-(covariate combination) montonicity statuses
  # then return the logical
  return(all(mono.status))
}
