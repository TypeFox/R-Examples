################################################################################
## Functions to run a histogram based naive bayes
################################################################################

## naiveBayes
##==============================================================================
#'  Naive bayes classifier using histograms and shrinkage
#'
#'  Fits a naive bayes model to continous and categorical/factor predictors.
#'    Continous predictors are first binned, then estimates shrunk towards zero.
##  Inputs:
#'  @param formula an object of class \code{\link{formula}} (or one that can be 
#'    coerced to that class): a symbolic description of the model to be fitted. 
#'    Only main effects (not interactions) are allowed. 
#'  @param data data.frame of predictors, can include continuous and
#'    categorical/factors along with a response vector (1 = linked, 0 = unlinked),
#'    and (optionally) observation weights (e.g., \code{weight} column). 
#'    The column names of data need to include the terms specified in \code{formula}.
#'  @param X data frame of categorical and/or numeric variables
#'  @param y binary vector indicating linkage (1 = linked, 0 = unlinked) or 
#'         logical vector (TRUE = linked, FALSE = unlinked)
#'  @param weights a vector of observation weights or the column name in 
#'    \code{data} that corresponds to the weights.
#'  @param df the degrees of freedom for each component density. if vector, each
#'    predictor can use a different df
#'  @param nbins the number of bins for continuous predictors
#'  @param partition for binning; indicates if breaks generated from quantiles
#'    or equal spacing
##  Outputs:
#'  @return BF a bayes factor object; list of component bayes factors
##  Notes:
#' @description After binning, this adds pseudo counts to each bin count to give
#'  df approximate degrees of freedom. If partition=quantile, this does not
#'  assume a continuous uniform prior over support, but rather a discrete uniform
#'  over all (unlabeled) observations points.
#' @seealso \code{\link{predict.naiveBayes}}, \code{\link{plot.naiveBayes}}
#'  @examples
#'  # See vignette: "Statistical Methods for Crime Series Linkage" for usage.  
#'  @export
#'  @name naiveBayes 
##==============================================================================
NULL

## naiveBayes
##==============================================================================
## Formula interface for naive Bayes classifier
#' @rdname naiveBayes
##==============================================================================
naiveBayes <- function(formula,data,weights,df=20,nbins=30,partition=c('quantile','width')){
  if (missing(data)) data = environment(formula)
  mf = match.call(expand.dots = FALSE)
  m = match(c("formula", "data", "weights"), names(mf), 0L)
  mf = mf[c(1L, m)]
  mf[[1L]] = quote(stats::model.frame)
  mf = eval(mf, parent.frame())
  mt = attr(mf, "terms")
  if(any(attr(mt, "order") > 1)) stop("only main effects allowed. Modify formula.")  
  x = mf[, attr(mt,"term.labels"), drop = FALSE]
  y = model.response(mf)
  weights = as.vector(model.weights(mf))
  if (!is.null(weights) && !is.numeric(weights)) 
      stop("'weights' must be a numeric vector")
  if (!is.null(weights) && any(weights < 0)) 
      stop("negative weights not allowed")
  NB = naiveBayes.fit(x,y,weights,df=df,nbins=nbins,partition=partition)
return(NB)
}


## naiveBayes.fit
##==============================================================================
## Direct call to naive bayes classifier
#' @rdname naiveBayes
##==============================================================================
naiveBayes.fit <- function(X,y,weights,df=20,nbins=30,partition=c('quantile','width')){
  partition = match.arg(partition)
  X = data.frame(X)
  y = as.integer(y)
  vars = colnames(X)
  nvars = length(vars)
  df = rep(df,length=nvars)
  BF = vector("list",nvars)
  if(missing(weights)) weights = rep(1,length(y))
  for(j in 1:nvars){
    var = vars[j]
    x = X[,var,drop=TRUE]
    if(class(x) %in% c('numeric','integer')){
      bks = make.breaks(x,partition,nbins=nbins)
    } else  bks = NULL
    BF[[j]] = getBF(x,y,weights,breaks=bks,df=df[j]) 
  }
  names(BF) = vars
  class(BF) = "naiveBayes"
return(BF)
}



## predict.naiveBayes
##==============================================================================
#'  Generate prediction (sum of log bayes factors) from a \code{naiveBayes} object
##  Inputs:
#'  @param object a naive bayes object from \code{\link{naiveBayes}}
#'  @param newdata data frame of new predictors, column names must match NB names
#'  @param components (logical) return the log bayes factors from each component
#'    or return the sum of log bayes factors
#'  @param vars the names or column numbers of specific predictors. If NULL, then
#'    all predictors will be used
#'  @param \ldots not currently used    
##  Outputs:
#'  @return BF if \code{components = FALSE}, the sum of log bayes factors, if
#'    \code{components = TRUE} the component bayes factors (useful for plotting). 
#'    
#'    It will give a warning, but still produce output if X is missing predictors. 
#'    The output in this situation will be based on the predictors that are in X. 
##  Notes:
#'  @description This does not include the log prior odds, so will be off by a
#'    constant. 
#'  @seealso \code{\link{naiveBayes}}, \code{\link{plot.naiveBayes}}  
#'  @examples
#'  # See vignette: "Statistical Methods for Crime Series Linkage" for usage.     
#'  @export
##==============================================================================
predict.naiveBayes <- function(object,newdata,components=FALSE,vars=NULL,...){
  if(missing(newdata)) stop("newdata must be provided")
  X = data.frame(newdata)
  NB = object
  if(is.null(vars))  vars = names(NB)
  missing.vars = !(vars %in% colnames(X))
  if(any(missing.vars)) 
    warning('The columns: ',paste(vars[missing.vars],collapse=','),' are missing from newdata')
  vars = vars[!missing.vars]
  nvars = length(vars)
  BF = matrix(NA,nrow(X),nvars)
  for(j in 1:nvars){
    var = vars[j]
    BF[,j] = predictBF(NB[[var]],X[,var],log=TRUE)
  }
  colnames(BF) = vars
  if(components) return(BF)
  BF = rowSums(BF)
  return(BF)
}




## getBF
##==============================================================================
#'  Estimates the bayes factor for continous and categorical predictors.
#'
#'  Continous predictors are first binned, then estimates shrunk towards zero.
##  Inputs:
#'  @param x predictor vector (continuous or categorical/factors)
#'  @param y binary vector indicating linkage (1 = linked, 0 = unlinked) or 
#'    logical vector (TRUE = linked, FALSE = unlinked)  
#'  @param weights a vector of observation weights or the column name in 
#'    \code{data} that corresponds to the weights. 
#'  @param breaks set of break point for continuous predictors or NULL for
#'    categorical or discrete
#'  @param df the effective degrees of freedom for the cetegorical density
#'    estimates
##  Outputs:
#'  @return data.frame containing the levels/categories with estimated Bayes factor
##  Notes:
#'  @description This adds pseudo counts to each bin count to give df effective
#'    degrees of freedom. Must have all possible factor levels and must be of
#'    factor class.
#'  @note Give linked and unlinked a different prior according to sample size
#'  @examples
#'  # See vignette: "Statistical Methods for Crime Series Linkage" for usage.     
#'  @export
##==============================================================================
getBF <- function(x,y,weights,breaks=NULL,df=5){
  replaceNA <- function(x,r=0) as.numeric(ifelse(is.na(x),r,x))
  if(is.data.frame(x)) stop("x must be a vector, not data frame")
  linked = (as.integer(y)==1L)
  var.type = class(x)      
  if(missing(weights)) weights = rep(1,length(linked))
  if(var.type %in% c('numeric','integer')){
    n.bks = length(breaks)
    x = cut(x,breaks=breaks,include.lowest=TRUE)
    x = addNA(x,ifany=TRUE)
    t.linked = tapply(weights[linked],x[linked],sum)
    t.unlinked = tapply(weights[!linked],x[!linked],sum)
    fromto = data.frame(from=c(breaks[-n.bks]),to=c(breaks[-1]))
    if(nrow(fromto)<length(t.linked)) fromto = rbind(fromto,c(NA,NA))
    E = data.frame(fromto,
                   value=levels(x),
                   N.linked=replaceNA(t.linked),
                   N.unlinked=replaceNA(t.unlinked))
        var.type = 'numeric'
  }
  if(var.type %in% c('factor','character')){
    x = as.factor(x)
    x = addNA(x,ifany=TRUE)
    t.linked = tapply(weights[linked],x[linked],sum)
    t.unlinked = tapply(weights[!linked],x[!linked],sum)
    E = data.frame(value=names(t.linked),
                   N.linked=replaceNA(t.linked),
                   N.unlinked=replaceNA(t.unlinked))
    attr(E,'levels') = levels(x)
    var.type = 'categorical'
  }
  #- get 'a' based on degrees of freedom
  df2a <- function(df,k,N)  (N/k) * ((k-df)/(df-1)) # k is # levels
  a2df <- function(a,k,N) k*(N+a)/(N+k*a)
  nlevs = nrow(E)
  df = min(df,nlevs-1e-8)
  a.linked = df2a(df,k=nlevs,N=sum(E$N.linked))
  a.unlinked = df2a(df,k=nlevs,N=sum(E$N.unlinked))

  getP <- function(N,a) (N+a)/sum(N+a)
  E$p.linked = getP(E$N.linked,a.linked)
  E$p.unlinked = getP(E$N.unlinked,a.unlinked)
  E$BF = E$p.linked/E$p.unlinked
  E[is.na(E$BF),'BF'] = 1           # Set 0/0=1
  attr(E,'breaks') = breaks         # add breaks (or NULL)
  attr(E,'a') = c(linked=a.linked,unlinked=a.unlinked) # add shrinkage parameters
  attr(E,'df') = df
  attr(E,'df2') = c(linked=a2df(a.linked,k=nlevs,N=sum(E$N.linked)),
                                unlinked=a2df(a.unlinked,k=nlevs,N=sum(E$N.unlinked)))
  attr(E,'type') = var.type
  return(E)
}





## predictBF
##==============================================================================
#'  Generate prediction of a component bayes factor
##  Inputs:
#'  @param BF bayes factor data.frame from \code{\link{getBF}}
#'  @param x vector of new predictor values
#'  @param log (logical) if \code{TRUE}, return the \bold{log} bayes factor estimate
##  Outputs:
#'  @return estimated (log) bayes factor from a single predictor
#'  @description This does not include the log prior odds, so will be off by a constant
#'  @examples
#'  # See vignette: "Statistical Methods for Crime Series Linkage" for usage.     
#'  @export
##==============================================================================
predictBF <- function(BF,x,log=TRUE){
  breaks = attr(BF,'breaks')
  if(!is.null(breaks)){
    x = cut(x,breaks=breaks,include.lowest=TRUE)
  }
  ind = match(as.character(x),as.character(BF$value))
  bf = BF$BF[ind]
  bf[is.na(bf)] = 1    # if data in x is outside of breaks, then set bf to 1
  if(log) bf = log(bf)
  attr(bf,'log') = log
  return(bf)
}


## make.breaks
##==============================================================================
#'  Make break points for binning continuous predictors
##  Inputs:
#'  @param x observed sample
#'  @param type one of \code{width} (fixed width) or \code{quantile} binning
#'  @param nbins number of bins
#'  @param binwidth bin width; corresponds to quantiles if type='quantile'
##  Outputs:
#'  @return set of unique break points for binning
#'  @keywords internal
##==============================================================================
make.breaks <- function(x,type='quantile',nbins=NULL,binwidth=NULL){
  type <- match.arg(type, c("quantile", "width"))
  if ((!is.null(nbins) && !is.null(binwidth)) || (is.null(nbins) &&
                                                    is.null(binwidth))) {
    stop("Specify exactly one of nbins or width")
  }
  if (type == "width") {
    rng <- range(x, na.rm = TRUE, finite = TRUE)
    if (!is.null(binwidth)) bks = unique(c(rng[1],seq(rng[1],rng[2],by=binwidth),rng[2]))
    else                    bks = seq(rng[1], rng[2], length = nbins + 1)
  }
  if (type == "quantile") {
    if (!is.null(binwidth)) probs <- seq(0, 1, by = binwidth)
    else                    probs <- seq(0, 1, length = nbins + 1)
    bks = quantile(x, probs, na.rm = TRUE)
  }
  return(sort(unique(bks)))
}


#==============================================================================#
# Plotting functions
#==============================================================================#

## plotBKG
##==============================================================================
#'  Generate a background plot
#'
#'  This facilitates a common plot background
#'  @param xlim range of x-axis
#'  @param ylim range of y-axis
#'  @param background (logical) should a background color be used
#'  @param x.minor values for minor axis
#'  @param x.major values for major axis
#'  @param grid.lines (logical) should grid lines be added to plot
#'  @param glwd1 linewidth of gridlines
#'  @param glwd2 linewidth of gridlines
#'  @param bkg.col color of plot background
#'  @param grid.col color of grid lines
#'  @param boxed (logical) should plot be boxed
#'  @param \ldots other arguments passed to plot
#'  @return makes a plot
#'  @keywords internal
##==============================================================================
plotBKG <- function(xlim,ylim,background=TRUE,x.minor,x.major,grid.lines=TRUE,
                    glwd1=2,glwd2=1,bkg.col='grey90',grid.col='white',
                    boxed=TRUE,...){
  plot(xlim,ylim,typ='n',
       ylab='',xlab='',
       col.axis='grey50',cex.axis=.8,tcl=-.35,las=1,...)
  if(background){    # Add background color
    rng = par('usr')
    if(par('ylog')) rng[c(3,4)] = 10^(rng[c(3,4)])
    if(par('xlog')) rng[c(1,2)] = 10^(rng[c(1,2)])
    rect(rng[1],rng[3],rng[2],rng[4],col=bkg.col,border="transparent")
  }
  if(grid.lines){
    abline(h=axTicks(2),col = grid.col,lty=1,lwd=glwd1)
    if(missing(x.major))  abline(v=axTicks(1),col = grid.col,lty=1,lwd=glwd1)
    else                  abline(v=x.major   ,col = grid.col,lty=1,lwd=glwd1)
    get.minor <- function(side=c(1,2)){  # return equal spaced minor axis
      xt = axTicks(side)
      nt = length(xt)
      axlog = ifelse(side==1,par('xlog'),par('ylog'))
      if(axlog) xt = log10(xt)
      minor = diff(xt[1:2])/2 + xt
      if(axlog) minor = 10^minor
    return(minor)
    }
    if(missing(x.minor)) abline(v=get.minor(1),col = grid.col,lty=1,lwd=glwd2)
    else abline(v=x.minor,col=grid.col,lty=1,lwd=glwd2)
    abline(h=get.minor(2),col = grid.col,lty=1,lwd=glwd2)
    if(boxed) box(col='grey50')
  }
}


## plot.naiveBayes
##==============================================================================
#'  Plots for Naive Bayes Model
#'  
#'  Plots (component) bayes factors from naiveBayes()
#'  @param x a \code{\link{naiveBayes}} object
#'  @param vars name or index of naive Bayes components to plot. Will plot all 
#'    if blank.
#'  @param log.scale (logical)
#'  @param show.legend either a value or values indicating which plot to show 
#'    the legend, or TRUE/FALSE to show or not show the legend on all plots.
#'  @param cols Colors for plotting. First element is for linkage, second unlinked
#'  @param \ldots arguemnts passed into \code{\link{plotBF}}
#'  @return plots of Bayes factor from a naive Bayes model
#'  @description This function attempts to plot all of the component plots in 
#'    one window by using the mfrow argument of par. If more control is desired
#'    then use \code{\link{plotBF}} to plot individual Bayes factors.
#'  @seealso \code{\link{plotBF}}, \code{\link{naiveBayes}}, 
#'    \code{\link{predict.naiveBayes}}  
#'  @examples
#'  # See vignette: "Statistical Methods for Crime Series Linkage" for usage.     
#'  @export
##==============================================================================
plot.naiveBayes <- function(x,vars,log.scale=TRUE,show.legend=1,
  cols = c(color('darkred',alpha=.75),color('darkblue',alpha=.75)),...){
  if(missing(vars)) vars = names(NB)
  NB = x
  vars = intersect(vars,names(NB))
  nvars = length(vars)
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  if(nvars > 1){
    par(mar=c(2,4,2,2),mgp = c(3, .75, 0))
    if(nvars <= 2)       mfrow=c(1,2)
    else if(nvars <= 4)  mfrow=c(2,2)
    else if(nvars <= 6)  mfrow=c(2,3)
    else if(nvars <= 9)  mfrow=c(3,3)
    else if(nvars <= 12) mfrow=c(3,4)
    else {mfrow=c(4,4);  warning("Too many components; use plotBF()")}
    par(mfrow=mfrow)
  }
  for(j in 1:nvars){
    var = vars[j]
    show.leg = ( j %in% show.legend || isTRUE(show.legend))
    plotBF(NB[[var]],log.scale=log.scale,show.legend=show.leg,cols=cols,...)
    title(paste(var))
  }
invisible()  
}

## plotBF
##==============================================================================
#'  plots 1D bayes factor
#'  @param BF Bayes Factor
#'  @param log.scale (logical)
#'  @param show.legend (logical)
#'  @param xlim range of x-axis
#'  @param ylim range of y-axis
#'  @param cols Colors for plotting. First element is for linkage, second unlinked
#'  @param \ldots arguemnts passed into \code{\link{plotBKG}}
#'  @return plot of Bayes factor
#'  @seealso \code{\link{plot.naiveBayes}}, \code{\link{plotBKG}}
#'  @examples
#'  # See vignette: "Statistical Methods for Crime Series Linkage" for usage.     
#'  @export
##==============================================================================
plotBF <- function(BF,log.scale=TRUE,show.legend=TRUE,xlim,ylim=NULL,
  cols = c(color('darkred',alpha=.75),color('darkblue',alpha=.75)),...){
  if(missing(ylim) || is.null(ylim)){
    ylim = range(BF$BF,na.rm=TRUE,finite=TRUE)
    if(log.scale) ylim = c(-1,1)*min(12,max(abs(log(ylim))))
  }
  if(attr(BF,'type') %in% "numeric"){
    BF = BF[!is.na(BF$to),]     # Remove NAs
    n = nrow(BF)
    xx = c(BF$from[1],BF$to)
    yy = c(BF$BF,BF$BF[n]) # Add extra line for step plotting
    if(log.scale) yy = log(yy)
    if(missing(xlim))  xlim = range(xx,na.rm=TRUE,finite=TRUE)
    plotBKG(xlim,ylim,...)
    title(ylab=ifelse(log.scale,'log(BF)','BF'))
    baseline = ifelse(log.scale,0,1)
    segments(xlim[1],baseline,xlim[2],baseline)
    y.thres = ifelse(log.scale,0,1)
    rect(xx[-(n+1)], y.thres, xx[-1], yy[-(n+1)],
      col=ifelse(yy>y.thres,cols[1],cols[2]),border=NA)
  }
  if(attr(BF,'type') %in% "categorical"){
    if(nrow(BF)>10) BF = rankBF(BF,n=10)
    else BF = transform(BF,logBF=log(BF))
    mp = barplot(BF$logBF,names.arg=BF$value,plot=FALSE)
    xadd = diff(mp)[1]/2 # .2
    plotBKG(c(min(mp)-xadd,max(mp)+xadd),ylim,xaxt='n',yaxt='n',x.major=NULL,x.minor=mp,...)
    barplot(BF$logBF,names.arg=BF$value,
            cex.names=.9,
            col=ifelse(BF$logBF>0,cols[1],cols[2]),
            col.axis='grey50',cex.axis=.8,tcl=-.25,
            ylim=ylim, las=1, add=TRUE)
    title(ylab=ifelse(log.scale,'log(BF)','BF'))
  }
  if(show.legend){
    leg.names = c(expression(paste('Favors ',H[L])),expression(paste('Favors ',H[U])))
    legend('topright',leg.names,#col=1:2,lwd=2,
       fill = cols,
       bty = 'n',border = NA)
  }
}

## rankBF
##==============================================================================
##  Orders Category levels according to absolute value log bayes factor
##==============================================================================
rankBF <- function(BF,n=10,thres=NULL){
  n = min(n,nrow(BF))
  logBF = log(BF$BF)  # log-Bayes factor
  ord = order(abs(logBF),decreasing=TRUE)
  if(is.null(thres)) ind = ord[1:n]
  else ind = ord[abs(logBF[ord])>=thres]
  data.frame(BF[ind,],logBF=logBF[ind])
}

## color
##==============================================================================
#'  Creates transparent colors
#'  @param col Color that \R recognizes (names or number)
#'  @param alpha transparency value
#'  @seealso \code{\link{col2rgb}}
#'  @importFrom grDevices col2rgb rgb 
#'  @keywords internal
##==============================================================================
color <- function(col,alpha=NULL){
  col = col2rgb(col,alpha=TRUE)/255
  if(!is.null(alpha) &  alpha>=0 & alpha <=1){
    col[4] = alpha
  }
  rgb.col = rgb(col[1],col[2],col[3],alpha=col[4])
return(rgb.col)
}





