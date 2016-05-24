#' Model stability and variable importance plots for glmnet
#'
#'
#' @param mf a fitted 'full' model, the result of a call
#'   to lm or glm.
#' @param nlambda how many penalty values to consider.  Default = 100.
#' @param lambda manually specify the penalty values (optional).
#' @param B number of bootstrap replications
#' @param n.cores number of cores to be used when parallel
#'   processing the bootstrap (Not yet implemented.)
#' @param force.in the names of variables that should be forced
#'   into all estimated models. (Not yet implemented.)
#' @param penalty.factor Separate penalty factors can be applied to each
#'   coefficient. This is a number that multiplies lambda to allow
#'   differential shrinkage. Can be 0 for some variables, which implies
#'   no shrinkage, and that variable is always included in the model.
#'   Default is 1 for all variables (and implicitly infinity for variables
#'   listed in exclude). Note: the penalty factors are internally rescaled
#'   to sum to nvars, and the lambda sequence will reflect this change.
#' @param screen logical, whether or not to perform an initial
#'   screen for outliers.  Highly experimental, use at own risk.
#'   Default = FALSE.
#' @param ... further arguments (currently unused)
#' @details The result of this function is essentially just a
#'   list. The supplied plot method provides a way to visualise the
#'   results.
#' @export
#' @seealso \code{\link{plot.bglmnet}}


bglmnet = function(mf, nlambda = 100, lambda = NULL, B = 100,
                   penalty.factor, screen = FALSE,
                   n.cores = NULL,
                   force.in = NULL) {
  m = mextract(mf, screen = screen)
  fixed = m$fixed
  yname = m$yname
  family = m$family
  fam = family$family
  if(!is.element(fam,c("gaussian", "binomial", "poisson","multinomial", "cox", "mgaussian"))){
    stop(paste("family is",fam,
               "but it needs to be one of gaussian, binomial, poisson, multinomial, cox, mgaussian"),
         call. = FALSE)
  }
  
  Xy = m$X
  kf = m$k
  X = Xy[,1:(kf - 1)]
  Y = Xy[,kf]
  n = m$n
  X = scale(X) * sqrt(n)/sqrt(n - 1)
  #X[which(is.na(X))] = 0
  X = cbind(1, X)
  if (missing(penalty.factor)) {
    # link this with force.in
    penalty.factor = c(0, rep(1, kf))
  }
  if (!is.null(lambda)) {
    nlambda = length(lambda)
  }
  temp = glmnet::glmnet(X, Y, alpha = 1, nlambda = nlambda,
                        lambda = lambda,
                        penalty.factor = penalty.factor,
                        weights = m$wts)
  mat = NULL
  # redefine lambda explicitly
  lambda = temp$lambda
  nlambda = length(lambda)
  compteur = matrix(0, kf, nlambda)
  mfstar = do.call("glm",list(fixed, data = Xy, family = family,weights = m$wts))
  ystar = stats::simulate(object = mfstar, nsim = B)
  #ystar[is.na(ystar)] = Xy[is.na(ystar),yname]
  
  betaboot = array(0,dim = c(kf,nlambda,B))
  rownames(betaboot) = names(mfstar$coef)
  for (j in 1:B) {
    for (i in 1:nlambda) {
      temp = glmnet::glmnet(X, ystar[,j], alpha = 1,
                            lambda = lambda[i],
                            #penalty.factor = penalty.factor,
                            family = fam,
                            weights = m$wts)
      betaboot[,i,j] = (temp$beta[, 1] != 0)
    }
  }
  compteur2 = apply(betaboot,c(1,2),sum)
  probavariable = compteur2/B
  mods = list()
  get.names = function(x) paste(names(x)[x == 1],collapse = "+")
  for (k in 1:length(lambda)) {
    mods[[k]] = table(apply(betaboot[,k,],2,get.names))
  }
  all.mods = unique(names(unlist(mods)))
  all.mods[all.mods == ""] = "1"
  all.ll = rep(0,length(all.mods))
  all.k = rep(0,length(all.mods))
  for (k in 1:length(all.mods)) {
    # don't need to do this for models that include REDUNDANT.VARIABLE
    all.ll[k] = -2*stats::logLik(stats::glm(stats::as.formula(paste(yname,"~",all.mods[k])),
                                            data = Xy,
                                            family = family,
                                            weights = m$wts))
    # number of variables including intercept
    all.k[k] = length(unlist(strsplit(all.mods[[k]],
                                      split = "+",fixed = TRUE))) + 1
  }
  all.k[all.mods=="1"] = 1
  colnames(probavariable) = round(lambda,3)
  mod.sum = data.frame(mod.names = all.mods,ll=all.ll,k=all.k)
  blarout = list(frequency = probavariable,
                 lambda = lambda,
                 mods = mods,
                 mod.sum = mod.sum,
                 screen = screen,
                 vars = names(mf$coef),
                 call = match.call())
  class(blarout) = "bglmnet"
  return(blarout)
}


#' Plot diagnostics for a bglmnet object
#'
#' A plot method to visualise the results of a \code{bglmnet} object.
#'
#' @param x \code{bglmnet} object, the result of \code{\link{bglmnet}}
#' @param highlight the name of a variable that will be highlighted.
#' @param classic logical.  If \code{classic=TRUE} a
#'   base graphics plot is provided instead of a googleVis plot.
#'   Default is \code{classic=FALSE}.
#' @param tag Default NULL. Name tag of the objects to be extracted 
#' from a gvis (googleVis) object. 
#' 
#' The default tag for is NULL, which will 
#' result in R opening a browser window.  Setting \code{tag='chart'} 
#' or setting \code{options(gvis.plot.tag='chart')} is useful when 
#' googleVis is used in scripts, like knitr or rmarkdown. 
#' 
#' @param shiny Default FALSE. Set to TRUE when using in a shiny interface.
#' 
#' @param which a vector specifying the plots to be output. Variable
#'   inclusion type plots \code{which="vip"} or model description loss against
#'   penalty parameter \code{which="boot"}.
#' @param width Width of the googleVis chart canvas area, in pixels.
#'   Default: 800.
#' @param height Height of the googleVis chart canvas area, in pixels.
#'   Default: 400.
#' @param chartWidth googleVis chart area width.
#'   A simple number is a value in pixels;
#'   a string containing a number followed by \code{\%} is a percentage.
#'   Default: \code{"60\%"}
#' @param chartHeight googleVis chart area height.
#'   A simple number is a value in pixels;
#'   a string containing a number followed by \code{\%} is a percentage.
#'   Default: \code{"80\%"}
#' @param fontSize font size used in googleVis chart.  Default: 12.
#' @param left space at left of chart (pixels?).  Default: "50".
#' @param top space at top of chart (pixels?).  Default: "30".
#' @param axisTitlesPosition Where to place the googleVis axis titles,
#'   compared to the chart area. Supported values:
#'   "in" - Draw the axis titles inside the the chart area.
#'   "out" - Draw the axis titles outside the chart area.
#'   "none" - Omit the axis titles.
#' @param dataOpacity The transparency of googleVis data points,
#'   with 1.0 being completely opaque and 0.0 fully transparent.
#' @param options a list to be passed to the googleVis function giving
#'   complete control over the output.  Specifying a value for
#'   \code{options} overwrites all other plotting variables.
#' @param backgroundColor The background colour for the main area
#'   of the chart. A simple HTML color string,
#'   for example: 'red' or '#00cc00'.  Default: 'transparent'
#' @param plb lower bound on the probability of a model being selected. If
#'   a model has a selection probability lower than plb it will not be
#'   plotted.
#' @param hAxis.logScale logical, whether or not to use a log scale on
#'   the horizontal axis. Default = TRUE.
#' @param ... further arguments (currently unused)
#' @export
#' @seealso \code{\link{bglmnet}}


plot.bglmnet = function(x, highlight, classic = FALSE, tag = NULL, shiny = FALSE,
                        which=c("boot","vip"),
                        width=800, height=400, fontSize=12,
                        left=50, top=30,
                        chartWidth="60%",
                        chartHeight="80%",
                        axisTitlesPosition="out",
                        dataOpacity=0.5,
                        options=NULL,
                        hAxis.logScale = TRUE,
                        backgroundColor = 'transparent',
                        plb = 0.01, ...) {
  if (backgroundColor == "transparent") {
    backgroundColor = "{stroke:null, fill:'null', strokeSize: 0}"
  } else {
    backgroundColor = paste("{stroke:null, fill:'",backgroundColor,
                            "', strokeSize: 0}", sep = "")
  }
  B = sum(x$mods[[1]])
  gvis.hAxis = paste("{title:'Penalty parameter',
                    logScale:'",hAxis.logScale,"' ,
                    baseline:",0," ,
                     maxValue:",max(x$lambda)*1.1," ,
                     minValue:",min(x$lambda),"}",sep="")
  
  if("boot"%in%which){
    l.vec = rep(x$lambda,times = lapply(x$mods,length))
    mod.vec = unlist(x$mods)
    mod.names = names(mod.vec)
    mod.names[mod.names==""] = "1"
    mod.vec.counts = as.numeric(mod.vec)
    mod.vec.prob = mod.vec.counts/B
    df.temp = data.frame(l.vec,mod.vec.counts,mod.vec.prob,mod.names)
    # remove redundant variables
    df = df.temp[-grep("REDUNDANT.VARIABLE",df.temp$mod.names),]
    df.full = merge(df,x$mod.sum,all.x = TRUE)
    df.sub = subset(df.full,df.full$mod.vec.prob > plb)
    df.sub$mod.names = as.character(df.sub$mod.names)
    if (classic) {
      warning("Classic plot not implemented.")
      #graphics::plot(df.sub$ll~df.sub$l.vec,cex=df.sub$mod.vec.counts/5,
      #     xlim=c(min(x$lambda),max(x$lambda)))
    }
    if (missing(highlight)) { # highlight best bivariate variable
      no.highlight = TRUE
      if(sum(df.sub$k==2)>0){
        dfk2 = unique(df.sub[df.sub$k==2,c(1,5)])
        highlight = dfk2$mod.names[which.min(dfk2$ll)]
      } else highlight =  x$vars[2]
    }
    
    mod.parts = lapply(df.sub$mod.names,FUN = strsplit,"+",fixed=TRUE)
    find.var = function(x,highlight){
      is.element(highlight,unlist(x))
    }
    var.ident = unlist(lapply(mod.parts,find.var,highlight=highlight))
    var.ident[var.ident==TRUE] =  paste("With",highlight)
    var.ident[var.ident==FALSE] =  paste("Without",highlight)
    df.sub$var.ident = var.ident
    gvis.title = paste("Model stability plot for glmnet",sep="")
    #x.ticks=paste(1:max(x$lk$k),collapse=",")
    chartArea = paste("{left:",left,
                      ",top:",top,
                      ",width:'",chartWidth,
                      "',height:'",chartHeight,"'}",sep="")
    bubble = paste("{opacity:",dataOpacity,
                   ", textStyle: {color: 'none'}}",sep="")
    
    if(is.null(options)){
      use.options=list(title=gvis.title,
                       fontSize = fontSize,
                       vAxis="{title:'-2*Log-likelihood'}",
                       hAxis=gvis.hAxis,
                       sizeAxis = "{minValue: 0, minSize: 1,
                                maxSize: 20, maxValue:1}",
                       axisTitlesPosition=axisTitlesPosition,
                       bubble = bubble,
                       chartArea=chartArea,
                       width=width, height=height,
                       backgroundColor=backgroundColor,
                       explorer= "{axis: 'vertical',
                               keepInBounds: true,
                               maxZoomOut: 1,
                               maxZoomIn: 0.01,
                               actions: ['dragToZoom',
                                         'rightClickToReset']}")
    } else {use.options = options}
    
    fplot = googleVis::gvisBubbleChart(data = df.sub,idvar = "mod.names",xvar = "l.vec",
                                       yvar = "ll", colorvar = "var.ident",
                                       sizevar = "mod.vec.prob",
                                       options = use.options)
    if(shiny){
      return(fplot)
    } else {
      graphics::plot(fplot, tag = tag)
    }
  }
  if ("vip" %in% which) {
    var.names = x$vars[x$vars != "(Intercept)"]
    p.var = t(x$freq)
    p.var = p.var[,colnames(p.var) %in% var.names]
    sortnames = names(sort(apply(p.var, 2, mean), decreasing = TRUE))
    vip.df = p.var[,sortnames]
    vip.df = data.frame(lambda = x$lambda, vip.df)
    #tid = c(1,2,4,6:dim(vip.df)[2])
    #vip.df[, tid] = sapply(vip.df[, tid], as.numeric)
    gvis.title = "Variable inclusion plot"
    chartArea = paste("{left:",left,
                      ",top:",top,
                      ",width:'",chartWidth,
                      "',height:'",chartHeight,"'}", sep = "")
    if (is.null(options)) {
      use.options = list(title = gvis.title,
                       fontSize = fontSize,
                       vAxis = "{title:'Bootstrapped probability'}",
                       hAxis = gvis.hAxis,
                       sizeAxis = "{minValue: 0, minSize: 1,
                       maxSize: 20, maxValue:1}",
                       axisTitlesPosition = axisTitlesPosition,
                       chartArea = chartArea,
                       width = width, height = height,
                       backgroundColor = backgroundColor,
                       annotations = "{style:'line'}")
    } else {use.options = options}
    fplot = googleVis::gvisLineChart(data = vip.df,
                                     xvar = "lambda",
                                     yvar = sortnames,
                                     options = use.options)
    if(shiny){
      return(fplot)
    } else {
      return(graphics::plot(fplot, tag = tag))
    }
  } else return(invisible())
}
