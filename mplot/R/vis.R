#' Model stability and variable inclusion plots
#'
#' Calculates and provides the plot methods for standard
#' and bootstrap enhanced model stability plots (\code{lvk} and
#' \code{boot}) as well as variable inclusion plots (\code{vip}).
#'
#' @param mf a fitted 'full' model, the result of a call
#'   to lm or glm (and in the future lme or lmer)
#' @param nvmax size of the largest model that can still be
#'   considered as a viable candidate
#' @param B number of bootstrap replications
#' @param lambda.max maximum penalty value for the vip plot,
#'   defaults to 2*log(n)
#' @param nbest maximum number of models at each model size
#'   that will be considered for the lvk plot. Can also take
#'   a value of \code{"all"} which displays all models.
#' @param n.cores number of cores to be used when parallel
#'   processing the bootstrap
#' @param force.in the names of variables that should be forced
#'   into all estimated models. (Not yet implemented.)
#' @param screen logical, whether or not to perform an initial
#'   screen for outliers.  Highly experimental, use at own risk.
#'   Default = \code{FALSE}.
#' @param redundant logical, whether or not to add a redundant
#'   variable.  Default = \code{TRUE}.
#' @param ... further arguments (currently unused)
#' @details The result of this function is essentially just a
#'   list. The supplied plot method provides a way to visualise the
#'   results.
#' @seealso \code{\link{plot.vis}}
#' @references Mueller, S. and Welsh, A. H. (2010), On model
#'   selection curves. International Statistical Review, 78:240-256.
#'   doi: 10.1111/j.1751-5823.2010.00108.x
#'
#'   Murray, K., Heritier, S. and Mueller, S. (2013), Graphical
#'   tools for model selection in generalized linear models.
#'   Statistics in Medicine, 32:4438-4451. doi: 10.1002/sim.5855
#' @export
#' @import foreach
#' @import parallel
#' @examples
#' n = 100
#' set.seed(11)
#' e = rnorm(n)
#' x1 = rnorm(n)
#' x2 = rnorm(n)
#' x3 = x1^2
#' x4 = x2^2
#' x5 = x1*x2
#' y = 1 + x1 + x2 + e
#' dat = data.frame(y,x1,x2,x3,x4,x5)
#' lm1 = lm(y~.,data=dat)
#' \dontrun{
#' v1 = vis(lm1)
#' plot(v1,highlight="x1")
#' }

vis=function(mf, nvmax, B=100, lambda.max, nbest=5,
             n.cores, force.in=NULL, screen=FALSE,
             redundant=TRUE,...){
  
  m = mextract(mf,screen=screen,redundant=redundant)
  fixed = m$fixed
  yname = m$yname
  family = m$family
  X = m$X
  kf = m$k
  n = m$n
  initial.weights = m$wts
  if(missing(nvmax)) nvmax = kf
  if(nbest=="all") {
    nbest = max(choose((kf-1),0:(kf-1)))
  }
  if(!is.numeric(nbest)){
    stop("nbest should be numeric or 'all'")
  }
  
  add.intercept.row = function(em,rs.which,rs.stats){
    # add an intercept row
    rs.which = rbind(c(1,rep(0,dim(rs.which)[2]-1)),rs.which)
    # used for testing
    # f0 = stats::as.formula(paste(yname," ~ 1"))
    # lm0 = stats::lm(formula=f0,data=X)
    i=1 # special case for the null model
    intercept = TRUE
    n1 = em$nn-intercept
    sigma2 = em$sserr/(n1+intercept-em$last)
    nullrss = em$nullrss
    ress = em$nullrss
    # ress = sum(resid(lm0)^2) # used for testing
    vr = ress/nullrss
    rsq = 1-vr
    adjr2 = 1-vr*n1/(n1+intercept-i)
    cp = ress/sigma2-(n1+intercept-2*i)
    # the way leaps calculates the bic:
    bic = (n1+intercept)*log(vr)+i*log(n1+intercept)
    logLikelihood = -(em$nn + em$nn*log(2*pi) + em$nn*log(em$nullrss/em$nn))/2
    rs.stats = rbind(c(logLikelihood,ress,bic,cp,rsq,adjr2,i),rs.stats)
    return(cbind(rs.which,rs.stats))
  }
  
  ## Initial single pass
  ## (gives the minimum envelopping set of models)
  if (any(class(mf) == "glm") == TRUE) {
    em = bestglm::bestglm(Xy = X,
                          family = family,
                          IC = "BIC",
                          TopModels = 1,
                          nvmax = nvmax)
    # starts with intercept row
    rs.which = em$Subsets[,1:kf] + 0
    rs.stats = em$Subsets[, -c(1:kf)]
    k = rowSums(rs.which)
    rs.all = cbind(rs.which, rs.stats, k)
    # in bestglm rs.all$logLikelihood comes from
    # stats::logLik(model) unless Gaussian in which case
    # -(n/2) * log(sum(resid(ans)^2)/n) is used
    # note the bic in bestglm is calculated as:
    # -2*rs.all$logLikelihood + log(n)*(rs.all$k-1)
  } else {
    em = leaps::regsubsets(x = fixed,
                           data = X,
                           nbest = nbest,
                           nvmax = nvmax,
                           intercept = TRUE,
                           force.in = force.in,
                           really.big = TRUE)
    rs = summary(em)
    # does not start with intercept row
    rs.which = data.frame(rs$which + 0,row.names = NULL)
    k = rowSums(rs.which)
    # assuming Gaussian errors:
    logLikelihood = -(n + n*log(2*pi) + n*log(rs$rss/n))/2
    rs.stats = cbind(logLikelihood,rs$rss,rs$bic,
                     rs$cp,rs$rsq,rs$adjr2,k)
    colnames(rs.stats) = c("logLikelihood","rss","bic",
                           "cp","rsq","adjr2","k")
    rs.all = add.intercept.row(em,rs.which,rs.stats)
    # note that the BIC in leaps (and add.intercept.row funtion)
    # differs from the bestglm BIC buy a constant
  }
  
  nms = colnames(rs.all)[2:kf]
  nm = apply(rs.all[,2:kf],1,function(x) {
    Reduce(paste,paste(nms[x == 1],sep = "",collapse = "+"))
  })
  nm[nm == ""] = "1"
  nm = paste(yname,"~",nm,sep = "")
  nm = gsub(pattern = " ",replacement = "",x = nm)
  nm = gsub(pattern = "REDUNDANT.VARIABLE","RV",x = nm)
  rs.all$name = nm
  res.single.pass = rs.all
  
  ## Bootstrap version:
  if (missing(n.cores)) n.cores = max(detectCores() - 1, 1)
  cl.visB = makeCluster(n.cores)
  doParallel::registerDoParallel(cl.visB)
  res = foreach(b = 1:B, .packages = c("bestglm")) %dopar% {
    wts = stats::rexp(n = n, rate = 1)
    if (any(class(mf) == "glm") == TRUE) {
      em = bestglm::bestglm(Xy = X,
                            family = family,
                            IC = "BIC",
                            TopModels = 1,
                            weights = wts,
                            nvmax = nvmax)
      # starts with intercept row
      rs.which = em$Subsets[, 1:kf] + 0
      rs.stats = em$Subsets[, -c(1:kf)]
      k = rowSums(rs.which)
      rs.all = cbind(rs.which, rs.stats, k)
      # in bestglm rs.all$logLikelihood comes from
      # stats::logLik(model) unless Gaussian in which case
      # -(n/2) * log(sum(resid(ans)^2)/n) is used
      # note the bic in bestglm is calculated as:
      # -2*rs.all$logLikelihood + log(n)*(rs.all$k-1)
    } else {
      em = leaps::regsubsets(x = fixed,
                             data = X,
                             nbest = nbest,
                             nvmax = nvmax,
                             intercept = TRUE,
                             force.in = force.in,
                             weights = wts,
                             really.big = TRUE)
      rs = summary(em)
      # does not start with intercept row
      rs.which = data.frame(rs$which + 0,row.names = NULL)
      k = rowSums(rs.which)
      # assuming Gaussian errors:
      logLikelihood = -(n + n*log(2*pi) + n*log(rs$rss/n))/2
      rs.stats = cbind(logLikelihood,rs$rss,rs$bic,
                       rs$cp,rs$rsq,rs$adjr2,k)
      colnames(rs.stats) = c("logLikelihood","rss","bic",
                             "cp","rsq","adjr2","k")
      rs.all = add.intercept.row(em,rs.which,rs.stats)
      # note that the BIC in leaps (and add.intercept.row funtion)
      # differs from the bestglm BIC buy a constant
    }
  }
  stopCluster(cl.visB)
  
  ### Variable inclusion Plot Calculations
  if (missing(lambda.max)) lambda.max = 2*log(n)
  lambdas = seq(0,lambda.max,0.01)
  var.in = matrix(NA,ncol=kf,nrow=length(lambdas))
  colnames(var.in) = colnames(res[[1]][1:kf])
  rownames(var.in) = lambdas
  for(i in 1:length(lambdas)){
    resl = lapply(res, function(x) -2*x$logLikelihood+lambdas[i]*x$k)
    min.pos = unlist(lapply(resl,which.min))
    temp.best = mapply(function(x,row) {x[row,1:kf]}, res, min.pos, SIMPLIFY=TRUE)
    temp.best = matrix(as.numeric(temp.best),nrow=dim(temp.best)[1],dimnames=dimnames(temp.best))
    var.in[i,] = rowSums(temp.best)
  }
  
  #### lvk where bubbles reflect frequencey of choice
  best.within.size = function(x){
    rankings.within.size = unlist(stats::aggregate(-2*x$logLikelihood,rank,by=list(x$k))$x)
    return(x[rankings.within.size==1,])
  }
  res.best = lapply(res,best.within.size)
  
  res.df = do.call(rbind.data.frame,res.best)
  res.df = plyr::count(df = res.df,vars = 1:kf)
  res.df$logLikelihood = NA
  res.df$name = NA
  res.df$k = rowSums(res.df[,1:kf])
  
  # iterate over all required models on the original date
  # to obtain logLiks (in future could extract other stats)
  for(i in 1:dim(res.df)[1]){
    ff = paste(yname," ~ ",
               paste(colnames(res.df[2:kf])[res.df[i,2:kf]==1],collapse="+"),sep="")
    if(ff == paste(yname," ~ ",sep="")){
      ff = stats::as.formula(paste(yname,"~1"))
    } else {
      ff = stats::as.formula(ff)
    }
    if(any(class(mf)=="glm")==TRUE){
      em = stats::glm(formula=ff, data=X, family=family,weights=initial.weights)
    } else {
      em = stats::lm(formula=ff, data=X,weights=initial.weights)
    }
    res.df$logLikelihood[i] = as.numeric(stats::logLik(em))
    nm = gsub(pattern=" ",replacement = "",x=Reduce(paste,deparse(ff)))
    nm = gsub(pattern="REDUNDANT.VARIABLE","RV",x=nm)
    res.df$name[i] = nm
  }
  
  res.df = res.df[with(res.df,order(k,-freq)),]
  
  output = list(res.df = res.df,
                res.single.pass = res.single.pass,
                var.in = var.in,
                screen = screen,
                mextract = m,
                lambdas = lambdas,
                B=B)
  
  class(output) = "vis"
  return(output)
}




#' Plot diagnostics for a vis object
#'
#' A plot method to visualise the results of a \code{vis} object.
#'
#' @param x \code{vis} object, the result of \code{\link{vis}}
#' @param highlight the name of a variable that will be highlighted
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
#' @param which a vector specifying the plots to be output.  Variable
#'   inclusion plots \code{which="vip"}; description loss against model
#'   size \code{which="lvk"}; bootstrapped description loss against
#'   model size \code{which="boot"}.
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
#'   for example: 'red' or '#00cc00'.  Default: 'null' (there is an
#'   issue with GoogleCharts when setting 'transparent' related to the
#'   zoom window sticking - once that's sorted out, the default
#'   will change back to 'transparent')
#' @param text logical, whether or not to add text labels to classic
#'   boot plot. Default = \code{FALSE}.
#' @param min.prob when \code{text=TRUE}, a lower bound on the probability of
#'   selection before a text label is shown.
#' @param srt when \code{text=TRUE}, the angle of rotation for the text labels.
#'   Default = -30.
#' @param print.full.model logical, when \code{text=TRUE} this determines if the full
#'   model gets a label or not.  Default=\code{FALSE}.
#' @param max.circle  circles are scaled to make largest dimension this size in inches.
#'   Default = 0.35.
#' @param jitterk amount of jittering of the model size in the lvk and boot plots.
#'   Default = 0.1.
#' @param ylim the y limits of the lvk and boot plots.
#' @param ... further arguments (currently unused)
#' @seealso \code{\link{vis}}
#' @references Mueller, S. and Welsh, A. H. (2010), On model
#'   selection curves. International Statistical Review, 78:240-256.
#'   doi: 10.1111/j.1751-5823.2010.00108.x
#'
#'   Murray, K., Heritier, S. and Mueller, S. (2013), Graphical
#'   tools for model selection in generalized linear models.
#'   Statistics in Medicine, 32:4438-4451. doi: 10.1002/sim.5855
#' @export
#' @examples
#' n = 100
#' set.seed(11)
#' e = rnorm(n)
#' x1 = rnorm(n)
#' x2 = rnorm(n)
#' x3 = x1^2
#' x4 = x2^2
#' x5 = x1*x2
#' y = 1 + x1 + x2 + e
#' dat = data.frame(y,x1,x2,x3,x4,x5)
#' lm1 = lm(y~.,data=dat)
#' \dontrun{
#' v1 = vis(lm1)
#' plot(v1,highlight="x1",which="lvk")
#' }

plot.vis = function(x, highlight, classic = FALSE, tag = NULL, shiny = FALSE,
                    which = c("vip","lvk","boot"),
                    width = 800, height = 400, fontSize = 12,
                    left = 50, top = 30, chartWidth = "60%", chartHeight = "80%",
                    axisTitlesPosition = "out", dataOpacity = 0.5,
                    options=NULL, ylim,
                    backgroundColor = 'transparent',
                    text=FALSE, min.prob = 0.4, srt = -30, max.circle = 0.35,
                    print.full.model = FALSE, jitterk=0.1, ...){
  if (backgroundColor == "transparent") {
    backgroundColor = "{stroke:null, fill:'null', strokeSize: 0}"
  } else {
    backgroundColor = paste("{stroke:null, fill:'",backgroundColor,
                            "', strokeSize: 0}",sep = "")
  }
  
  if (missing(highlight)) {
    # highlight first variable in the coefficient list
    highlight = x$m$exp.vars[1]
    vars = make.names(x$m$exp.vars)
  } else {
    vars = highlight
  }
  if ("lvk" %in% which) {
    var.ident = n.var.ident = NA
    m2ll = -2*x$res.single.pass$logLikelihood
    spk = x$res.single.pass$k
    jitter = stats::runif(length(spk),0 - jitterk,0 + jitterk)
    spk = spk + jitter
    if (classic) {
      for (i in 1:length(vars)) {
        col_high = grDevices::rgb(1, 0, 0, alpha = 0)
        col_nohigh = grDevices::rgb(0, 0, 1, alpha = 0.5)
        colbg = grDevices::rgb(1, 0, 0, alpha = 0.5)
        var.ident = which(x$res.single.pass[,vars[i]] == 1)
        n.var.ident = which(x$res.single.pass[,vars[i]] == 0)
        graphics::par(mar = c(3.4,3.4,0.1,0.1),mgp = c(2.0, 0.75, 0))
        if (missing(ylim)) ylim = c(min(m2ll),max(m2ll))
        graphics::plot(m2ll[n.var.ident] ~ spk[n.var.ident],
                       pch = 19, cex = 1.3, col = col_nohigh,
                       bg = colbg,
                       xlab = "Number of parameters",
                       ylab = "-2*Log-likelihood",
                       ylim = ylim,
                       xlim = c(min(spk),max(spk)))
        graphics::points(m2ll[var.ident] ~ spk[var.ident],
                         pch = 24, bg = colbg,
                         col = col_high, cex = 1.2)
        graphics::legend("topright",legend = c(paste("With",vars[i]),paste("Without",vars[i])),
                         col = c(col_high,col_nohigh),
                         pt.bg = colbg, pch = c(24,19))
        if (length(vars) > 1) {
          graphics::par(ask = TRUE)
        }
      }
    } else {# googleVis version
      var.ident = which(x$res.single.pass[,vars[1]] == 1)
      n.var.ident = which(x$res.single.pass[,vars[1]] == 0)
      with.var = without.var = rep(NA,dim(x$res.single.pass)[1])
      with.var[var.ident] = m2ll[var.ident]
      without.var[n.var.ident] = m2ll[n.var.ident]
      mods = x$res.single.pass$name
      dat = data.frame(k = spk,
                       without.var, without.var.html.tooltip = mods,
                       with.var, with.var.html.tooltip = mods)
      colnames(dat)[4] = paste("With",highlight)
      colnames(dat)[2] = paste("Without",highlight)
      gvis.title = paste("Model stability plot",sep = "")
      x.ticks = paste(1:max(spk), collapse = ",")
      gvis.hAxis = paste("{title:'Number of parameters', ticks: [",
                         x.ticks,"]}")
      chartArea = paste("{left:",left,
                        ",top:",top,
                        ",width:'",chartWidth,
                        "',height:'",chartHeight,"'}",sep = "")
      if (is.null(options)) {
        use.options = list(title = gvis.title,
                           fontSize = fontSize,
                           vAxis = "{title:'-2*Log-likelihood'}",
                           hAxis = gvis.hAxis,
                           axisTitlesPosition = axisTitlesPosition,
                           chartArea = chartArea,
                           width = width, height = height,
                           dataOpacity = dataOpacity,
                           backgroundColor = backgroundColor,
                           series = "{0:{color: 'gray', visibleInLegend: true}, 1:{color: 'blue', visibleInLegend: true}}",
                           explorer = "{axis: 'vertical',  keepInBounds: true, maxZoomOut: 1, maxZoomIn: 0.01, actions: ['dragToZoom', 'rightClickToReset']}")
      } else {use.options = options}
      fplot = googleVis::gvisScatterChart(data = dat, options = use.options)
      if(shiny) {
        return(fplot)
      } else {
        graphics::plot(fplot, tag = tag)
      }
    }
  }
  if ("boot" %in% which) {
    var.ident = x$res.df[,vars[1]] == 1
    vi = var.ident
    var.ident[var.ident==TRUE] = paste("With",vars[1])
    var.ident[var.ident==FALSE] = paste("Without",vars[1])
    jitter = stats::runif(length(vi),0-jitterk,0+jitterk)
    dat = data.frame(mods = x$res.df$name,
                     k =  x$res.df$k + jitter,
                     LL = -2*x$res.df$logLikelihood,
                     prob = x$res.df$freq/x$B,
                     var.ident = var.ident)
    if(classic){
      graphics::par(mar = c(3.4,3.4,0.1,0.1),mgp = c(2.0, 0.75, 0))
      if(missing(ylim)) ylim = NULL
      graphics::symbols(dat$k,dat$LL,sqrt(dat$prob),inches=max.circle,
                        bg = ifelse(vi,grDevices::rgb(1, 0, 0, alpha=0.5),grDevices::rgb(0, 0, 1, alpha=0.5)),
                        fg = "white",
                        ylim = ylim,
                        xlab = "Number of parameters",
                        ylab = "-2*Log-likelihood")
      graphics::legend("topright",legend = c(paste("With",vars[1]),paste("Without",vars[1])),
                       col = c(grDevices::rgb(1, 0, 0, alpha=0.5),grDevices::rgb(0, 0, 1, alpha=0.5)),pch=19)
      if(text){
        bdat = dat[dat$prob>min.prob,]
        if(!print.full.model){
          bdat = bdat[-dim(bdat)[1],]
        }
        text(bdat$k,bdat$LL,bdat$mods,cex=0.9,pos=2,offset=0.5,srt=srt)
      }
    } else {
      gvis.title = paste("Model stability plot",sep="")
      x.ticks=paste(1:max(dat$k),collapse=",")
      gvis.hAxis = paste("{title:'Number of parameters',
                         maxValue:",max(dat$k)+0.5," ,
                         minValue:",0.5," ,
                         ticks: [",x.ticks,"]}")
      chartArea = paste("{left:",left,
                        ",top:",top,
                        ",width:'",chartWidth,
                        "',height:'",chartHeight,"'}", sep = "")
      bubble = paste("{opacity:",dataOpacity,
                     ", textStyle: {color: 'none'}}", sep = "")
      if (is.null(options)) {
        use.options = list(title = gvis.title,
                           fontSize = fontSize,
                           vAxis = "{title:'-2*Log-likelihood'}",
                           hAxis = gvis.hAxis,
                           sizeAxis = "{minValue: 0, minSize: 1,
                         maxSize: 20, maxValue:1}",
                           axisTitlesPosition = axisTitlesPosition,
                           bubble = bubble,
                           chartArea = chartArea,
                           width = width, height = height,
                           backgroundColor = backgroundColor,
                           explorer = "{axis: 'vertical',
                         keepInBounds: true,
                         maxZoomOut: 1,
                         maxZoomIn: 0.01,
                         actions: ['dragToZoom',
                         'rightClickToReset']}")
      } else {use.options = options}
      fplot = googleVis::gvisBubbleChart(data = dat, idvar = "mods", xvar = "k",
                                         yvar = "LL", colorvar = "var.ident",
                                         sizevar = "prob",
                                         options = use.options)
      if(shiny){
        return(fplot)
      } else {
        graphics::plot(fplot, tag = tag)
      }
    }
  }
  if ("vip" %in% which) { # variable inclusion plot
    var.names = make.names(x$m$exp.vars)
    var.names = gsub(pattern = ":",replacement = ".",x = var.names)
    B = x$B
    p.var = x$var.in[,is.element(colnames(x$var.in),var.names)]
    colnames(p.var) = gsub("REDUNDANT.VARIABLE","RV",colnames(p.var))
    sortnames = names(sort(apply(p.var,2,mean),decreasing = TRUE))
    vip.df = p.var[,sortnames]
    vip.df = data.frame(lambda = x$lambdas, AIC = NA, AIC.annotation = NA,
                        BIC = NA, BIC.annotation = NA, vip.df/B)
    aicline = rbind(c(2,0,NA,NA,NA,rep(NA,length(var.names))),
                    c(2,1,"AIC",NA,NA,rep(NA,length(var.names))),
                    c(log(x$mextract$n),NA,NA, 0,NA,rep(NA,length(var.names))),
                    c(log(x$mextract$n),NA,NA,1,"BIC",rep(NA,length(var.names))))
    colnames(aicline) = colnames(vip.df)
    vip.df = rbind(vip.df,aicline)
    tid = c(1,2,4,6:dim(vip.df)[2])
    vip.df[, tid] = sapply(vip.df[, tid], as.numeric)
    if (classic) {
      classic.lambda = vip.df$lambda
      classic.vip.df = subset(vip.df,select = -c(get("AIC"),get("AIC.annotation"),
                                                 get("BIC"),get("BIC.annotation"),
                                                 get("lambda")))
      lwds = log((1:length(var.names)) + 1)
      lwds = rev(2*lwds/max(lwds))
      graphics::par(mar = c(3.4,3.4,0.1,0.1), mgp = c(2.0, 0.75, 0))
      graphics::matplot(x = classic.lambda,jitter(as.matrix(classic.vip.df)),type = "l",
                        ylab = "Bootstrapped probability", xlab = "Penalty", lwd = lwds)
      leg.nm = names(classic.vip.df)
      graphics::legend("topright", leg.nm, bg = "transparent", bty = "n", inset = c(0.015),
                       # these are the defaults for matplot:
                       lty = 1:5, col = 1:6, lwd = lwds)
    } else {
      gvis.title = "Variable inclusion plot"
      #lineDashStyle = paste("[",paste(1:2,collapse=","),"]",sep="")
      lineseries = "[{lineDashStyle: [2,2], lineWidth: 2, color:'gray',
      visibleInLegend: false},
      {lineDashStyle: [2,2], lineWidth: 2, color:'gray',
      visibleInLegend: false}]"
      chartArea = paste("{left:",left,
                        ",top:",top,
                        ",width:'",chartWidth,
                        "',height:'",chartHeight,"'}", sep = "")
      if (is.null(options)) {
        use.options = list(title = gvis.title,
                           fontSize = fontSize,
                           vAxis = "{title:'Bootstrapped probability'}",
                           hAxis = "{title:'Penalty'}",
                           sizeAxis = "{minValue: 0, minSize: 1,
                           maxSize: 20, maxValue:1}",
                           axisTitlesPosition = axisTitlesPosition,
                           series = lineseries,
                           #lineDashStyle = lineDashStyle,
                           chartArea = chartArea,
                           width = width, height = height,
                           backgroundColor = 'transparent',
                           annotations = "{style:'line'}")
      } else {use.options = options}
      fplot = googleVis::gvisLineChart(data = vip.df,
                                       xvar = "lambda",
                                       yvar = c("AIC","AIC.annotation",
                                                "BIC","BIC.annotation",
                                                sortnames),
                                       options = use.options)
      if(shiny){
        return(fplot)
      } else {
      return(graphics::plot(fplot, tag = tag))
      }
    }
  } else return(invisible())
  
  
  
  
  
  
  
  
  
}


#' Print method for a vis object
#'
#' Prints basic output of the bootstrap results of an
#' vis object.
#'
#' @param x a \code{vis} object, the result of \code{\link{vis}}
#' @param min.prob a lower bound on the probability of
#'   selection before the result is printed
#' @param print.full.model logical, determines if the full
#'   model gets printed or not.  Default=\code{FALSE}.
#' @param ... further arguments (currently unused)
#' @export
# S3 print method for class 'vis'
print.vis = function (x, min.prob=0.3, print.full.model=FALSE, ...) {
  dat = x$res.df
  dat$prob = dat$freq/x$B
  print.obj = dat[dat$prob>min.prob,c("name","prob","logLikelihood")]
  if(!print.full.model) print.obj = print.obj[-dim(print.obj)[1],]
  print.obj[,2:3] = round(print.obj[,2:3],2)
  rownames(print.obj) = NULL
  print(print.obj,row.names=FALSE)
  invisible(x)
}
