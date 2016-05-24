wtd.boxplot.default <-
function(x, weights=NULL, ..., range = 1.5, width = NULL, varwidth = FALSE,
         notch = FALSE, outline = TRUE, names, plot = TRUE,
         border = par("fg"), col = NULL, log = "",
         pars = list(boxwex = 0.8, staplewex = 0.5, outwex = 0.5),
         horizontal = FALSE, add = FALSE, at = NULL)
{
args <- list(x, ...)
    namedargs <-
if(!is.null(attributes(args)$names))
    attributes(args)$names != ""
else
    rep(FALSE, length.out = length(args))
    pars <- c(args[namedargs], pars)
    groups <- if(is.list(x)) x else args[!namedargs]
    if (!is.null(weights)){
    if(!is.list(weights)) weights<-list(weights)
    datasize<-sapply(groups,length)
    wtsize<-sapply(weights,length)
    if (length(datasize)!=length(wtsize))
      stop("number of groups for data and weights are different")
    if (any(datasize != wtsize))
        stop("group sizes for data and weights are different")
    groupwts<-weights
    }
    else groupwts<-NULL
    if(0 == (n <- length(groups)))
stop("invalid first argument")
    if(length(class(groups)))
groups <- unclass(groups)
    if(!missing(names))
attr(groups, "names") <- names
    else {
if(is.null(attr(groups, "names")))
    attr(groups, "names") <- 1:n
        names <- attr(groups, "names")
    }
    for(i in 1:n) {
  if(is.null(groupwts[[i]]))
groups[i] <- list(wtd.boxplot.stats(groups[[i]],
   weights=NULL,
   coef=range)) # do.conf=notch)
   else
groups[i] <- list(wtd.boxplot.stats(groups[[i]],
  weights=groupwts[[i]],
   coef=range)) # do.conf=notch)
   }
    stats <- matrix(0,nrow=5,ncol=n)
    conf  <- matrix(0,nrow=2,ncol=n)
    ng <- out <- group <- numeric(0)
    ct <- 1
    for(i in groups) {
stats[,ct] <- i$stats
        conf [,ct] <- i$conf
        ng <- c(ng, i$n)
        if((lo <- length(i$out))) {
            out   <- c(out,i$out)
            group <- c(group, rep.int(ct, lo))
        }
        ct <- ct+1
    }
    z <- list(stats = stats, n = ng, conf = conf, out = out, group = group,
              names = names)
    if(plot) {
bxp(z, width, varwidth = varwidth, notch = notch, log = log,
#            border = border, col = col, pars = pars,
            border = border, boxfill = col, pars = pars,
            outline = outline, horizontal = horizontal, add = add, at = at)
invisible(z)
    }
    else z
}

