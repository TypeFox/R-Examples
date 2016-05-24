#-->> BEGIN general summary and plotting routines for the NADA package

## Generics

setGeneric("censummary", 
           function(obs, censored, groups) standardGeneric("censummary"))

## Classes

setClass("censummary", "list")

## Methods

setMethod("censummary", 
signature(obs="numeric", censored="logical", groups="missing"),
function (obs, censored, groups) 
{
    ret = cohn(obs, censored)
    ret$n = length(obs)
    ret$n.cen = length(obs[censored])

    props = c(ret$n, ret$n.cen, pctCen(obs, censored), min(obs), max(obs))
    names(props) = c("n", "n.cen", "pct.cen", "min", "max")

    limits = data.frame(ret$limit, ret$C, ret$A, ret$P)
    colnames(limits) = c("limit", "n", "uncen", "pexceed")

    return(new("censummary", list(all=props, limits=limits)))
})

setMethod("censummary", 
signature(obs="numeric", censored="logical", groups="factor"),
function (obs, censored, groups) 
{
    ret = list()
    for (i in levels(groups))
      {
        x  = obs[groups == i]
        xc = censored[groups == i]
        ret[[i]] = censummary(x, xc)
        names(ret[i]) = i
      }
    new("NADAList", ret)
})

setMethod("censummary", 
signature(obs="numeric", censored="logical", groups="numeric"),
function (obs, censored, groups) 
{
    censummary(obs, censored, as.factor(groups))
})

setMethod("show", signature("censummary"), function(object)
{
    for (i in names(object))
      {
        cat(i, ":\n", sep="")
        print(object[[i]])
        cat("\n")
      }
})


censtats =
function(obs, censored) 
{
    skm  = cenfit(obs, censored)
    sros = cenros(obs, censored)
    smle = cenmle(obs, censored)

    med  = c(median(skm), median(sros), median(smle))
    sd   = c(sd(skm), sd(sros), sd(smle)[1])
    mean = c(mean(skm)[1], mean(sros), mean(smle)[1])

    len = c(length(obs), length(which(censored == T)), pctCen(obs, censored))
    names(len) = c("n", "n.cen", "pct.cen")
    print(len)

    data.frame(median=med, mean=mean, sd=sd, row.names=c("K-M", "ROS", "MLE"))
}

## BEGIN general plots routines

# Dennis' censored boxplots 
cenboxplot =
function(obs, cen, group, log=TRUE, range=0, ...) 
{
  if (log) log="y"
  else     log=""

  if (missing(group)) 
      ret = boxplot(cenros(obs, cen), log=log, range=range, ...)
  else
    {
      modeled = numeric()
      groups  = character()
      for (i in levels(as.factor(group)))
        {
            mod = suppressWarnings(
                    cenros(obs[group == i], cen[group == i])$modeled)
            grp = rep(i, length(mod))

            modeled = c(modeled, mod)
            groups = c(groups, grp)
        }
      boxplot(modeled~as.factor(groups), log=log, range=range, ...)
      ret = data.frame(ros.model=modeled, group=groups)
    }

  # Draw horiz line at max censored value
  abline(h=max(obs[cen])) 

  invisible(ret)
}

# Dennis' censored xy plots -- need to fix this
cenxyplot =
function (x, xcen, y, ycen, log="", lty="dashed", ...) 
{
    plot(x, y, log = log, type = "n", ...)

    points(x[!ycen], y[!ycen], ...)

    x0 = x1 = x[ycen]
    y0 = rep(cenpar.usr(log)[3], length(x[ycen]))
    y1 = y[ycen]

    segments(x0, y0, x1, y1, lty=lty, ...)
}

## END general plots routines


#-->> END general summary and plotting routines for the NADA package
