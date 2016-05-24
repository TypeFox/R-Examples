pyramidlattice <-
function (x, data = parent.frame(), panel = panel, prepanel = prepanel, 
    strip = TRUE, box.ratio = 2, groups = NULL, beside = FALSE, 
    horizontal = NULL, subset = TRUE, subscripts = !is.null(groups), 
    ...) 
{
    dots <- list(...)
    groups <- eval(substitute(groups), data, parent.frame())
    subset <- eval(substitute(subset), data, parent.frame())
    if (!is.function(panel)) 
        panel <- eval(panel)
    if (!is.function(strip)) 
        strip <- eval(strip)
    prepanel <- if (is.function(prepanel)) 
        prepanel
    else if (is.character(prepanel)) 
        get(prepanel)
    else eval(prepanel)
    x <- do.call("pyramid2", c(list(x = x, data = data, horizontal = horizontal, 
          beside = beside, groups = groups, subscripts = subscripts, 
          subset = subset, panel = panel, prepanel = prepanel, 
          strip = strip, box.ratio = box.ratio), dots))
          
    limax <- max(abs(x$x.limits)) + max(abs(x$x.limits)) * 0.05
    labmax <- limax - limax / 5
    
    hilo <- signif(seq(labmax,0,length.out=4),2)
    lohi <- signif(seq(0,labmax,length.out=4),2)
    
    xlimits <- rep(list(c(-limax,0),c(0,limax)),times=prod(x$layout)/2)
    xlabels <- rep(list(hilo,lohi),times=prod(x$layout)/2)
    xat <- rep(list(-hilo,lohi),times=prod(x$layout)/2)
    
    update(x,scales=list(x=list(relation="free",limits=xlimits,labels=xlabels,at=xat)))
}
