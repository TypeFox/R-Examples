#in development code
#[TBC - NUMBER] functions 

#trianglePlot 
#panel.trianglePlot
#panel.trianglePlotFrame
#panel.trianglePlotGrid
#panel.trianglePlotAxes
#triLimsReset
#triABC2XY
#triXY2ABC


#others now 
#removed



###############################
#new trianglePlot
###############################

#this uses panelPal
#allows conditioning like in standard lattice


trianglePlot <- function(x, data = NULL, ..., ref.cols = TRUE){

#trianglePlot ver 0.2
#panelPal update

    extra.args <- list(...)

    if (is.logical(ref.cols)) 
        ref.cols <- if (all(ref.cols)) 
            c("darkgreen", "darkred", "darkblue")
        else "black"
    ref.cols <- rep(ref.cols, length.out=3)

    extra.args <- listUpdate(list(x=x, data=data, formula.type="z~a0+b0+c0|cond", 
                                  coord.conversion=triABC2XY, panel=panel.trianglePlot, ref.cols=ref.cols), 
                             extra.args)

    do.call(loaPlot, extra.args)

}


##############################
#new panel.trianglePlot
##############################

#uses local.scales methods to handle axes
#local.scales.panel = panel.trianglePlotFrame
#data.panel = panel.loaPlot

#think about getting and adding to default settings in default data.panel?
#this would turn off key of key not set in data.panel?



panel.trianglePlot <- 

function(x = NULL, y = NULL, a0 = NULL, b0 = NULL, c0 = NULL, ..., 
         loa.settings = FALSE, plot = TRUE, process = TRUE){

################
#this is based on panel.polarPlot
################

    if(loa.settings)
        return(list(group.args= c("col"),
                    zcase.args= c("pch"),
                    common.args = c("alim", "blim", "clim", "grid", "axes", "ticks", "annotation", "ref.cols"),
                    default.settings = list(local.scales = TRUE, local.scales.panel = panel.trianglePlotFrame,
                                            data.panel = panel.loaPlot, grid = TRUE, axes = TRUE, 
                                            allowed.scales =c("a0", "b0", "c0"), disallowed.scales = c("x", "y"), 
                                            aspect = "loa.iso", reset.xylims = triLimsReset, 
                                            load.lists = c("grid", "axes", "ticks", "annotation"),                                            
                                            key.fun = "draw.loaPlotZKey")))

    if(process){

###################
#to think about
###################

#make if missing
#x, y, a0, b0, c0,

#see below re polarplot code

           if(!plot) return(list(x = x, y = y, a0=a0, b0=b0, c0=c0))

    }

##
#######################
## the if missing would 
## be like in polarPlot
#######################
##        
##        r <- x
##        theta <- y
##        x <- r * sin(pi * theta/180)
##        y <- r * cos(pi * theta/180) 
##        if(!plot) return(list(x = x, y = y, r = r, theta = theta))
##    }

    if(plot){

        extra.args <- listUpdate(list(...), 
                                 list(x = x, y = y, a0 = a0, b0=b0, c0=c0, 
                                      plot = plot, process = process))

###################
#might want to allow user to set this
#might want to be able to turn this off
#so not local.scales.panels
###################

###################
#might want to tidy
#the do.call 
#data.panel does not 
#need to be sent to
#data.panel...
###################

        if(extra.args$local.scales)
            if(is.function(extra.args$local.scales.panel)){
                do.call(extra.args$local.scales.panel, extra.args)
                extra.args$grid <- NULL
            }
        do.call(extra.args$data.panel, extra.args)
    }
}




##############################
#panel.trianglePlotFrame
##############################

panel.trianglePlotFrame <- 

function (...,  grid = NULL, axes = NULL) 
{
    extra.args <- list(...)

#could tidy this?

    if (isGood4LOA(grid)){
        if(is.list(grid) && is.function(grid$panel))
            do.call(grid$panel, listUpdate(extra.args, list(grid=grid), ignore="panel")) else
            do.call(panel.trianglePlotGrid, listUpdate(extra.args, list(grid=grid), ignore="panel"))
    }
    if (isGood4LOA(axes)){
        if(is.list(axes) && is.function(axes$panel))
            do.call(axes$panel, listUpdate(extra.args, list(axes=axes), ignore="panel")) else
            do.call(panel.trianglePlotAxes, listUpdate(extra.args, list(axes=axes), ignore="panel"))
    }

}


############################
#panel.trianglePlotGrid
############################

panel.trianglePlotGrid <- 

function (alim = NULL, blim = NULL, clim = NULL, 
          ..., grid = TRUE, panel.scales = NULL) 
{
    extra.args <- list(...)

    if (!is.list(panel.scales)) 
        panel.scales <- list()
    if (!is.list(grid)) 
        grid <- list()
    panel.scales <- listUpdate(list(draw = TRUE, arrows = FALSE, 
        lty=3,
        tick.number = 5, abbreviate = FALSE, minlength = 4, tck = 1, 
        col= "lightgrey", cex = 0.8), panel.scales)

    grid <- do.call(listLoad, listUpdate(grid, list(load="a0")))
    grid <- do.call(listLoad, listUpdate(grid, list(load="b0")))
    grid <- do.call(listLoad, listUpdate(grid, list(load="c0")))
    temp <- list(a0=list(col=extra.args$ref.cols[1]),
                 b0=list(col=extra.args$ref.cols[2]),
                 c0=list(col=extra.args$ref.cols[3]))
    if(is.null(grid$col)) grid <- listUpdate(temp, grid)

    grid.pars <- getPlotArgs("axis.line", local.resets = panel.scales, 
                             user.resets = grid, elements = c("a0", "b0", "c0"), 
                             is.scale=TRUE, defaults.only = FALSE)

    at.loc <- function(par, lim){
        temp <- listUpdate(list(tick.number=5), par, use=c("at", "tick.number"))
        temp <- if (!is.null(temp$at)) temp$at else
                                       pretty(lim, temp$tick.number)
#remove any out of range 
#if range too big 
        temp[temp >= min(lim, na.rm=T) & temp <= max(lim, na.rm=T)]
    }

    #a axis
    a.at <- at.loc(grid.pars$a, alim)
    temp <- function(x) {
        a1 <- triABC2XY(c(x, x), c(blim[1], blim[2] - (x - alim[1])), 
            c(clim[2] - (x - alim[1]), clim[1]))
        if (isGood4LOA(grid.pars$a)){
            grid.pars$a$x <- a1$x
            grid.pars$a$y <- a1$y
            do.call(llines, grid.pars$a)
        } 
    }
    for (i in a.at) temp(i)


    #b axis
    b.at <- at.loc(grid.pars$b, blim)
    temp <- function(x) {
        b1 <- triABC2XY(c(alim[1], alim[2] - (x - blim[1])), 
            c(x, x), c(clim[2] - (x - blim[1]), clim[1]))
        if (isGood4LOA(grid.pars$b)){
            grid.pars$b$x <- b1$x
            grid.pars$b$y <- b1$y
            do.call(llines, grid.pars$b)
        } 
    }
    for (i in b.at) temp(i)

    #c axis
    c.at <- at.loc(grid.pars$c, clim)
    temp <- function(x) {
        c1 <- triABC2XY(c(alim[1], alim[2] - (x - clim[1])), 
            c(blim[2] - (x - clim[1]), blim[1]), c(x, x))
        if (isGood4LOA(grid.pars$c)){
            grid.pars$c$x <- c1$x
            grid.pars$c$y <- c1$y
            do.call(llines, grid.pars$c)
        } 
    }
    for (i in c.at) temp(i)

}




############################
#panel.trianglePlotAxes
############################

panel.trianglePlotAxes <- 

function (alim = NULL, blim = NULL, clim = NULL, ..., 
          axes = TRUE, ticks=TRUE, annotation=TRUE, 
          panel.scales = NULL) 
{

#reposition a and b labs so parallel to axes
#look and management of font size in local labels handler
#

    extra.args <- list(...)

    alab <- if(is.null(extra.args$alab)) extra.args$a0lab else extra.args$alab
    blab <- if(is.null(extra.args$blab)) extra.args$b0lab else extra.args$blab
    clab <- if(is.null(extra.args$clab)) extra.args$c0lab else extra.args$clab

    if (!is.list(panel.scales)) 
        panel.scales <- list()

    if (!is.list(axes)) 
        axes <- list()
#    if (!is.list(ticks)) 
#        ticks <- list()
#    if (!is.list(annotation)) 
#        annotation <- list()

    temp <- prod(dim(trellis.currentLayout()),na.rm=T)
    text.cex <- 1
    if(temp>1) text.cex <- 0.8 
    if(temp>3) text.cex <- 0.7
    if(temp>3) text.cex <- 0.6
    if(temp>9) text.cex <- 0.5 

    temp <- list(a0=list(col=extra.args$ref.cols[1]),
                 b0=list(col=extra.args$ref.cols[2]),
                 c0=list(col=extra.args$ref.cols[3]))



    axes <- do.call(listLoad, listUpdate(axes, list(load="a0")))
    axes <- do.call(listLoad, listUpdate(axes, list(load="b0")))
    axes <- do.call(listLoad, listUpdate(axes, list(load="c0")))

    if(isGood4LOA(ticks)){
       if(!is.list(ticks)) ticks <- list()
       ticks <- do.call(listLoad, listUpdate(ticks, list(load="a0")))
       ticks <- do.call(listLoad, listUpdate(ticks, list(load="b0")))
       ticks <- do.call(listLoad, listUpdate(ticks, list(load="c0")))
       if(is.null(ticks$col)) ticks <- listUpdate(temp, ticks)
    } else ticks <- list(col=NA)

    if(isGood4LOA(annotation)){
       if(!is.list(annotation)) annotation <- list()
       annotation <- do.call(listLoad, listUpdate(annotation, list(load="a0")))
       annotation <- do.call(listLoad, listUpdate(annotation, list(load="b0")))
       annotation <- do.call(listLoad, listUpdate(annotation, list(load="c0")))
       if(is.null(annotation$cex)) annotation$cex <- (text.cex*0.8)
       if(is.null(annotation$col)) annotation <- listUpdate(temp, annotation)
       
    } else annotation <- list(col=NA)

#    panel.scales <- listUpdate(list(draw = TRUE, arrows = FALSE, 
#        tick.number = 5, abbreviate = FALSE, minlength = 4, tck = 1, 
#        col = "red", col.line = 1, cex = 0.8), panel.scales)


    axis.loc <- function(n, lim) (n * (max(lim, na.rm = TRUE) - 
        min(lim, na.rm = TRUE))) + min(lim, na.rm = TRUE)
    at.loc <- function(par, axes, ticks, lim){
        temp <- listUpdate(par, axes, use=c("at", "tick.number"))
        temp <- listUpdate(temp, ticks, use=c("at", "tick.number"))
        if(!is.null(temp$at)) temp$at else pretty(lim, temp$tick.number)
    }
    axes.pars <- getPlotArgs(default.as = "axis.line", source = panel.scales, 
        elements = c("a0", "b0", "c0"), is.scales = TRUE, user.resets = axes)
    tick.pars <- getPlotArgs(default.as = "axis.line", source = panel.scales, 
        elements = c("a0", "b0", "c0"), is.scales = TRUE, user.resets = ticks)

    ann.pars <- getPlotArgs(default.as = "axis.text", source = panel.scales, 
        elements = c("a0", "b0", "c0"), is.scales = TRUE, user.resets = annotation)

#this fixes the current issue with getPlotArgs
    ann.pars$a0 <- listUpdate(annotation, ann.pars$a0, ignore.a=c("a0", "b0", "c0"))
    ann.pars$b0 <- listUpdate(annotation, ann.pars$b0, ignore.a=c("a0", "b0", "c0"))
    ann.pars$c0 <- listUpdate(annotation, ann.pars$c0, ignore.a=c("a0", "b0", "c0"))

#need to fix this

###currently changes to ticks but not others
###this gives us full axes control but 
###using n0 axis name not n 
###print(tick.pars)

    panel.scales$a <- listUpdate(list(at = at.loc(list(tick.number=5), axes.pars$a, tick.pars$a, alim)),
                                 panel.scales$a)
 #   tick.pars$a$isGood4LOA <- TRUE
 #   ann.pars$a$isGood4LOA <- TRUE
    temp <- triABC2XY(c(axis.loc(0.5, alim), axis.loc(0.5, alim), 
        axis.loc(0, alim), axis.loc(1, alim)), c(axis.loc(0.1, 
        blim), axis.loc(0, blim), axis.loc(0, blim), axis.loc(0, 
        blim)), c(axis.loc(0.4, clim), axis.loc(0.5, clim), axis.loc(1, 
        clim), axis.loc(0, clim)))
    panel.localScale(panel.scale = panel.scales$a, x.loc = temp$x[3:4], 
        y.loc = temp$y[3:4], lim = alim, x.offset = temp$x[2] - 
            temp$x[1], y.offset = temp$y[2] - temp$y[1], axis = axes.pars$a0, 
            ticks=tick.pars$a0, annotation=ann.pars$a0)
#    ltext(x = temp$x[1] - (3 * (temp$x[2] - temp$x[1])) - (3 * 
#        (temp$x[1] - temp$x[3])), y = temp$y[1] + (2 * (y.offset = temp$y[2] - 
#        temp$y[1])), alab, adj = c(1, 0.5), srt=60)

    ltext(x = temp$x[1] + (4 * (temp$x[2] - temp$x[1])), y = temp$y[1] + 
        (3 * (y.offset = temp$y[2] - temp$y[1])), alab, adj = c(0.5, 
        0.5), srt=60, cex = text.cex)



    panel.scales$b <- listUpdate(list(at = at.loc(list(tick.number=5), axes.pars$b, tick.pars$b, blim)),
                                 panel.scales$b)

#    tick.pars$b$isGood4LOA <- TRUE
#    ann.pars$b$isGood4LOA <- TRUE

    temp <- triABC2XY(c(axis.loc(0.4, alim), axis.loc(0.5, alim), 
        axis.loc(1, alim), axis.loc(0, alim)), c(axis.loc(0.5, 
        blim), axis.loc(0.5, blim), axis.loc(0, blim), axis.loc(1, 
        blim)), c(axis.loc(0.1, clim), axis.loc(0, clim), axis.loc(0, 
        clim), axis.loc(0, clim)))
    panel.localScale(panel.scale = panel.scales$b, x.loc = temp$x[3:4], 
        y.loc = temp$y[3:4], lim = blim, x.offset = temp$x[2] - 
            temp$x[1], y.offset = temp$y[2] - temp$y[1], axis = axes.pars$b0, 
        ticks = tick.pars$b0, annotation = ann.pars$b0)

    ltext(x = temp$x[1] + (4 * (temp$x[2] - temp$x[1])), y = temp$y[1] + 
        (3 * (y.offset = temp$y[2] - temp$y[1])), blab, adj = c(0.5, 
        0.5), srt=300, cex = text.cex)


    panel.scales$c <- listUpdate(list(at = at.loc(list(tick.number=5), axes.pars$c, tick.pars$c, clim)),
                                 panel.scales$c)

 #   tick.pars$c$isGood4LOA <- TRUE
 #   ann.pars$c$isGood4LOA <- TRUE

    temp <- triABC2XY(c(axis.loc(0.1, alim), axis.loc(0, alim), 
        axis.loc(0, alim), axis.loc(0, alim)), c(axis.loc(0.4, 
        blim), axis.loc(0.5, blim), axis.loc(1, blim), axis.loc(0, 
        blim)), c(axis.loc(0.5, clim), axis.loc(0.5, clim), axis.loc(0, 
        clim), axis.loc(1, clim)))
    panel.localScale(panel.scale = panel.scales$c, x.loc = temp$x[3:4], 
        y.loc = temp$y[3:4], lim = clim, x.offset = temp$x[2] - 
            temp$x[1], y.offset = temp$y[2] - temp$y[1], axis = axes.pars$c0, 
        ticks = tick.pars$c0, annotation = ann.pars$c0)

    ltext(x = temp$x[2], y = temp$y[1] + (3 * (y.offset = temp$y[2] - 
        temp$y[1])), clab, adj = c(0.5, 0.5), cex=text.cex)




}




#############################
#############################
##data handlers
#############################
#############################




###############################
#triLimsReset
###############################

triLimsReset <- 

function(ans){

    #what is not a, b, clim values
    #does function in preprocess make them?

#messy 

    temp <- ans$panel.args.common
    temp <- triABC2XY(a = c(temp$alim[1], temp$alim[1], temp$alim[2]), 
        b = c(temp$blim[1], temp$blim[2], temp$blim[1]), c = c(temp$clim[2], 
            temp$clim[1], temp$clim[1]), verbose = FALSE)
    xlim <- range(temp$x, na.rm = TRUE)
    ylim <- range(temp$y, na.rm = TRUE)

    ans$panel.args.common$xlim <- xlim
        ans$x.limits <- xlim
        ans$panel.args.common$ylim <- ylim
        ans$y.limits <- ylim

    temp <- function(lim, q1, q2) {
        if (diff(lim) == 0) 
            lim + q1
        else lim + c(-(diff(lim)/q2[1]), (diff(lim)/q2[2]))
    }

#messy

   extra.args <- ans$panel.args.common

    alab <- if(is.null(extra.args$alab)) extra.args$a0lab else extra.args$alab
    blab <- if(is.null(extra.args$blab)) extra.args$b0lab else extra.args$blab
    clab <- if(is.null(extra.args$clab)) extra.args$c0lab else extra.args$clab


   if (is.null(clab)) {
        extra.args$xlim <- temp(extra.args$xlim, c(-0.5, 0.5), 
            c(5, 5))
        extra.args$ylim <- temp(extra.args$ylim, c(-0.5, 0.5), 
            c(5, 5))
    }
    else if (is.character(clab) && clab == "") {
        extra.args$xlim <- temp(extra.args$xlim, c(-0.5, 0.5), 
            c(5, 5))
        extra.args$ylim <- temp(extra.args$ylim, c(-0.5, 0.5), 
            c(5, 5))
    }
    else {
        extra.args$xlim <- temp(extra.args$xlim, c(-0.5, 0.5), 
            c(5, 5))
        extra.args$ylim <- temp(extra.args$ylim, c(-0.3, 0.5), 
            c(3, 5))
    }

    ans$panel.args.common <- extra.args

  
        ans$x.limits <- ans$panel.args.common$xlim
        ans$y.limits <- ans$panel.args.common$ylim 

    ans

}




##############################
##############################
##triABC2XY
##############################
##############################

triABC2XY <- function(a, b=NULL, c=NULL, ..., force.abc=TRUE, 
              if.na="remove.row", if.neg="remove.row", verbose=FALSE){

    #############
    #setup
    #############

    #extra.args
    extra.args <- list(...)

#############
#new
#############

#if a,b,c not there 
#use a0, b0, and c0 if there 

    if(missing(a) && "a0" %in% names(extra.args)) a <- extra.args$a0
    if(is.null(b) && "b0" %in% names(extra.args)) b <- extra.args$b0
    if(is.null(c) && "c0" %in% names(extra.args)) c <- extra.args$c0

#############
#

    #make a,b,c a data.frame

################
#could standardise this next bit and front end of triXY2ABC
#and make common function dataHandler
#could also put the logs in there
################

    data.abc <- if(is.data.frame(a)) a else
                    if(is.list(a)) as.data.frame(a) else
                       if(is.vector(a)) data.frame(a=a) else
                           stop("unable to handle supplied data", call. = FALSE)
###################
#possible issue if
#a shorter than b,c
###################
    if(is.vector(b))
        data.abc$b <- if(length(b) < nrow(data.abc))
                          rep(b, ceiling(nrow(data.abc)/length(b)))[1:nrow(data.abc)] else 
                              b[1:nrow(data.abc)]
    if(is.vector(c))
        data.abc$c <- if(length(c) < nrow(data.abc))
                          rep(c, ceiling(nrow(data.abc)/length(c)))[1:nrow(data.abc)] else 
                              c[1:nrow(data.abc)]
    #check dim
    if(ncol(data.abc) < 3)
        stop("insufficient data for 'abc' assignment", call. = FALSE)

    #force.abc/rescale
    temp <- data.abc[,1:3]
    if(force.abc){
        if("a" %in% names(data.abc)) temp[,1] <- data.abc$a
        if("b" %in% names(data.abc)) temp[,2] <- data.abc$b
        if("c" %in% names(data.abc)) temp[,3] <- data.abc$c
    } 
    data.abc <- temp
    abc.status <- rep(0, nrow(data.abc)) #abc.status

    ###########
    #if.neg and if.na
    ###########

##############
#need keep.as.is catcher
##############
    na.log <- apply(data.abc, 1, function(x) any(is.na(x))) #na
    neg.log <- apply(data.abc, 1, function(x) any(!is.na(x) & x<0)) #negs

    if(any(na.log)) {
        if(if.na == "remove.row")
            data.abc[na.log, 1:3] <- c(NA,NA,NA)
        if(if.na == "make.zero")
            data.abc[is.na(data.abc)] <- 0
    }

    if(any(neg.log)) {
        if(if.neg == "remove.row")
            data.abc[neg.log, 1:3] <- c(NA,NA,NA)
        if(if.neg == "make.zero")
            data.abc[data.abc<0] <- 0
        if(if.neg == "rescale.col")
            if(nrow(data.abc)==1)
                data.abc[data.abc<0] <- 0 else 
                    data.abc <- as.data.frame(apply(data.abc, 2, function(x) 
                        if(min(x, na.rm=TRUE)<0) x-min(x, na.rm=TRUE) else x))            
    }

#############
#below needs documenting in help
#############

    #abc 2 prop(abc)
    #function used again later
    prop.abc <- function(d){
        temp <- d[,1] + d[,2] + d[,3]
        d/temp
    }
    data.abc <- prop.abc(data.abc)

################################
################################

##############
#lim forcings
##############

    ##################
    #following based on ade4/Daniel Chessel method
    #################
    ##selection of sensible lim range

#could move this to separate/as data handlers
#

    #lim forcings
    if(is.null(extra.args$lims)){
        lims <- apply(data.abc, 2, range, na.rm = TRUE, finite = TRUE)
        lims[1,] <- (floor(lims[1,]/0.1))/10
        lims[2,] <- (floor(lims[2,]/0.1) + 1)/10
    } else {
              lims <- as.data.frame(matrix(rep(extra.args$lims[1:2], 3), 
                                    ncol = 3, nrow = 2))
    }
    if(!is.null(extra.args$alim))
        lims[,1] <- extra.args$alim
    if(!is.null(extra.args$blim))
        lims[,2] <- extra.args$blim
    if(!is.null(extra.args$clim))
        lims[,3] <- extra.args$clim
    if(!is.null(extra.args$abc.mins))
        lims[1,] <- extra.args$abc.mins
    if(!is.null(extra.args$abc.maxs))
        lims[2,] <- extra.args$abc.maxs

    #stop the bad
    lims[1,][lims[1,] < 0] <- 0
    lims[2,][lims[2,] > 1] <- 1

    #lims range
    #ade4 method

##############
#rethink this?
##############

    lim.lite <-function(lims){
        temp <- lims[2,] - lims[1,]
        temp2 <- max(temp)
        if (!all(temp == temp2)) {
            for (j in 1:3) {
                k <- temp2 - temp[j]
                while (k > 0) {
                    if ((k > 0) & (lims[2,j] < 1)) {
                        lims[2,j] <- lims[2,j] + 0.1
                        k <- k - 1
                    }
                    if ((k > 0) & (lims[1,j] > 0)) {
                        lims[1,j] <- lims[1,j] - 0.1
                        k <- k - 1
                    }
                }
            }
        }
        cal <- matrix(0, 9, 3)
        cal[1, 1] <- lims[1,1]
        cal[1, 2] <- lims[1,2]
        cal[1, 3] <- 1 - cal[1, 1] - cal[1, 2]
        cal[2, 1] <- lims[1,1]
        cal[2, 2] <- lims[2,2]
        cal[2, 3] <- 1 - cal[2, 1] - cal[2, 2]
        cal[3, 1] <- lims[2,1]
        cal[3, 2] <- lims[1,2]
        cal[3, 3] <- 1 - cal[3, 1] - cal[3, 2]
        cal[4, 1] <- lims[1,1]
        cal[4, 3] <- lims[1,3]
        cal[4, 2] <- 1 - cal[4, 1] - cal[4, 3]
        cal[5, 1] <- lims[1,1]
        cal[5, 3] <- lims[2,3]
        cal[5, 2] <- 1 - cal[5, 1] - cal[5, 3]
        cal[6, 1] <- lims[2,1]
        cal[6, 3] <- lims[1,3]
        cal[6, 2] <- 1 - cal[6, 1] - cal[6, 3]
        cal[7, 2] <- lims[1,2]
        cal[7, 3] <- lims[1,3]
        cal[7, 1] <- 1 - cal[7, 2] - cal[7, 3]
        cal[8, 2] <- lims[1,2]
        cal[8, 3] <- lims[2,3]
        cal[8, 1] <- 1 - cal[8, 2] - cal[8, 3]
        cal[9, 2] <- lims[2,2]
        cal[9, 3] <- lims[1,3]
        cal[9, 1] <- 1 - cal[9, 2] - cal[9, 3]
        lims[1,] <- apply(cal, 2, min)
        lims[1,] <- round(lims[1,], digits = 4)
        lims[2,] <- apply(cal, 2, max)
        lims[2,] <- round(lims[2,], digits = 4)

###################
#new bit
#bad lims

        temp <- lims[2,] -lims[1,]
        if(length(unique(temp))>1)
           lims[2,] <- lims[1,] + max(temp, na.rm=TRUE)

#####################

        #stop the bad
        lims[1,][lims[1,] < 0] <- 0
        lims[2,][lims[2,] > 1] <- 1
        lims
    }
    lims <- lim.lite(lims)

#check if next bit needed

    temp <- lims[2,] - lims[1,]
    temp2 <- max(temp)
    if (!all(temp == temp2))
          lims <- lim.lite(lims)
    if (!all(temp == temp2))
          lims <- lim.lite(lims)



################################
################################

    #check for out of range values
    ##compare with lims
    oor.log <- rep(FALSE, nrow(data.abc))
    oor.log <- ifelse(data.abc[,1] < min(lims[,1], na.rm=TRUE) |  
                      data.abc[,1] > max(lims[,1], na.rm=TRUE), TRUE, oor.log)
    oor.log <- ifelse(data.abc[,2] < min(lims[,2], na.rm=TRUE) |  
                      data.abc[,2] > max(lims[,2], na.rm=TRUE), TRUE, oor.log)
    oor.log <- ifelse(data.abc[,3] < min(lims[,3], na.rm=TRUE) |  
                      data.abc[,3] > max(lims[,3], na.rm=TRUE), TRUE, oor.log)
#catch na's
    oor.log[is.na(oor.log)]<- FALSE

    data.abc[oor.log, 1:3] <- c(NA,NA,NA)

    #################
    ##following based on Mu Zhu, 
    ##Statistical Computing and Graphics 19,1, 2008
    #################
    #eta <- 0.08
    #anchor <- diag(3)
    #alpha1 <- c(1, -1, 0)/sqrt(2)
    #alpha2 <- c(-0.5, -0.5, 1)/sqrt(1.5)
    #vert.x <- anchor %*% alpha1
    #vert.y <- anchor %*% alpha2
 
    #x <- as.matrix(data.abc) * alpha1
    #y <- as.matrix(data.abc) * alpha2

    #ans <- list(x = x, 
    #            y = y,
    #            alim = lims[,1], blim = lims[,2], clim = lims[,3])

    ##################
    #following based on ade4/Daniel Chessel method
    #################
    ##with this abc -> xy scaling is local so triangle size remains constant
    ##if this is adopted we could simplify panel... functions
    ##rescale axis
    ##rescale.axis <- function(n, lim)
    ##                    (n - lim[1])/(lim[2]-lim[1])
    ##temp <- data.abc
    ##temp[,1] <- rescale.axis(temp[,1], lims[,1])
    ##temp[,2] <- rescale.axis(temp[,2], lims[,2])
    ##temp[,3] <- rescale.axis(temp[,3], lims[,3])
    ##temp <- prop.abc(temp)
    #basic return x,y a/b/clim
    ##ans <- list(x = (temp[,1] - temp[,3])/sqrt(2),
    ##            y = (2 * temp[,2] - temp[,1] - temp[,3])/sqrt(6),
    ##            alim = lims[,1], blim = lims[,2], clim = lims[,3])

    ##alternative method from Leic method?
    ##with this abc -> xy scaling fixed and triangle size changes
    ##with this panel... stay as is but triXY2ABC is greatly simplified
    ##also lim setting much less flexible
    
    ans <- list(x =  data.abc[,2]+(0.5*data.abc[,1]), 
                y = ((data.abc[,1]*0.866)*1.1)/1,      #confirm 1.1
                alim = lims[,1], blim = lims[,2], clim = lims[,3])

    if(!verbose) return(ans) 

    #full return
    #may want to rethink structure
    #re passing a,b,c to xyplot.formula...

    listUpdate(ans, list(a=data.abc[,1], b=data.abc[,2], c=data.abc[,3],
                         report = list(nas=na.log, negs=neg.log, oor=oor.log)))    
}



##############################
##############################
##triXY2ABC
##############################
##############################

triXY2ABC <- function(x, y=NULL, ..., force.xy=TRUE, verbose=FALSE){

    #############
    #setup
    #############
    #make xy a data.frame
    data.xy <- if(is.data.frame(x)) a else
                    if(is.list(x)) as.data.frame(x) else
                       if(is.vector(x)) data.frame(x=x) else
                           stop("unable to handle supplied data", call. = FALSE)
###################
#possible issue if
#x shorter than y - as above
###################

#need a lim checker/out of range handler
#need 

    if(is.vector(y))
        data.xy$y <- if(length(y) < nrow(data.xy))
                          rep(y, ceiling(nrow(data.xy)/length(y)))[1:nrow(data.xy)] else 
                              y[1:nrow(data.xy)]

    #check dim
    if(ncol(data.xy) < 2)
        stop("insufficient data for 'xy' assignment", call. = FALSE)

    #force.abc/rescale
    temp <- data.xy[,1:2]
    if(force.xy){
        if("x" %in% names(data.xy)) temp[,1] <- data.xy$x
        if("y" %in% names(data.xy)) temp[,2] <- data.xy$y
    } 
    data.xy <- temp

    a <- data.xy[,2]/(1.1 *0.866) * 1   #comfirm 1.1
    b <- (data.xy[,1] - (a * 0.5))
    c <- 1 - (a+b)

    ans <- list(a=a, b=b, c=c)

    if(!verbose) return(ans)    

    listUpdate(ans, list(x=data.xy[,1], y=data.xy[,2]))    
}


