#in development code
#[TBC - NUMBER] functions 

#stackPlot
#stackXZ2XYZ 
#panel.stackPlot


###############################
#stackPlot
###############################

#this uses panelPal
#allows conditioning like in standard lattice


stackPlot <- function (x, data = NULL, ...){

    extra.args <- list(...)
    extra.args <- listUpdate(list(x = x, data = data, formula.type = "z~x|cond", 
        coord.conversion = stackXZ2XYZ, panel = panel.stackPlot), 
        extra.args)

    do.call(loaPlot, extra.args)
}



###############################
#stackXZ2XYZ
###############################

stackXZ2XYZ <- function(col=NULL, col.regions=NULL, par.settings=NULL, scheme=NULL, border=NULL, key.handling = NULL, force.key=NULL,...){

     #might rethink this
     #might not export it

     #set up
     extra.args <- list(...)

     #y data and labels
     extra.args$y <- extra.args$z
     if(!"ylab" %in% names(extra.args))
         extra.args$ylab <- extra.args$zlab

     #defualt one col
     cols <- 1

     if(!"zcases" %in% names(extra.args)) {
         zcases <- rep("default", length(extra.args$y))
         zcase.ids <- "default"
     } else {
         zcases <- extra.args$zcases
         zcase.ids <- extra.args$zcase.ids
         cols <- 1:length(zcase.ids)
     }

     if("panel.condition" %in% names(extra.args)){
         zcases <- extra.args$panel.condition$zcases
         zcase.ids <- unique(extra.args$panel.condition$zcases)
     }

     ref.x <- extra.args$x[zcases == zcase.ids[1]]
     ref.y <- rep(0, length(ref.x))
     extra.args$y0 <- extra.args$y

     for(i in zcase.ids){
         temp <- extra.args$y[zcases == i]
         temp <- temp - min(temp, na.rm=T) #range[1] might be more robust?
         temp[is.na(temp)] <- 0
         extra.args$y[zcases == i] <- temp + ref.y
         extra.args$y0[zcases == i] <- ref.y
         if(!"panel.condition" %in% names(extra.args)) ref.y <- ref.y + temp
     }


#this is the cheat to recolour the plot
#this works because it is only applied to 
#in the colHandler call below

     if(is.null(scheme))
         scheme <- "kr.web"

#need to look into why include par.settings = NULL 
#stop colHandler looking at scheme

     extra.args$col <- colHandler(z=cols, ref=cols, col=col, col.regions=col.regions, scheme=scheme)
     if(is.null(border))
         extra.args$border <- extra.args$col

#this is the cheat for the new key handling 

     if(is.null(force.key))
         extra.args$ycase.key.method2 <- TRUE

     extra.args
}




###############################
#panel.stackPlot
###############################

panel.stackPlot <- function (..., process=TRUE, plot=TRUE, loa.settings = FALSE) 
{
    if (loa.settings) 
        return(list(group.args = c("pch"), zcase.args = c("col", "border"), 
            common.args = c("col", "border"), 
            default.settings = list(key.fun = "draw.ycasePlotKey",
                grid = TRUE, order.x=TRUE)))
    extra.args <- list(...)

#can shift process but...

#    if(process){
#         #extra.args$y <- extra.args$y - min(extra.args$y, na.rm=TRUE)
#    }

#need to return something it can handle

    if(!plot){
        return(list(y = extra.args$y))
    }
    if(plot){

    #grid layer
    if (isGood4LOA(extra.args$grid)) 
        panel.loaGrid(panel.scales = extra.args$panel.scales, 
            grid = extra.args$grid)

    #layers of stack
    plot.fun <- function(...) {
        extra.args <- list(...)

###################
#new bit
#to order by x to make panels neater
#
        #z and y0 might not exist 
        #if this does not come from stack plot
        if("order.x" %in% names(extra.args) && is.logical(extra.args$order.x) && extra.args$order.x){
            df <- as.data.frame(extra.args[names(extra.args) %in% c("x", "y", "z", "y0")])
            df <- df[order(df$x),] 
            extra.args <- listUpdate(extra.args, as.list(df))
        }
#ends here
###################

        extra.args$x <- c(extra.args$x, rev(extra.args$x))
        if(!"y0" %in% names(extra.args))
            extra.args$y0 <- rep(0, length(extra.args$x))
        extra.args$y <- c(extra.args$y, rev(extra.args$y0)) 
        do.call(panel.polygon, extra.args)
    }
    do.call(groupsAndZcasesPanelHandler, listUpdate(list(...), 
        list(panel = plot.fun)))

    }
}
















