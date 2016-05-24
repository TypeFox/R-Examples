#in development code
#[TBC - NUMBER] functions 

#panel.zcasePieSegmentPlot
#panel.zcasePiePlot

#NOTE: much borrowed from... 

#to do

############################
#repairs
############################

#



###############################
###############################
##panel.zcasePiePlot
###############################
###############################

panel.zcasePiePlot <- function (..., zcase.rescale = TRUE, 
                                       loa.settings = FALSE){

    ################################
    #panel.zcasePieSegmentPlot v0.2 
    #kr 
    ################################
    #draws the pieSegment plot
    ################################
    #notes
    #
    #faster revision of v0.1
    #revision as part of more flexible 
    #general panel handler for these plots
    ################################

    if (loa.settings) 
        return(list(zcase.args = c("col"),  
            default.settings = list(key.fun = "draw.zcasePlotKey", 
                grid = FALSE)))

    extra.args <- list(...)

#temp col fix
#needs work so borders and segment fill 
#transparency sepearated

    if("alpha" %in% names(extra.args)){
        extra.args$alpha.regions <- extra.args$alpha
        extra.args$alpha <- NULL
    }

#compare these with getZcaseDimensions
#could do this as a make defaults option?

    if(!"z" %in% names(extra.args))
        extra.args$z <- rep(1, length(extra.args$x))

    if(!"zcases" %in% names(extra.args))
        extra.args$zcases <- rep(1, length(extra.args$x))

    if(!"zcase.ids" %in% names(extra.args))
        extra.args$zcase.ids <- unique(extra.args$zcases)

    if ("groups" %in% names(extra.args)) {
        if ("group.args" %in% names(extra.args) && length(extra.args$group.args) > 
            0) {
            temp <- as.numeric(factor(extra.args$groups, levels = extra.args$group.ids))
            for (i in extra.args$group.args) {
                extra.args[[i]] <- extra.args[[i]][temp]
            }
        }
        extra.args$groups <- NULL
    }

    if ("zcases" %in% names(extra.args)) {
        if ("zcase.args" %in% names(extra.args) && length(extra.args$zcase.args) > 
            0) {
            temp <- as.numeric(factor(extra.args$zcases, levels = extra.args$zcase.ids))
            for (i in extra.args$zcase.args) {
                extra.args[[i]] <- extra.args[[i]][temp]
            }
        }
#        extra.args$zcases <- NULL
    }
    
    #reorder to make better layers
    #want the 'lay down' order to be 
    # first to last groups by rows 
    # so later full pies lay on 
    # earlier pies if there is any overlap

    zcase.count <- length(extra.args$zcase.ids)
    zcase.len <- max(sapply(extra.args$zcase.ids, function(x) length(extra.args$x[extra.args$zcases == 
        x])))

    ref <- unlist(lapply(1:zcase.len, function(i) seq(i, length(extra.args$x), zcase.len))) 
    temp <- unique(c(extra.args$panel.elements, extra.args$zcase.args, extra.args$group.args))
    temp <- temp[temp %in% names(extra.args)]

    for(i in temp){
        extra.args[[i]] <- extra.args[[i]][ref]
    }


    zcase.starts <- (((1:zcase.len) - 1)*zcase.count) + 1

    zcase.sums <- as.vector(sapply(zcase.starts, function(x) sum(na.omit(extra.args$z[x:(x+zcase.count-1)]))))
    
#print(zcase.sums)

#hard reset
    if("z.rowsum.lim" %in% names(extra.args))
        extra.args$zlim <- extra.args$z.rowsum.lim else range(zcase.sums)

#make nice pie for single z case
    if(length(extra.args$zcase.ids) < 2 && !"centre" %in% names(extra.args))
        extra.args$center <- FALSE

    temp <- unique(c(extra.args$panel.elements, extra.args$group.args, "zcase.ref", "angle", "start"))
    temp <- temp[!temp %in% "z"]
    temp2 <- unique(c("z", extra.args$zcase.args))

    for(k in 1:zcase.len){

        new <- extra.args
        i <- zcase.starts[k]
        for(j in temp){
            new[[j]] <- new[[j]][i]
        }
        for(j in temp2){
            new[[j]] <- extra.args[[j]][i:(i+zcase.count-1)]
        }
        
        if(!"radius" %in% names(new)){
            new$radius <- do.call(cexHandler, listUpdate(new, list(z=sum(new$z)))) 
        }

        if(!"angle" %in% names(new)){
            zz <- cumsum(new$z)
            zz <- zz/(zz[length(zz)]) * 360
            zz <- c(zz[1], diff(zz))
            new$angle <- zz
        }

        if(!"start" %in% names(new)){
            zz <- cumsum(new$angle)
            zz <- c(0, zz[-length(zz)])
            zz <- zz - (new$angle[1]/2)
            new$start <- zz
        }

        temp2 <- unique(c(temp2, "angle", "start"))

        for(j in 1:zcase.count){
        new2 <- new
            for(l in temp2){
                new2[[l]] <- new[[l]][j]
            }
            do.call(loaPieSegment, new2)          
        }        
        
    }

}






##############################
##############################
##panel.zcasePieSegmentPlot
##############################
##############################

panel.zcasePieSegmentPlot <- function (..., zcase.rescale = TRUE, 
                                       loa.settings = FALSE){

    ################################
    #panel.zcasePieSegmentPlot v0.2 
    #kr 
    ################################
    #draws the pieSegment plot
    ################################
    #notes
    #
    #faster revision of v0.1
    #revision as part of more flexible 
    #general panel handler for these plots
    ################################

    if (loa.settings) 
        return(list(zcase.args = c("col"),  
            default.settings = list(key.fun = "draw.zcasePlotKey", 
                grid = FALSE)))

    extra.args <- list(...)

#temp col fix
#needs work so borders and segment fill 
#transparency sepearated

    if("alpha" %in% names(extra.args)){
        extra.args$alpha.regions <- extra.args$alpha
        extra.args$alpha <- NULL
    }

    if(!"z" %in% names(extra.args))
        extra.args$z <- rep(1, length(extra.args$x))

    if(!"zcases" %in% names(extra.args))
        extra.args$zcases <- rep(1, length(extra.args$x))

    if(!"zcase.ids" %in% names(extra.args))
        extra.args$zcase.ids <- unique(extra.args$zcases)

    if ("groups" %in% names(extra.args)) {
        if ("group.args" %in% names(extra.args) && length(extra.args$group.args) > 
            0) {
            temp <- as.numeric(factor(extra.args$groups, levels = extra.args$group.ids))
            for (i in extra.args$group.args) {
                extra.args[[i]] <- extra.args[[i]][temp]
            }
        }
        extra.args$groups <- NULL
    }

    if ("zcases" %in% names(extra.args)) {
        if ("zcase.args" %in% names(extra.args) && length(extra.args$zcase.args) > 
            0) {
            temp <- as.numeric(factor(extra.args$zcases, levels = extra.args$zcase.ids))
            for (i in extra.args$zcase.args) {
                extra.args[[i]] <- extra.args[[i]][temp]
            }
        }
#        extra.args$zcases <- NULL
    }
    
    #reorder to make better layers
    #want the 'lay down' order to be 
    # first to last groups by rows 
    # so later full pies lay on 
    # earlier pies if there is any overlap


#make nice pie for single z case
    if(length(extra.args$zcase.ids) < 2 && !"centre" %in% names(extra.args))
        extra.args$center <- FALSE

    zcase.count <- length(extra.args$zcase.ids)
    zcase.len <- max(sapply(extra.args$zcase.ids, function(x) length(extra.args$x[extra.args$zcases == 
        x])))





    ref <- unlist(lapply(1:zcase.len, function(i) seq(i, length(extra.args$x), zcase.len))) 
    temp <- unique(c(extra.args$panel.elements, extra.args$zcase.args, extra.args$group.args))
    temp <- temp[temp %in% names(extra.args)]


    for(i in temp){
        extra.args[[i]] <- extra.args[[i]][ref]
    }

    #make numeric reference for zcase
    extra.args$zcase.ref <- as.numeric(extra.args$zcases)

    extra.args$angle <- rep(360/zcase.count, length(extra.args$x))
    extra.args$start <- extra.args$angle * (extra.args$zcase.ref - 1)
    extra.args$start <- extra.args$start - (extra.args$angle/2)

    temp <- unique(c(temp, "zcase.ref", "angle", "start"))

    for(i in 1:length(extra.args$x)){

        new <- extra.args
        for(j in temp){
            new[[j]] <- new[[j]][i]
        }

#this is a little messy but it
#works fine!


        if (!"radius" %in% names(new)) {
            n2 <- new
            if (zcase.rescale) {
                if ("zcase.zlim" %in% names(extra.args)) 
                  new$zlim <- new$zcase.zlim[[new$zcase.ref[1]]]
            }
            new$radius <- do.call(cexHandler, new)
        }

        do.call(loaPieSegment, new)          

    }

}




