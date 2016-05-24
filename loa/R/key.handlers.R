#in development code
#[TBC - NUMBER] functions 

#key Handling

#do properly


#urgent 
##################
#tidy all of these
#

#




###########################
###########################
#keyHandler
###########################
###########################

keyHandler <- function(key = NULL, ..., output = "key"){


    #set up
    extra.args <- list(...)

    #output == other.args
    if(output == "other.args")
        return(extra.args[!names(extra.args) %in% names(extra.args)[grep("^key", names(extra.args))]])

    #make key... args list
    temp <- grep("^key[.]", names(extra.args))
    key.args <- if(length(temp)>0){
                    temp <- extra.args[temp]
                    names(temp) <- gsub("^key[.]", "", names(temp))
                    temp  
                } else list()


#could add | option below or above 
#to stop these being done for a 
#key2legend option

    key <- if(is.list(key)) key else 
               if(is.null(key)) list() else 
                   if(is.logical(key) && key) list() else 
                       if(is.character(key) | is.function(key))
                           list(fun=key) else
                               FALSE


    if(isGood4LOA(key)){

        if(is.logical(key))
            key <- list()

        if(is.character(key) | is.function(key))
            key <- list(fun = key) else 
                if(!is.list(key)){
#need warning
                    key <- list()
                }

        key <- listUpdate(key, key.args)

        ##checks 

        #if key.position used this overrides 
        #space settings everywhere
        if("position" %in% names(key))
            key$space <- key$position 

        #if fun not there
        #put there
        if(!"fun" %in% names(key))
            key$fun <- "draw.loaPlotZKey"

    }

    if(output == "key")
        return(key)

    #make like legend
    if(isGood4LOA(key)){

    temp <- list(fun = key$fun,
                 args = list(key = list(), draw = FALSE))

    temp$args$key <- listUpdate(temp$args$key, key)

    temp$args$key <- listUpdate(temp$args$key, key, ignore = c(names(temp), "draw"))

    if(!"space" %in% names(temp$args$key))
         temp$args$key$space <- "right"

    temp <- list(right = temp)
    names(temp) <- temp$right$args$key$space

    key <- temp

    }    


    if(output == "legend") 
        return(key)


    key

}


##############################
##############################
##draw.loaColorKey
##############################
##############################


draw.loaColorKey <- function (key = NULL, draw = FALSE, vp = NULL, ...){

    if (!"at" %in% names(key)) 
        key$at <- seq(min(key$zlim), max(key$zlim), length.out = 100)

######################
#catch col and alpha 
#in draw.loaColorKey
######################

#    if(!"col" %in% names(key)){

        temp <- listUpdate(list(...), key)

#new bit testing

        if ("isolate.col.regions" %in% names(key)){

            key$col <- NULL
            key$alpha <- NULL
        }



        key$col <- do.call(colHandler, listUpdate(key, 
                                       list(z=temp$at, ref=temp$at)))
#    }

#####################
#simplified 
#testing
#####################

#    if("alpha.regions" %in%names(key)){
        key$alpha <- NULL
        key$alpha.regions <- NULL
#    }

#    if(!"alpha" %in% names(key)){
#        key$alpha <- do.call(colHandler, listUpdate(key, list(z = key$zlim, 
#        output = "all")))$alpha.regions
#        key$alpha <- key$alpha.regions
#        key$alpha.regions <- NULL
#    }

    draw.colorkey(key, draw, vp)

}






##############################
##############################
##draw.loaColorRegionsKey
##############################
##############################


draw.loaColorRegionsKey <- function (key = NULL, draw = FALSE, vp = NULL, ...) 
{
    if (!"at" %in% names(key)) 
        key$at <- pretty(c(min(key$zlim), max(key$zlim)))

    if ("isolate.col.regions" %in% names(key)) 
        key$col <- NULL

    if (!"col" %in% names(key)) {
        temp <- listUpdate(list(...), key)
        key$col <- colHandler(1:(length(key$at) - 1), col.regions = temp$col.regions, 
            output = "col")
    }
    key <- listUpdate(key, list(labels = list(at = key$at)))
    if (!"col" %in% names(key)) {
        key$col <- do.call(colHandler, listUpdate(key, list(z = key$zlim, 
            output = "all")))$col.regions
    }
    if (!"alpha" %in% names(key)) {
        key$alpha <- do.call(colHandler, listUpdate(key, list(z = key$zlim, 
            output = "all")))$alpha.regions
    }
    draw.colorkey(key, draw, vp)
}




######################################
######################################
##draw.zcasePlotKey
######################################
######################################


draw.zcasePlotKey <- function (key = NULL, draw = FALSE, vp = NULL, ...) 
{
    extra.args <- list(...)

    if (!is.list(key)) {
        warning("suspect key ignored", call. = FALSE)
        return(nullGrob())
    }

#new version 0.2.28
#test for turning off zcase key if only one colour set
#only currently applied by stackPlot as a ycases argument
#overridden there using force.key = TRUE

    if("zcases.key.method2" %in% names(key) && key$zcases.key.method2)
        if("col" %in% names(key) && length(key$col)<2)
            return(nullGrob())

    
    if(!"zcase.ids" %in% names(key)){
        key$zcase.ids <- if("zlab" %in% names(key))
                             key$zlab else " "
    }

    #taken from draw.loaPlotZKey
    #could simplify this

    z.main <- nullGrob()
    z.main.ht <- unit(0.01, "cm")
    zcases.main.ht <- z.main.ht
    zcases.main.wd <- z.main.ht
    zcases.elements <- z.main
    zcases.elements.ht <- z.main.ht
    zcases.elements.wd <- z.main.ht
    zcases.labels <- z.main
    zcases.labels.ht <- z.main.ht
    zcases.labels.wd <- z.main.ht

    temp <- key$panel.elements
    z.check <- if (is.null(temp)) 
        FALSE
    else if ("z" %in% temp) 
        TRUE
    else FALSE
    z.check <- FALSE

    z <-list()


#this is where zcases starts
#in the original code

    zcases.temp <- key[grep("^zcases[.]", names(key))]
    names(zcases.temp) <- gsub("^zcases[.]", "", names(zcases.temp))
    zcases <- if (!"zcases" %in% names(key)) 
        list()
    else if (is.list(key$zcases)) 
        key$zcases
    else if (is.logical(key$zcases) && key$zcases) 
        list()
    else FALSE


    zcases <- if (is.logical(zcases) && !zcases) 
        zcases
    else if (is.list(zcases)) 
        listUpdate(zcases, zcases.temp)
    else zcases.temp
    if (is.list(zcases) && !"at" %in% names(zcases)) 
        zcases$at <- key$zcase.ids
    if (is.list(zcases) && length(zcases) > 0) {


#could simplify this
#currently needs 
#key.zcases.main to reset name
#but key.main should be enough

        if (!"main" %in% names(zcases)) 
            zcases$main <- "zcases"
        if (!"labels" %in% names(zcases)) 
            zcases$labels <- if (is.null(zcases$at)) 
                NULL
            else as.character(zcases$at)

 
        zcases$col <- if ("col" %in% names(key)) 
                key$col
            else if (is.null(z$col)) 
                do.call(colHandler, listUpdate(key, list(z = NULL, 
                  ref = 1:length(zcases$labels))))
            else z$col[1]

#might not track if multiple cols are set 
#not zcases as set.
       
        zcases$border <- if("border" %in% names(key))
                    key$border else 
                    getPlotArgs("plot.polygon")$border

        
#        if (!"cex" %in% names(zcases)) 
#            zcases$cex <- if ("cex" %in% key$zcase.args) 
#                key$cex
#            else if (is.null(z$cex)) 
#                do.call(cexHandler, listUpdate(key, list(z = NULL, 
#                  ref = 1:length(zcases$labels))))
#            else z$cex[1]

#        if (!"pch" %in% names(zcases)) 
#            zcases$pch <- if ("pch" %in% key$zcase.args) 
#                key$pch
#            else if (is.null(z$pch)) 
#                do.call(pchHandler, listUpdate(key, list(z = NULL, 
#                  ref = 1:length(zcases$labels))))
#            else z$pch[1]


         #don't need these

         zcases$cex<-3
         zcases$pch<-15

        if (isGood4LOA(zcases$main)) {

##########################
#temp fix
#txt <- zcases$main[[1]][[1]]

            txt <- if (is.list(zcases$main)) 
                       zcases$main[[1]] else zcases$main

#to catch expressions as main
#could simply/merge 
#this and next bit
###########################
            
            temp <- if (is.list(zcases$main)) 
                listUpdate(list(cex = 1.1), zcases$main)
            else list(cex = 1.1)
            txt.settings <- getPlotArgs("axis.text", user.resets = temp)
            zcases.main <- textGrob(txt, gp = do.call(gpar, txt.settings))
            zcases.main.ht <- unit(2, "grobheight", data = list(zcases.main))
            zcases.main.wd <- unit(1.1, "grobwidth", data = list(zcases.main))
        }
        if (isGood4LOA(zcases$col) & isGood4LOA(zcases$cex) & 
            isGood4LOA(zcases$pch)) {

            len <- length(key$zcase.ids)

            temp <- rep(zcases$cex, length.out = len) * 0.8
            y <- rep(temp/2, each = 2)
            y <- y + 0.5
            if (isGood4LOA(zcases$labels)) {
                y[y < 1] <- 1
            }
            ht <- sum(y)
            y <- cumsum(y)
            y <- y[seq(1, length(y), 2)]/ht
            x <- rep(0.5, length(y))

#this should be tidied

             x1 <- rep(c(0.2, 0.2, 0.8, 0.8), len)
             y1 <- rep(c(0.2, 0.8, 0.8, 0.2), len)
             y1 <- y1 + rep(0:(len-1), each=4)
             y1 <- y1/(max(y1)+0.2)

            zcases.elements <- polygonGrob(x = x1, y = y1,
                  id.lengths=rep(4,len),  
                  default.units = "npc", gp = gpar(fill = rep(zcases$col, 
                  length.out = len), col =rep(zcases$border, len)))


#            zcases.elements <- pointsGrob(x = x, y = y, pch = rep(zcases$pch, 
#                length.out = len), size = unit(par("cex"), "char"), 
#                default.units = "npc", gp = gpar(col = rep(zcases$col, 
#                  length.out = len), cex = temp * 0.8))

            zcases.elements.ht <- unit(ht/4, "cm")
            zcases.elements.wd <- unit(max(zcases$cex)/4, "cm")

        }
        if (isGood4LOA(zcases$labels)) {
            txt <- if (is.list(zcases$labels)) 
                zcases$labels[[1]]
            else zcases$labels
            temp <- if (is.list(zcases$labels)) 
                listUpdate(list(cex = 1), zcases$labels)
            else list(cex = 1)
            txt.settings <- getPlotArgs("axis.text", user.resets = temp)
            zcases.labels <- textGrob(txt, x = 0, y = y, just = c("left", 
                "centre"), gp = do.call(gpar, txt.settings), 
                default.units = "npc")
            zcases.labels.ht <- unit(1.1, "grobheight", data = list(zcases.labels))
            zcases.labels.wd <- unit(1.1, "grobwidth", data = list(zcases.labels))
        }
    }


    scales.ht <- unit.c(zcases.main.ht, zcases.elements.ht)

#simplify
    temp1 <- max(zcases.elements.wd)
    temp2 <- max(zcases.labels.wd)
    temp3 <- max(zcases.main.wd)
    if (as.numeric(convertX(temp1 + temp2, "cm")) > as.numeric(convertX(temp3, 
        "cm"))) {
        scales.wd <- unit.c(temp1, temp2)
    }
    else {
        scales.wd <- unit.c(temp1, temp3 - temp2)
    }
    key.layout <- grid.layout(nrow = 2, ncol = 2, heights = scales.ht, 
        widths = scales.wd, respect = TRUE, just = "right")
    key.gf <- frameGrob(layout = key.layout, vp = vp)
    key.gf <- placeGrob(key.gf, zcases.main, row = 1, col = 1:2)
    key.gf <- placeGrob(key.gf, zcases.elements, row = 2, col = 1)
    key.gf <- placeGrob(key.gf, zcases.labels, row = 2, col = 2)
    key.gf
}




######################################
######################################
##draw.ycasePlotKey
######################################
######################################


draw.ycasePlotKey <- function (key = NULL, draw = FALSE, vp = NULL, ...) 
{
    extra.args <- list(...)

    if (!is.list(key)) {
        warning("suspect key ignored", call. = FALSE)
        return(nullGrob())
    }

#new to version 0.2.28
#test for turning off key if only one colour set
#only currently applied by stackPlot 
#overridden there using force.key = TRUE

    if("ycase.key.method2" %in% names(key) && key$ycase.key.method2)
        if("col" %in% names(key) && length(key$col)<2)
            return(nullGrob())

    if(!"ycases.main" %in% names(key))
        key$ycases.main <- "ycases"

    #cheat to use zcasePlot for ycases
    names(key) <- gsub("ycases", "zcases", names(key))
    names(extra.args) <- gsub("ycases", "zcases", names(extra.args))
    do.call(draw.zcasePlotKey, listUpdate(list(key = key, draw = draw, vp = vp), extra.args))  
    
}   












#############################################################
#############################################################
##draw.loaPlotZKey
#############################################################
#############################################################



draw.loaPlotZKey <- function (key = NULL, draw = FALSE, vp = NULL, ...){

    #############
    #setup
    #############

    #key is args from key
    #... 

##############################################
#the key.z.labels does not seem to track properly
#need to fix this
##############################################

    extra.args <- list(...) 

    #check key useable
    if (!is.list(key)){ 
        warning("suspect key ignored", call. = FALSE)
        return(nullGrob())
    }

###pch not tracked if set in call 
###but not a group.arg


#might not need some of these

    #default key components
    z.main <- nullGrob() 
    z.main.ht <- unit(0.01, "cm")
    z.main.wd <- z.main.ht
    z.elements <- z.main 
    z.elements.ht <- z.main.ht
    z.elements.wd <- z.main.ht
    z.labels <- z.main 
    z.labels.ht <- z.main.ht
    z.labels.wd <- z.main.ht

    groups.main <- z.main 
    groups.main.ht <- z.main.ht
    groups.main.wd <- z.main.ht
    groups.elements <- z.main 
    groups.elements.ht <- z.main.ht
    groups.elements.wd <- z.main.ht
    groups.labels <- z.main 
    groups.labels.ht <- z.main.ht
    groups.labels.wd <- z.main.ht

    zcases.main <- z.main 
    zcases.main.ht <- z.main.ht
    zcases.main.wd <- z.main.ht
    zcases.elements <- z.main 
    zcases.elements.ht <- z.main.ht
    zcases.elements.wd <- z.main.ht
    zcases.labels <- z.main 
    zcases.labels.ht <- z.main.ht
    zcases.labels.wd <- z.main.ht


    #check for z
    temp <- key$panel.elements
    z.check <- if(is.null(temp)) FALSE else 
                  if("z" %in% temp) TRUE else FALSE

    #check for z info
    z.temp <- key[grep("^z[.]", names(key))]
    names(z.temp) <- gsub("^z[.]", "", names(z.temp))
    
    #make z from inputs
    z <- if(!"z" %in% names(key)) list() else
             if(is.list(key$z)) key$z else 
                 if(is.logical(key$z) && key$z) list() else 
                     FALSE
    z <- if(is.logical(z) && !z) z else
             if(is.list(z)) listUpdate(z, z.temp) else z.temp  

    #add at if list and at not there 
    if(is.list(z) && !"at" %in% names(z))
        z$at <- if(is.null(key$zlim))
                    NULL else {
                       temp <- pretty(key$zlim,
                                      if(is.null(z$n.ticks)) 5 else z$n.ticks)
                       temp[temp >= min(key$zlim) & temp <= max(key$zlim)]
                   }


if(is.list(z) && length(z) > 0){


    #text if plot controls have been transposed
    #if so warning that these may not map on key

    temp <- c("col", "pch", "cex")
    check.trans <- temp[temp %in% key$panel.elements]
    if(length(check.trans)>0){
        warning(paste(check.trans, collapse=", "), 
                " transposed; key may not automatically track", call. = FALSE)
    }

    if(length(check.trans)>0){
        
        temp.fun <- function(x, d=key){
                        if(paste(x, "lim", sep="") %in% names(d)) return(d[[paste(x, "lim", sep="")]][1])
                        if(paste(x, "unique", sep="") %in% names(d)) return(d[[paste(x, "unique", sep="")]][1])
                    }
        temp <- lapply(check.trans, temp.fun)
        names(temp) <- check.trans 

##################
#this might trip up if col suppled by user as key arg
#not sure because has to be in group.elements
##################

        z <- listUpdate(z, temp)
    }



    if(!"main" %in% names(z))
         z$main <- key$zlab
    if(!"labels" %in% names(z))
         z$labels <- if(is.null(z$at)) NULL else as.character(z$at)

    if(!"col" %in% names(z)){
         test <- if("group.ids" %in% names(key) && "col" %in% names(key)) TRUE else
                     if("zcase.ids" %in% names(key) && "col" %in% names(key)) TRUE else FALSE
         z$col <- if(test){
                      temp <- key$col
                      temp <- temp[ceiling(length(temp)/2)]
                      do.call(colHandler, listUpdate(key, list(col=temp, ref = 1:length(z$at)))) 
                  } else {
                      if("at" %in% names(z))
                          do.call(colHandler, listUpdate(key, list(z=z$at, ref=1:length(z$at)))) else
                              NULL
                  }         
    }

    if(!"cex" %in% names(z)){
         test <- if("group.ids" %in% names(key) && "cex" %in% names(key)) TRUE else
                     if("zcase.ids" %in% names(key) && "cex" %in% names(key)) TRUE else FALSE
         z$cex <- if(test){
                      temp <- key$cex
                      temp <- temp[ceiling(length(temp)/2)]
                      do.call(cexHandler, listUpdate(key, list(cex=temp, ref = 1:length(z$at)))) 
                  } else {
                      if("at" %in% names(z))
                          do.call(cexHandler, listUpdate(key, list(z=z$at, ref=1:length(z$at)))) else
                             NULL
                  }         
    }

    if(!"pch" %in% names(z)){
         test <- if("group.ids" %in% names(key) && "pch" %in% names(key)) TRUE else
                     if("zcase.ids" %in% names(key) && "pch" %in% names(key)) TRUE else FALSE
         z$pch <- if(test){
                      temp <- key$pch
                      temp <- temp[ceiling(length(temp)/2)]
                      do.call(pchHandler, listUpdate(key, list(pch=temp, ref = 1:length(z$at)))) 
                  } else {
                      if("at" %in% names(z))
                          do.call(pchHandler, listUpdate(key, list(z=NULL, ref=1:length(z$at)))) else
                              NULL
                  }          
    }

    if(isGood4LOA(z$main)){
    #handle character vector or list
##########################
#temp fix
#as above
#        txt <- z$main[[1]][[1]]
        txt <- if(is.list(z$main)) z$main[[1]] else z$main
###########################

        temp <- if(is.list(z$main))
                listUpdate(list(cex = 1.1), z$main) else list(cex = 1.1)
        txt.settings <- getPlotArgs("axis.text", user.resets = temp)
        
        z.main <- textGrob(txt, gp = do.call(gpar, txt.settings))
        z.main.ht <- unit(2, "grobheight", data = list(z.main))
        z.main.wd <- unit(1.1, "grobwidth", data = list(z.main))
    } 

    if(isGood4LOA(z$col) & isGood4LOA(z$cex) & isGood4LOA(z$pch)){

        len <- max(c(length(z$cex), length(z$col), length(z$pch)))
        temp <- rep(z$cex, length.out=len) * 0.8
        y <- rep(temp/2, each=2) 
        y <- y + 0.5
        if(isGood4LOA(z$labels)){
            y[y<1] <- 1 #reset for small columns
        }
        ht <- sum(y)
        y <- cumsum(y)
        y <- y[seq(1, length(y), 2)] / ht 
        x <- rep(0.5, length(y))       


#some of this might not be needed
#rep(...)?

        z.elements <- pointsGrob(x = x, y = y, pch=rep(z$pch, length.out = len), size = unit(par("cex"), "char"), 
                                 default.units = "npc", gp = gpar(col = rep(z$col, length.out = len), 
                                 cex=temp*0.8))

        z.elements.ht <- unit(ht/4, "cm")
        z.elements.wd <- unit(max(z$cex)/4, "cm")

    } 


    if(isGood4LOA(z$labels)){


#tidy next bit latter
#if.list about txt and temp setting

#####################
#update pchHandler
#rescale cex
#rescale spacing when cex small
#in elements
#####################
#transpose check
#group main, element, labels
#etc
#####################

        txt <- if(is.list(z$labels)) z$labels[[1]] else z$labels
        temp <- if(is.list(z$labels))
                listUpdate(list(cex = 1), z$labels) else list(cex = 1)
        txt.settings <- getPlotArgs("axis.text", user.resets = temp)
        
        z.labels <- textGrob(txt, x = 0, y = y, just = c("left", "centre"), 
                             gp = do.call(gpar, txt.settings), default.units = "npc")
        z.labels.ht <- unit(1.1, "grobheight", data=list(z.labels))
        z.labels.wd <- unit(1.1, "grobwidth", data=list(z.labels))

    } 

}



#groups info


#groups pch and cex missing
#if z not set



    #check for groups info
    groups.temp <- key[grep("^groups[.]", names(key))]
    names(groups.temp) <- gsub("^groups[.]", "", names(groups.temp))
    
    #make groups from inputs
    groups <- if(!"groups" %in% names(key)) list() else
                  if(is.list(key$groups)) key$groups else 
                      if(is.logical(key$groups) && key$groups) list() else 
                          FALSE
    groups <- if(is.logical(groups) && !groups) groups else
                  if(is.list(groups)) listUpdate(groups, groups.temp) else groups.temp  

    #add labels if list and at not there 
    if(is.list(groups) && !"at" %in% names(groups))
        groups$at <- key$group.ids

if(is.list(groups) && length(groups) > 0){

    if(!"main" %in% names(groups))
         groups$main <- "groups"
    if(!"labels" %in% names(groups))
         groups$labels <- if(is.null(groups$at)) NULL else as.character(groups$at)

    if(!"col" %in% names(groups)) 
         groups$col <- if("col" %in% key$group.args) 
                           do.call(colHandler, listUpdate(key, list(z = NULL, 
                                                          ref = 1:length(key$col)))) 
                               else if(is.null(z$col))
                                   do.call(colHandler, listUpdate(key, list(z=NULL, ref=1:length(groups$labels)))) else
                                       z$col[1]
    if(!"cex" %in% names(groups)) 
         groups$cex <- if("cex" %in% key$group.args) 
                           key$cex else if(is.null(z$cex))
                               do.call(cexHandler, listUpdate(key, list(z=NULL, ref=1:length(groups$labels)))) else
                                   z$cex[1]
    if(!"pch" %in% names(groups)) 
         groups$pch <- if("pch" %in% key$group.args) 
                           key$pch else if(is.null(z$pch))
                               do.call(pchHandler, listUpdate(key, list(z=NULL, ref=1:length(groups$labels)))) else
                                   z$pch[1]

    if(isGood4LOA(groups$main)){
    #handle character vector or list
        txt <- groups$main[[1]][[1]]
        temp <- if(is.list(groups$main))
                listUpdate(list(cex = 1.1), groups$main) else list(cex = 1.1)
        txt.settings <- getPlotArgs("axis.text", user.resets = temp)
        
        groups.main <- textGrob(txt, gp = do.call(gpar, txt.settings))
        groups.main.ht <- unit(2, "grobheight", data = list(groups.main))
        groups.main.wd <- unit(1.1, "grobwidth", data = list(groups.main))
    } 

    if(isGood4LOA(groups$col) & isGood4LOA(groups$cex) & isGood4LOA(groups$pch)){

        len <- max(c(length(groups$cex), length(groups$col), length(groups$pch)))
        temp <- rep(groups$cex, length.out=len) * 0.8
        y <- rep(temp/2, each=2) 
        y <- y + 0.5
        if(isGood4LOA(groups$labels)){
            y[y<1] <- 1 #reset for small columns
        }
        ht <- sum(y)
        y <- cumsum(y)
        y <- y[seq(1, length(y), 2)] / ht 
        x <- rep(0.5, length(y))       


#some of this might not be needed
#rep(...)?

        groups.elements <- pointsGrob(x = x, y = y, pch=rep(groups$pch, length.out = len), size = unit(par("cex"), "char"), 
                                 default.units = "npc", gp = gpar(col = rep(groups$col, length.out = len), 
                                 cex=temp*0.8))

        groups.elements.ht <- unit(ht/4, "cm")
        groups.elements.wd <- unit(max(groups$cex)/4, "cm")

    } 


    if(isGood4LOA(groups$labels)){


#tidy next bit latter
#if.list about txt and temp setting

#####################
#update pchHandler
#rescale cex
#rescale spacing when cex small
#in elements
#####################
#transpose check
#group main, element, labels
#etc
#####################

        txt <- if(is.list(groups$labels)) groups$labels[[1]] else groups$labels
        temp <- if(is.list(groups$labels))
                listUpdate(list(cex = 1), groups$labels) else list(cex = 1)
        txt.settings <- getPlotArgs("axis.text", user.resets = temp)
        
        groups.labels <- textGrob(txt, x = 0, y = y, just = c("left", "centre"), 
                             gp = do.call(gpar, txt.settings), default.units = "npc")
        groups.labels.ht <- unit(1.1, "grobheight", data=list(groups.labels))
        groups.labels.wd <- unit(1.1, "grobwidth", data=list(groups.labels))

    } 


##########################

}



##################################
##################################
##################################

#zcases info

#zcase pch and cex missing
#if z not set



    #check for zcases info
    zcases.temp <- key[grep("^zcases[.]", names(key))]
    names(zcases.temp) <- gsub("^zcases[.]", "", names(zcases.temp))
    
    #make zcases from inputs
    zcases <- if(!"zcases" %in% names(key)) list() else
                  if(is.list(key$zcases)) key$zcases else 
                      if(is.logical(key$zcases) && key$zcases) list() else 
                          FALSE
    zcases <- if(is.logical(zcases) && !zcases) zcases else
                  if(is.list(zcases)) listUpdate(zcases, zcases.temp) else zcases.temp  

    #add labels if list and at not there 
    if(is.list(zcases) && !"at" %in% names(zcases))
        zcases$at <- key$zcase.ids

if(is.list(zcases) && length(zcases) > 0){

    if(!"main" %in% names(zcases))
         zcases$main <- "zcases"
    if(!"labels" %in% names(zcases))
         zcases$labels <- if(is.null(zcases$at)) NULL else as.character(zcases$at)

    if(!"col" %in% names(zcases)) 
         zcases$col <- if("col" %in% key$zcase.args) 
                           do.call(colHandler, listUpdate(key, list(z = NULL, 
                                                          ref = 1:length(key$col)))) 
                               else if(is.null(z$col))
                                   do.call(colHandler, listUpdate(key, list(z=NULL, ref=1:length(zcases$labels)))) else
                                       z$col[1]
    if(!"cex" %in% names(zcases)) 
         zcases$cex <- if("cex" %in% key$zcase.args) 
                           key$cex else if(is.null(z$cex))
                               do.call(cexHandler, listUpdate(key, list(z=NULL, ref=1:length(zcases$labels)))) else
                                   z$cex[1]
    if(!"pch" %in% names(zcases)) 
         zcases$pch <- if("pch" %in% key$zcase.args) 
                           key$pch else if(is.null(z$pch))
                               do.call(pchHandler, listUpdate(key, list(z=NULL, ref=1:length(zcases$labels)))) else
                                   z$pch[1]

    if(isGood4LOA(zcases$main)){
    #handle character vector or list
        txt <- zcases$main[[1]][[1]]
        temp <- if(is.list(zcases$main))
                listUpdate(list(cex = 1.1), zcases$main) else list(cex = 1.1)
        txt.settings <- getPlotArgs("axis.text", user.resets = temp)
        
        zcases.main <- textGrob(txt, gp = do.call(gpar, txt.settings))
        zcases.main.ht <- unit(2, "grobheight", data = list(zcases.main))
        zcases.main.wd <- unit(1.1, "grobwidth", data = list(zcases.main))
    } 

    if(isGood4LOA(zcases$col) & isGood4LOA(zcases$cex) & isGood4LOA(zcases$pch)){

        len <- max(c(length(zcases$cex), length(zcases$col), length(zcases$pch)))
        temp <- rep(zcases$cex, length.out=len) * 0.8
        y <- rep(temp/2, each=2) 
        y <- y + 0.5
        if(isGood4LOA(zcases$labels)){
            y[y<1] <- 1 #reset for small columns
        }
        ht <- sum(y)
        y <- cumsum(y)
        y <- y[seq(1, length(y), 2)] / ht 
        x <- rep(0.5, length(y))       


#some of this might not be needed
#rep(...)?

        zcases.elements <- pointsGrob(x = x, y = y, pch=rep(zcases$pch, length.out = len), size = unit(par("cex"), "char"), 
                                 default.units = "npc", gp = gpar(col = rep(zcases$col, length.out = len), 
                                 cex=temp*0.8))

        zcases.elements.ht <- unit(ht/4, "cm")
        zcases.elements.wd <- unit(max(zcases$cex)/4, "cm")

    } 


    if(isGood4LOA(zcases$labels)){


#tidy next bit latter
#if.list about txt and temp setting

#####################
#update pchHandler
#rescale cex
#rescale spacing when cex small
#in elements
#####################
#transpose check
#group main, element, labels
#etc
#####################

        txt <- if(is.list(zcases$labels)) zcases$labels[[1]] else zcases$labels
        temp <- if(is.list(zcases$labels))
                listUpdate(list(cex = 1), zcases$labels) else list(cex = 1)
        txt.settings <- getPlotArgs("axis.text", user.resets = temp)
        
        zcases.labels <- textGrob(txt, x = 0, y = y, just = c("left", "centre"), 
                             gp = do.call(gpar, txt.settings), default.units = "npc")
        zcases.labels.ht <- unit(1.1, "grobheight", data=list(zcases.labels))
        zcases.labels.wd <- unit(1.1, "grobwidth", data=list(zcases.labels))

    } 


##########################

}

###############################
###############################
###############################
















    ###############
    #key output
    ###############

    scales.ht <- unit.c(z.main.ht, z.elements.ht, groups.main.ht, groups.elements.ht, zcases.main.ht, zcases.elements.ht)

##################
#this needs fixing
#wd need to know all
##################

    temp1 <- max(z.elements.wd, groups.elements.wd, zcases.elements.wd) 
    temp2 <- max(z.labels.wd, groups.labels.wd, zcases.labels.wd)    
    temp3 <- max(z.main.wd, groups.main.wd, zcases.main.wd)
    
    if(as.numeric(convertX(temp1 + temp2, "cm")) > as.numeric(convertX(temp3, "cm"))){
        scales.wd <- unit.c(temp1, temp2)
    } else {
#testing changes - temp2 to - temp1
        scales.wd <- unit.c(temp1, temp3 - temp1)
    }
    
###    scales.wd <- unit.c(z.elements.wd, z.labels.wd)

    key.layout <- grid.layout(nrow = 6, ncol = 2,
                              heights = scales.ht, 
                              widths = scales.wd, 
                              respect = TRUE, just="right")
    key.gf <- frameGrob(layout = key.layout, vp = vp)

    key.gf <- placeGrob(key.gf, z.main, row = 1, col = 1:2)
    key.gf <- placeGrob(key.gf, z.elements, row = 2, col = 1)
    key.gf <- placeGrob(key.gf, z.labels, row = 2, col = 2)

    key.gf <- placeGrob(key.gf, groups.main, row = 3, col = 1:2)
    key.gf <- placeGrob(key.gf, groups.elements, row = 4, col = 1)
    key.gf <- placeGrob(key.gf, groups.labels, row = 4, col = 2)

    key.gf <- placeGrob(key.gf, zcases.main, row = 5, col = 1:2)
    key.gf <- placeGrob(key.gf, zcases.elements, row = 6, col = 1)
    key.gf <- placeGrob(key.gf, zcases.labels, row = 6, col = 2)

    key.gf

}



######################################
#temp keys
######################################


draw.key.log10 <- function (key = NULL, draw = FALSE, vp = NULL, ...) {
    if (!"at" %in% names(key))
        key$at <- seq(min(key$zlim), max(key$zlim), length.out = 100)
   ticks <- if("tick.number" %in% names(key))
                     key$tick.number else 5
    if(!"labels" %in% names(key)){
        temp <- logTicks(10^c(min(key$zlim), max(key$zlim)), 1:9)
        temp <- temp[log10(temp)>=min(key$at) & log10(temp)<=max(key$at)]
        key$labels$at <- log10(temp)
        temp2 <- logTicks(10^c(min(key$zlim), max(key$zlim)), 1)
        key$labels$labels <- ifelse(temp %in% temp2, temp, "")
    }
   draw.loaColorKey(key = key, draw = draw, vp = vp, ...) }



