#in development code
#[1 -TBC] functions 


#loa shapes 
#building blocks for plots
#loaPolygon
#loaPieSegment
#loaCircle

#######################
##might want to make 
##own space for conds
#######################


#urgent
#fixes 
#to dos
#suggestions

###########################
###########################
#loaPolygon
###########################
###########################


loaPolygon <- function(x, y, ..., polygon = NULL, loa.scale = NULL){

    #####################
    #loaPolygon v0.1
    #####################

    #might want to think about the polygon and loa.scale positions 
    #in formals?

    #current.panel.limits() gives the x and ylims

#could have relative.x and relative.y options???
#could give default col (via settings) and 
#                   no border using this????
#could give different default for polygon
#diamond, circle, etc.

    #new test
    #have to do this because lpolygon hardcoded
    #to border=TRUE="black"
    extra.args <- list(...)
    if(!"border" %in% names(extra.args))
        extra.args$border <- do.call(getPlotArgs, 
                                     listUpdate(extra.args, 
                                     list(defaults.as="plot.polygon")))$border

    if(!"col" %in% names(extra.args)){
        extra.args$col <- do.call(getPlotArgs, 
                                  listUpdate(extra.args, 
                                  list(defaults.as="plot.symbol")))$col
    }


    if(is.null(polygon))
        polygon = list(x=c(1, 1, -1, -1), 
                       y=c(1, -1, -1, 1))

#could use listLoad on this to use
#loa.scale.fit="absolute", etc.

#could make the loa.scale handling clever
#numeric sets list(scale)
#character sets list(fit)
#could also set the units of output for x,y
#             e.g. "cm", etc.

    if(is.null(loa.scale))
        loa.scale <- list()
    loa.scale <- listUpdate(list(fit="relative", scale = 1/50),
                            loa.scale)

    if(loa.scale$fit=="relative"){

        #rescale objects

        #covert x and y scales to npc

#        x <- as.numeric(grid::convertX(grid::unit(x, "native"), 
#            "npc"))
#        y <- as.numeric(grid::convertY(grid::unit(y, "native"), 
#            "npc"))
#        temp.fun <- function(x) x[2] - x[1]
#        if ("scale" %in% names(loa.scale)) {
#            polygon$x <- polygon$x * loa.scale$scale
#            polygon$y <- polygon$y * loa.scale$scale
#        }
#        temp <- sapply(current.panel.limits("mm"), temp.fun)
#        temp <- temp[2]/temp[1]
#        x <- x + (polygon$x * temp)
#        y <- y + polygon$y
#        x <- as.numeric(grid::convertX(grid::unit(x, "npc"), 
#            "native"))
#        y <- as.numeric(grid::convertY(grid::unit(y, "npc"), 
#            "native"))

        #tested simplification

        temp <- current.panel.limits("native")
        ref1 <- (temp[[1]][2]-temp[[1]][1])/50
        ref2 <- (temp[[2]][2]-temp[[2]][1])/50
        temp <- current.panel.limits("mm")
        ref3 <- (temp[[2]][2]-temp[[2]][1])/
                (temp[[1]][2]-temp[[1]][1])
        x <- x + (polygon$x * ref1 * ref3)
        y <- y + (polygon$y * ref2)

    } else {

        #absolute scale objects

        x <- x + polygon$x
        y <- y + polygon$y

    }

    do.call(lpolygon, listUpdate(list(x=x, y=y), extra.args))

}


#######################################
#loaCircle
#######################################


loaCircle <- function(..., polygon = NULL, radius = 1){

    extra.args <- list(...)

    #this is just loaPieSegement with 
    #center = FALSE, start=0, angle=360 forced

    extra.args <- listUpdate(extra.args, 
                             list(polygon=polygon, radius=radius, 
                                  center=FALSE, start=0, angle=360))
    do.call(loaPieSegment, extra.args)

}


######################################
#loaPieSegment
######################################


loaPieSegment <- function(..., polygon = NULL, start = 0, angle=360, radius = 1, center=TRUE){


#might want to rewrite this 
#using extra args???

    a <- start:(start+angle)
    
    myx <- radius * sin(pi * a/180) 
    myy <- radius * cos(pi * a/180) 

    if(center){
        myx <- c(0, myx, 0)
        myy <- c(0, myy, 0)
    }

    polygon <- list(x = myx,
                    y = myy)

    loaPolygon(..., polygon=polygon)

}


#testing

#a <- 1:10
#b <- a*10

#panel.temp <- function(x,y,...){
#for(i in 1:length(x))
#    loaCircle(x[i],y[i],...)
#}

#loaPlot(~a*b, panel=panel.temp, col="red", radius=0.3)

