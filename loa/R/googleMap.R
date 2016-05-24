#in development code
#[8 -TBC] functions 

#GoogleMap - main function
#GoogleMap.old
#googlemap 
#quickMap - simplified version as code demo
#makeMapArg - map maker and modifier for above
#getMapArg - recover map from previous GoogleMap plot
#panel.GoogleMapsRaster - panels
#panel.GoogleMaps
#xscale.components.GoogleMaps - axis handlers
#yscale.components.GoogleMaps
#axis.components.GoogleMaps


##re GoogleMapsPlot
##GoogleMapsPlot(d, "lat", "lon", maptype="roadmap", 
##               path="&style=feature:water%7Chue:0xffffff|&style=feature:poi%7Cvisibility:off|&style=feature:all%7Celement:labels%7Cvisibility:off", 
##               map.cols=c("bisque", "white"))


#NEEDS urgent work

#GoogleMap code needs tidy

#IMPORTANT
#This needs archiving and tidying...


#GoogleMap - z, colorkey, etc not sort
#makeMapArg - recolor.map 
#           - think I have done this better elsewhere
#           - also need to see if this could be done quicker 
#             via the Google Maps call

#rename getMap because of GetMap in RgoogleMaps


###########################
###########################
#GoogleMap
###########################
###########################

googleMap <- function(...) GoogleMap(...)


GoogleMap <- function (x, data = NULL, panel = panel.loaPlot, map = NULL, 
    map.panel = panel.GoogleMapsRaster, recolor.map = FALSE, 
    ..., lon.lat = FALSE) {

#new rewrite using loaPlot 
#note changes
#z ~ lon * lat | cond
#

    extra.args <- list(...)

    #formular default z ~ lat * lon | cond
    if (!lon.lat) 
        if (!"formula.type" %in% names(extra.args)) 
            extra.args$formula.type = "z~y*x|cond"

    #make map based on data range
    temp <- do.call(formulaHandler, listUpdate(extra.args, list(x=x, data=data, output="lattice.like")))
    myx <- if("xlim" %in% names(extra.args))
                extra.args$xlim else temp$x
    myy <- if("ylim" %in% names(extra.args))
                extra.args$ylim else temp$y
    if (is.null(map)) {
        temp <- list(xlim = myx, ylim = myy, 
            recolor.map = recolor.map, aspect = NULL)
        temp <- listUpdate(temp, extra.args)
        map <- do.call(makeMapArg, temp)
    }

    #plot basic plot
    ans <- do.call(loaPlot, listUpdate(list(x = x, data = data, 
        panel = panel), extra.args))

    #insert and setup map
    ans$aspect.ratio <- map$aspect
    ans$panel.args.common$xlim <- map$xlim
    ans$x.limits <- map$xlim
    ans$panel.args.common$ylim <- map$ylim
    ans$y.limits <- map$ylim

    #update xy elements
    #think about better handling
    #think about doing in panelPal
    #and passing coversion program so all similar
    #think about keeping lat,lon?
    if(is.null(ans$panel.args.common$x.elements))
        ans$panel.args.common$x.elements <- "x"
    if(is.null(ans$panel.args.common$y.elements))
        ans$panel.args.common$y.elements <- gsub("x", "y", ans$panel.args.common$x.elements)
    for (i in 1:length(ans$panel.args)) {

        for(j in 1:length(ans$panel.args.common$y.elements)){

            temp <- RgoogleMaps::LatLon2XY.centered(map, 
                ans$panel.args[[i]][[ans$panel.args.common$y.elements[j]]], 
                ans$panel.args[[i]][[ans$panel.args.common$x.elements[j]]])
            ans$panel.args[[i]][[ans$panel.args.common$y.elements[j]]] <- temp$newY
            ans$panel.args[[i]][[ans$panel.args.common$x.elements[j]]] <- temp$newX

        }
    }
    panel <- ans$panel
    panel.with.map <- function(...) {
        map.panel(map)
        panel(...)
    }
    map.axis.comps <- axis.components.GoogleMaps(map)
    map.axis <- function(components, ...) axis.default(components = map.axis.comps, 
        ...)
    ans <- update(ans, panel = panel.with.map, aspect = map$aspect, 
        axis = map.axis)
    ans$panel.args.common$map <- map

    #output plot
    return(ans)
}






#####################################
#####################################





GoogleMap.old <- function(x, data = NULL, map = NULL,
    map.panel = panel.GoogleMapsRaster, 
    panel = panel.xyplot, 
    recolor.map = FALSE, ...){

    #noaa.sunrise GoogleMapsPlot revision
    #karl 2011-12-14

    ##############
    #notes 
    ##############
    #use US spelling - more consistent with 'parent' lattice.
    #if using colorkey as arg need col in there
    #(or col arg gets assigned as colorkey if set!)
    #map.panel for later
    #panel = panel.levelplot 
    #

    ##################
    #to do
    ##################
    #
    #
    #
    ###################
    #better legend/key control blocks or scale handling
    #tidier subscript handling
    #rethink aspect
    #

    ##################
    #to think about
    ##################
    #sort data before map?
    #if problem best to stop before getting map?
    ##################
    #currently ~ lat * lon
    #might want to rethink as ~ lon * lat
    #

    #####################
    #uses 
    #####################
    #lattice levelplot, latticeParseFormula
    #RgoogleMaps MapBackground
    #RColorBrewer (recoloring?)

    #not needed 'cos installation requirements

#    stopifnot(require("grid"))
#    stopifnot(require("lattice"))
#    stopifnot(require("RgoogleMaps"))
#    stopifnot(require("RColorBrewer"))

    ####################
    #setups
    ####################

    #extra.args
    extra.args <- list(...)
    extra.args$recolor.map <- recolor.map
 
    #get lat lon for formula
    d1 <- try(latticeParseFormula(x, data, dimension = 3, 
                                  multiple = TRUE),
              silent = TRUE)
    if(is(d1)[1] == "try-error")
        stop("problem with x/data combination", call. = FALSE)

    ######################
    #set plot defaults 
    #######################
    
    #add x/ylims, aspect, pch
    #pass to extra.args 
    #add recolor.map 
    #user listUpdate to user update

    #note: these are taken from z ~ x * y 
    #      so get flipped

#################
#this next bit could be simplified
#################

    temp <- list(xlim = d1$right.y, ylim = d1$right.x,
                 aspect = NULL, pch = 20)
    extra.args <- listUpdate(temp, extra.args)

    ###############
    #set up cex and col
    #using z ...handlers
    ###############

#################
#add ref just for this
#################

#but still an issue about wrapping

#both cex and col handlers need ref in
#then need a way to scale data smaller than x/y range
#to wrap terms like pch

#add ref = d1$right.y to list

#    extra.args <- listUpdate(list(z = d1$left), extra.args)
    extra.args <- listUpdate(list(z = d1$left, ref = d1$right.x), 
                             extra.args)

    extra.args$cex <- do.call(cexHandler, extra.args)
    extra.args <- listUpdate(extra.args, 
                             do.call(colHandler, listUpdate(extra.args, list(output="all"))))

###############
    
    ###############
    #get map
    ###############

    #get map if not suppled
    if(is.null(map))
        map <- do.call(makeMapArg, extra.args)

    #update xlim, ylim, aspect, etc because 
    #may not be exactly what was asked for
    extra.args$aspect <- map$aspect
    extra.args$xlim <- map$xlim
    extra.args$ylim <- map$ylim


    ##############################
    #rescale axis
    ##############################

    #scale axis for map projection
    map.axis.comps <- axis.components.GoogleMaps(map)
    map.axis <- function(components, ...) 
                   axis.default(components = map.axis.comps, ...)

    ############################
    #rescale data
    ############################

    #scale data for map projection
    temp <- LatLon2XY.centered(map, d1$right.x, d1$right.y)
    ..lat.rescale <- temp$newY
    ..lon.rescale <- temp$newX

    #rearrange formula for xyplot
    #keeping any conditioning
    temp <- as.character(x)
    temp <- temp[length(temp)]
    temp <- strsplit(temp, "[|]")
    x <- "..lat.rescale~..lon.rescale"
    if(length(temp[[1]])>1)
        x <- paste(x, temp[[1]][2], sep ="|")
    x <- as.formula(x)


#rescale wrappable data
#temp <- extra.args[names(extra.args) %in% c("z", "pch", "col", "cex")]

#expand problem cases
##    extra.args <- lapply(extra.args, function(x){
##                            x <- if(is.vector(x)){
##                                if(length(x) > 1 & length(x) < length(extra.args$ref))
##                                    rep(x, length(extra.args$ref))[1:length(extra.args$ref)] else x
##                            } else x
##                         })

    extra.args <- listUpdate(extra.args, 
                      listExpand(extra.args, ignore=c("xlim", "ylim", "at", "col.regions", 
                                                      "aspect", "ref", "layout"), 
                          ref = extra.args$ref)
                  )

    #########################
    #plot data
    #########################

    #use update method to allow
    #user fine control
    temp <- list(x = x, data = data, 
                 xlab = d1$right.y.name,
                 ylab = d1$right.x.name,
                 axis = map.axis,
##                 panel = function(x, y, subscripts, at, col.regions, ...){
##                             map.panel(map)
##
##                             #############
##                             #got to be better way of doing this
##                            #make all panel-specific
##                             temp <- list(...)
##                             if(!is.null(subscripts)){
##                                 temp <- lapply(temp, function(x)
##
#####################
#possible fix
#####################
##
##if length x < the real x!!!
##extrapolated?
##or supply pch in update?
##not convinced col.regions does anything in update
## 
##                                             x <- if(length(x)>1) x[subscripts] else x )
##                                 subscripts <- 1:length(subscripts)
##                             }
##                             temp <- listUpdate(
##                                         list(x = x, y = y, z = temp$z, at = at, 
##                                              col.regions=col.regions, subscripts=subscripts),
##                                              temp)
##                                         do.call(panel, temp)
##                             })
       panel = function(..., subscripts){
                   map.panel(map)
                   panel = panelPal.old(..., subscripts=subscripts, panel=panel)
               }

       )

    extra.args <- listUpdate(temp, extra.args)
    ans <- do.call(xyplot, extra.args) 
   
    ans$panel.args.common$map <- map

    #messy but means <-... still updates
    #might rethink
    plot(ans)
    invisible(ans)       

}










###########################
###########################
#quickMap
###########################
###########################

quickMap <- function(lat, lon, show.data = FALSE, ...){

    #get map
    map <- makeMapArg(lat, lon, ...)

    #scale axis for map projection
    map.axis.comps <- axis.components.GoogleMaps(map)
    map.axis <- function(components, ...) 
                   axis.default(components = map.axis.comps, ...)

    #scale data for map projection
    temp <- LatLon2XY.centered(map, lat, lon)
    lat <- temp$newY
    lon <- temp$newX

    #plot data on map
    xyplot(lat~lon, 
           xlim = map$xlim, ylim =map$ylim,
           aspect = map$aspect, 
           axis = map.axis,
           panel = function(...){
               panel.GoogleMapsRaster(map)
               if(show.data)
                   panel.xyplot(...)
           }, ...)
}




#########################
#########################
##makeMapArg
#########################
#########################


makeMapArg <- function(ylim, xlim, 
     aspect = NULL, recolor.map = FALSE, 
     ...
){

    #RgoogleMaps wrapper 
    #make maps are easier for panel work 
    #karl 2011-12-14

    ##############
    #notes 
    ##############
    #

    ##################
    #to think about
    ##################
    #move local function localMapManager in code body
    #

    #####################
    #uses 
    #####################

    #extra.args
    extra.args <- list(...)

    #local function
    localMapManager <- function(map){

#################
#destfile name fix needed?
#################

        #######################
        #native raster handler
        #######################

#######################
##        if("nativeRaster" %in% class(map$myTile) & require(png)){
##dropped require png 
##while testing namespace change
######################

        if("nativeRaster" %in% class(map$myTile)){
            #do to png native output
            writePNG(map$myTile, "XtempX.png")
            map$myTile <- readPNG("XtempX.png", native = FALSE)
            attr(map$myTile, "type") <- "rgb"
        }

        #set up
        ra <- dim(map$myTile)
    
        #version test
        tempfun <- function(x, pck = "RgoogleMaps") 
                      compareVersion(packageDescription(pck)$Version, x)

        #if RgoogleMaps version between 1.1.5 and 1.1.9.13 
        #currently don't know structure
        if(tempfun("1.1.9.13") < 0){
            warning(paste("GoogleMaps may not be able to support this version of 'RgoogleMaps'",
                          "\n\t[You may encounter problems]", 
                          "\n\t[If so, suggest updating RgoogleMaps]",
                    sep=""), call.=FALSE)

            #NOTE: this assumes
            #imagematrix, png
            #NOT tested for jgp
            map$myTile <- matrix(attr(map$myTile, "COL")[map$myTile],
                                 nrow = ra[1], ncol = ra[2]    
                          )
            map$myTile <- t(map$myTile)
            map$myTile <- map$myTile[nrow(map$myTile):1,]
            attr(map$myTile, "type") <- "local"
            return(map)
        }

        if(length(ra) > 2){
            if(ra[3] == 4 & attr(map$myTile, "type") == "rgb"){
                map$myTile <- rgb(map$myTile[, , 1], map$myTile[, , 2], 
                                  map$myTile[, , 3], map$myTile[, , 4])
                dim(map$myTile) <- ra[1:2]
                attr(map$myTile, "type") <- "local"
                return(map)
            }

            if(ra[3] == 3 & attr(map$myTile, "type") == "rgb"){
                map$myTile <- rgb(map$myTile[, , 1], map$myTile[, , 2], 
                                  map$myTile[, , 3])
                dim(map$myTile) <- ra[1:2]
                attr(map$myTile, "type") <- "local"
                return(map)
            }

            if(ra[3] == 1 & attr(map$myTile, "type") == "grey"){
                map$myTile <- grey(map$myTile[, , 1])
                dim(map$myTile) <- ra[1:2]
                attr(map$myTile, "type") <- "local"
                return(map)
            }
        } 

        if(length(ra) == 2){

            if(is.character(attr(map$myTile, "type")) && attr(map$myTile, "type") == "grey"){
                map$myTile <- grey(map$myTile[,])
                dim(map$myTile) <- ra[1:2]
                attr(map$myTile, "type") <- "openair"
                return(map)
            }
        }

        warning(paste("GoogleMaps encountered unexpected 'RgoogleMaps' output",
                      "\n\t[You may encounter problems or some options may not be supported]", 
                      "\n\t[If so, suggest updating RgoogleMaps]",
                sep=""), call.=FALSE)
        return(map)
    }

    #####################
    #main
    #####################

    #get names of args that MapBackground can handle
    temp <- unique(c(names(formals(MapBackground)), 
                     names(formals(GetMap.bbox)),
                     names(formals(GetMap))))

    #get suitable ranges
    temp2 <- try(qbbox(lat = ylim, lon = xlim), silent = TRUE)
    if(is(temp2)[1] == "try-error")
        stop(paste("\tCould not apply supplied lat, lon combination",
                   "\n\t[check call settings and data source]", sep = ""),
             call.=FALSE)

    #set size
    my.y <- diff(range(temp2$latR, na.rm=TRUE))
    my.x <- diff(range(temp2$lonR, na.rm=TRUE))

    #size was c(640, 640)
    my.size <- if(my.y > my.x)
                   c(ceiling((my.x/my.y) * 640), 640) else 
                       c(640, ceiling((my.y/my.x) * 640))

    #override some RgoogleMaps defaults
    map <- list(lon = temp2$lonR, lat = temp2$latR, destfile = "XtempX.png",
                 maptype = "terrain", size = my.size)

    if(my.x==0 & my.y==0){
        if(is.null(map$zoom))
            map$zoom <- 15
        map$size <- c(640,640)
    }

    if(any(is.na(map$size)))
        map$size[is.na(map$size)] <- 64
    map$size[map$size < 1] <- 64

    ##update my defaults with relevant ones in call
    map <- listUpdate(map, extra.args, use.b = temp)

    #use MapBackground and list of allowed args
    map <- try(do.call(GetMap.bbox, map), silent = TRUE)

    if(is(map)[1] == "try-error")
        stop(paste("\tCould not apply supplied lat, lon and RgoogleMap args",
                   "\n\t[check call settings and data source]", sep = ""),
             call.=FALSE)



#    #temp til RgoogleMaps Fixed
#
#    temp2 <- try(qbbox(lat = ylim, lon = xlim), silent = TRUE)
#    if(is(temp2)[1] == "try-error")
#        stop("could not apply supplied lat, lon combination", call.=FALSE)

#    #override some other RgoogleMaps defaults
#    #map <- list(lon = temp2$lonR, lat = temp2$latR, destfile = "XtempX.png", 
#    #            maptype = "terrain")

#    #fix
#    #get square range
#    temp.2 <- c(min(ylim, na.rm = TRUE) + (diff(range(ylim, na.rm=TRUE))/2),
#                min(xlim, na.rm = TRUE) + (diff(range(xlim, na.rm=TRUE))/2))
#    zoom <- min(MaxZoom(range(ylim, na.rm=TRUE), 
#                        range(xlim, na.rm=TRUE)), 
#                na.rm=TRUE)
#    if(is.infinite(zoom))
#        zoom <- 15
#    map <- list(center = temp.2, zoom = zoom, destfile = "XtempX.png", 
#                maptype = "terrain")

#    #set best size for map
#    #NOTE: superseded by extra.args if user forced
#    size <- diff(temp2$lonR)/diff(temp2$latR)
#    size <- sqrt(size^2)
#    size <- c(size, 1/size)
#    size <- ceiling((size/max(size)) * 640)

#    #force size until fixed
#    #size <- c(640, 640)

#    map$size <- size 

#    ##update my defaults with relevant ones in call
#    map <- listUpdate(map, extra.args, use.b = temp)

#    #use MapBackground and list of allowed args
#    #map <- try(do.call(MapBackground, map), silent = TRUE)
#    map <- try(do.call(GetMap, map), silent = TRUE)

#    if(is(map)[1] == "try-error")
#        stop("problem with supplied lat, lon or RgoogleMap args.", call. = FALSE)
        
#    #save size for later
#    map$size <- size

    #get ranges for panel
    #use RgoogleMaps function
    temp <- LatLon2XY.centered(map, c(map$BBOX$ll[1], map$BBOX$ur[1]), 
                                    c(map$BBOX$ll[2], map$BBOX$ur[2]))
    map$xlim <- temp$newX
    map$ylim <- temp$newY    

    ####################
    #map management
    ####################
    map <- localMapManager(map)

    #aspect
    #fix this when fixed
    if(is.null(aspect)){
        aspect <- diff(map$ylim)/diff(map$xlim) 
        aspect <- sqrt(aspect^2)
    }
    map$aspect <- aspect

    ##############
    #recolor map 
    #############

    if(is.logical(recolor.map) && recolor.map)
        recolor.map <- c("white", "grey")

    if(!is.null(recolor.map) & !is.logical(recolor.map)){

        ra <- dim(map$myTile)

        #if single expand to col range if possible
        if(length(recolor.map) == 1)
             recolor.map <- brewer.pal(9, recolor.map)

        #make an intensity scale
        temp <- apply(col2rgb(map$myTile), 2, prod)

        #reset cols in frame
        map$myTile <- level.colors(temp, pretty(temp, 200), colorRampPalette(recolor.map)(200))
        dim(map$myTile) <- ra[1:2]
    }

    #return map
    map
}



#################################
#################################
##getMapArg
#################################
#################################

getMapArg <- function(object = trellis.last.object()){
    #recovers map from previous GoogleMap plot
    object$panel.args.common$map
}






#################################
#################################
##panel handlers
#################################
#################################





#################################
#################################
##panel.GoogleMapsRaster
#################################
#################################
#raster map panel
#


panel.GoogleMapsRaster <- function(map){
    grid.raster(map$myTile,
         x = unit(map$xlim[1], "native"), y = unit(map$ylim[1], "native"),
         width = unit(map$xlim[2] - map$xlim[1], "native"),
         height = unit(map$ylim[2] - map$ylim[1], "native"),
         just = c("left", "bottom")
    )
}


#################################
#################################
##panel.GoogleMaps
#################################
#################################

#non-raster map panel
#

panel.GoogleMaps <- function(map){

    #there has to be a better way

    #both the rect handling
    #and the map.col generation
    #need thinking about

#############
#think about this
#############

#    if(attr(map$myTile, "type") != "local")
#        map <- localMapManager(map)

    ra <- dim(map$myTile)
    map.col <- map$myTile

    map.lon <- rep(seq(map$xlim[1], map$xlim[2],
                       length.out = ra[1]),
                   each = ra[2])
    map.lat <- rep(seq(map$ylim[2], map$ylim[1],
                       length.out = ra[2]),
                   time = ra[1])
    width <- (map$xlim[1] - map$xlim[2]) / ra[1]
    height <- (map$ylim[1] - map$ylim[1]) / ra[2]

    panel.rect(x = map.lon, y = map.lat,
               width = width, height = height,
               col = map.col, border = map.col)
}






###############################
###############################
##axis handling
###############################
###############################

##############################
#to do
##############################
#import pass-on for axis function
#any value in on-passing?
#



###############################
###############################
##xscale.components.GoogleMaps
###############################
###############################

xscale.components.GoogleMaps <- function(lim, ..., map = map){
        ans <- xscale.components.default(c(map$BBOX$ll[2], map$BBOX$ur[2]), ...)
        ans$num.limit <- map$xlim
        temp <- LatLon2XY.centered(map, map$BBOX$ll[1],
                                        ans$bottom$ticks$at)
        temp.2 <- LatLon2XY.centered(map, map$BBOX$ur[1],
                                          ans$bottom$ticks$at)
        ans$bottom$ticks$at <- temp$newX
        ans$bottom$labels$at <- temp$newX
        ans$top <- ans$bottom
        ans$top$ticks$at <- temp.2$newX
        ans$top$labels$at <- temp.2$newX
        ans
    }




###############################
###############################
##yscale.components.GoogleMaps
###############################
###############################

yscale.components.GoogleMaps <- function(lim, ..., map = map){
        ans <- yscale.components.default(c(map$BBOX$ll[1], map$BBOX$ur[1]), ...)
        ans$num.limit <- map$ylim

        temp <- LatLon2XY.centered(map, ans$left$ticks$at,
                                        map$BBOX$ll[2])
        temp.2 <- LatLon2XY.centered(map, ans$left$ticks$at,
                                          map$BBOX$ur[2])
        ans$left$ticks$at <- temp$newY
        ans$left$labels$at <- temp$newY
        ans$right <- ans$left
        ans$right$ticks$at <- temp.2$newY
        ans$right$labels$at <- temp.2$newY
        ans
    }






###############################
###############################
##axis.components.GoogleMaps
###############################
###############################

axis.components.GoogleMaps <- function(map, xlim = NULL, ylim = NULL, ...){

#the dots currently go no further!

    #default xlim, ylim to map size
    #if not suppied
    if(is.null(xlim))
        xlim <- map$xlim
    if(is.null(ylim))
        ylim <- map$ylim   
    #get and combine
    ans <- xscale.components.GoogleMaps(xlim, map = map)
    ans <- listUpdate(ans, yscale.components.GoogleMaps(ylim, map = map)) 
    ans
}















