## Emilio Torres Manzanera
## University of Oviedo
## Time-stamp: <2016-01-13 18:08 emilio on emilio-despacho>
## ============================================================

##' Plot the axis
##' 
##' This function plots the axis
##'
##' It plots the axis of the graph. 
##'
##' @param xrange The range of the X axe.
##' @param yrange The range of the Y axe.
##' @param ... Other arguments.
##' @return A layer with the axis. 
##' @import ggplot2
##' @import extrafont
##' @importFrom Hmisc bezier
##' @importFrom stats runif
##' @export
##' @examples
##' xrange <- range(mtcars$mpg)
##' yrange <- range(mtcars$wt)
##' p <- ggplot() +
##'      geom_point(aes(mpg, wt), data=mtcars) +
##'      xkcdaxis(xrange,yrange)
##' p
xkcdaxis <- function(xrange, yrange, ...) {
  if( is.null(xrange) | is.null(yrange) )
    stop("Arguments are: xrange, yrange")
  xjitteramount <- diff(xrange)/50
  yjitteramount <- diff(yrange)/50
  ## This cause R CMD check to give the note
  ## 'no visible binding for global variable'
  ## Notes do not forbbiden the submission
  ## I will follow this suggestion:
  ## http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
  ##mappingsegment <- aes(x=x,y=y,xend=xend,yend=yend) ## Put it within a with!!!
 
  dataaxex <- data.frame(xbegin=xrange[1]-xjitteramount,
                          ybegin=yrange[1]-yjitteramount,
                          xend=xrange[2]+xjitteramount,
                          yend=yrange[1]-yjitteramount)
  mappingsegment <- with(dataaxex, aes(xbegin=xbegin,ybegin=ybegin,xend=xend,yend=yend))
  axex <- xkcdline(mappingsegment, dataaxex, yjitteramount = yjitteramount, mask = FALSE, ... )
  
 
  dataaxey <- data.frame(xbegin=xrange[1]-xjitteramount,
                          ybegin=yrange[1]-yjitteramount,
                          xend=xrange[1]-xjitteramount,
                          yend=yrange[2]+yjitteramount)
  mappingsegment <- with(dataaxey, aes(xbegin=xbegin,ybegin=ybegin,xend=xend,yend=yend))
  axey <- xkcdline(mappingsegment, dataaxey, xjitteramount = xjitteramount, mask = FALSE, ... )
  coordcarte <- coord_cartesian(xlim = xrange + 1.5*c(-xjitteramount,xjitteramount),
                                ylim = yrange + 1.5*c(-yjitteramount,yjitteramount))
  list(c(axex,axey), coordcarte,theme_xkcd())
}




theme_xkcd <- function(){
  if( "xkcd" %in% extrafont::fonts() ) {
    theme(panel.grid.major = element_blank(),
          ##axis.ticks = element_blank(),
          axis.ticks = element_line(colour = "black"),
          panel.background = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key = element_blank(),
          strip.background = element_blank(),
          text = element_text(size = 16, family = "xkcd"))
  } else {
    warning("Not xkcd fonts installed! See vignette(\"xkcd-intro\")")
    theme(panel.grid.major = element_blank(),
          ##axis.ticks = element_blank(),
          axis.ticks = element_line(colour = "black"),
          panel.background = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key = element_blank(),
          strip.background = element_blank(),
          text = element_text(size = 16))} 
}


## Therefore, we get a data frame with the default names of the mapping
## and a new mapping with the names by default
## For instance
## mapping <- aes(x= x1 +y1, y = y1) -> mapping <- aes(x= x, y = y)
## data[ , c("x1","y1","color")]  -> data[, c("x","y","x1","y1","color")]
createdefaultmappinganddata <- function(mapping, data, mandatoryarguments =c("x","y")) {
  
  ## Check the names of the aes
  nm <- names(mapping)
  positionswithoutname <- (1:length(nm))[nm==""]
  failsthisarguments <- mandatoryarguments[ !(mandatoryarguments %in% nm) ]
  if(length(failsthisarguments) != length(positionswithoutname))
    stop(paste("Argumenst of aes are ", paste(mandatoryarguments, collapse=", "),".",sep=""))
  names(mapping)[positionswithoutname] <- failsthisarguments
  ## New names
    namesmapping <- names(mapping)

  ## Create a new data
  ## For each name of the mapping, evaluate it and create a new data base
  ## with the names of the mapping.
  dataaes <- as.data.frame(lapply(mapping, function(xnamedataxkcdveryrare.327) with(data, eval(xnamedataxkcdveryrare.327))))
  ## Add the rest of variables of the data base

  variablestocbind <- names(data)[!(names(data) %in% namesmapping)] 
  dataaes[, variablestocbind] <- data[,variablestocbind]
  ## Now, it creates a new mapping with the default variables x=x, y=x, yend=yend, and so on.
  ## See the definition of the function ggplot2::aes_string
  parsed <- lapply(namesmapping, function(x) parse(text = x)[[1]])
  names(parsed) <- namesmapping
  newmapping <- structure(parsed, class = "uneval")
  list(mapping = newmapping, data = dataaes)
  }


## Apply a FUN to each row of the DATA
## If the arguments of the FUN are in the DATA and in the ELLIPSIS
## only use the variable of the DATA
##
## If doitalsoforoptargs = TRUE, then try to get the row of the ELLIPSIS variable
## when applying the function FUN to each row of the DATA.
## Otherwise, use the original ELLIPSIS variable
## when calling the function FUN
doforeachrow <- function(data, fun, doitalsoforoptargs, ...) {
  ## Do not pass the variables of the ELLIPSIS
  ## that are they are in the DATA
  argList <- list(...)
  for( i in intersect(names(data), names(argList) ) )
    argList[i] <- NULL
  if(doitalsoforoptargs) {
    ## If there are variable of ELLIPSIS with the same length than the data base
    ## copy them to the data base and delete them from the ELLIPSIS
    for( i in  names(argList) ) {
      if(!(is.null(argList[i])==TRUE)){
        if(length(argList[[i]]) == 1
           | length(argList[[i]]) == dim(data)[1] ) {
          data[,i] <-  argList[[i]] 
          argList[i] <- NULL
        }
      }
    }
  }
  ##print(data)
  ## Now, apply for each row the FUN
  lapply(1:(dim(data)[1]),
         function(i, data, fun,  argList) {
           largstopass <- as.list(data[i,,])
           mylistofargs <- c(largstopass, argList)
           ## Arguments of the function?
           fcn <- get(fun, mode = "function")
           argsfcntt <-  names(formals(fcn))
           if( "..." %in% argsfcntt ) do.call(fun, mylistofargs)
           else { ## we can only pass the arguments of the function
             for( i in names(mylistofargs)[ !(names(mylistofargs) %in% argsfcntt )])
               mylistofargs[i] <- NULL
             do.call(fun, mylistofargs)
           }
         },
         data = data,
         fun = fun,
         argList = unlist(argList)
         )
}

mappingjoin <- function(x,y) {
  nm1 <- names(x)
  nm2 <- names(y)
  for( i in intersect(nm1,nm2)) y[[i]] <- NULL
  parsed <- lapply(c(x,y), function(x) parse(text = x)[[1]])
  structure(parsed, class = "uneval")
}



pointscircunference <- function(x =0, y=0, diameter = 1, ratioxy=1, npoints = 16, alpha=  runif(1, 0, pi/2)){
    ##require(Hmisc) # bezier
    center <- c(x,y)
    r <- rep( diameter / 2, npoints )
    tt <- seq(alpha,2*pi + alpha,length.out = npoints)
    r <- jitter(r)
    sector <-  tt > alpha & tt <= ( pi/ 2 + alpha)
    r[ sector ] <- r[sector] * 1.05
    sector <-  tt > ( 2 * pi/2 + alpha)  & tt < (3* pi/ 2 +alpha)
    r[ sector ] <- r[sector] * 0.95    
    xx <- center[1] + r * cos(tt) * ratioxy
    yy <- center[2] + r * sin(tt) 
    ##return(data.frame(x = xx, y = yy))
    return(data.frame(bezier(x = xx, y =yy,evaluation=60)))
}


pointssegment <- function(xbegin, ybegin, xend, yend, npoints = 10, xjitteramount= 0, yjitteramount=0, bezier = TRUE) {
  ##require(Hmisc) # bezier
  if(npoints < 2 )
    stop("npoints must be greater than 1")
  ## If there are no jitters, do not interpolate
  if( xjitteramount == 0 & yjitteramount == 0) npoints <- 2 
  x <- seq(xbegin,xend,length.out = npoints)
  if( (xend - xbegin) != 0 ) {
    y <- (yend - ybegin) * ( x - xbegin ) / (xend - xbegin) + ybegin
  } else {
    y <-  seq(ybegin, yend, length.out = npoints)
  }
  if(xjitteramount !=0) x <- jitter(x, amount=xjitteramount)
  if(yjitteramount !=0) y <- jitter(y, amount=yjitteramount)
  x[1] <- xbegin
  y[1] <- ybegin
  x[length(x)] <- xend
  y[length(y)] <- yend
  if(bezier & length(x)>2 & (xjitteramount != 0 | yjitteramount != 0)) {
    data <- data.frame(bezier(x=x, y=y, evaluation=30))
  }
  else data <- data.frame(x=x,y=y)   
  data
}
