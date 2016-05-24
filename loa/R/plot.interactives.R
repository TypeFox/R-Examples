#in development code
#[TBC - NUMBER] functions 

#getXY
#getLatLon

#screenLatticePlot

#NOTE: much borrowed from lattice 


#####################
#####################
##getXY
#####################
#####################

#locate points on plot
#with correction for 
#[in development]

#to do

#think about getXY argument order
#think about getLatLon order
#might want to hide unit 
#drop plot.list control?

#overplotting issue with using alpha
#see e.g. getXY(n=5, col="black", cex=10, alpha=0.5)


getXY <- function (n = -1, ..., unit = "native", scale.correction = NULL) 
{

    #check
    if(n == 0)
        stop("'n = 0' ends call without action. (See ?getXY)", call. = FALSE)

    #setup
    extra.args <- list(...)
    ans <- NULL
    n <- n - 1

    #get safe trellis focus
    focus.list <- listUpdate(extra.args, list(), 
                             use.a = names(formals(trellis.focus)))   
    do.call(trellis.focus, focus.list)

    #check for plot arguments
    plot.list <- listUpdate(extra.args, list(), 
                             use.a = c("type", "pch", "lty", "col", 
                                          "cex", "lwd", "font", "fontfamily", 
                                          "fontface", "col.line", "col.symbol", 
                                          "alpha", "fill"))

#####################
#this loop could be tidier
#something better than while?
#####################

    #loop grid.locator 
    while (!is.null(temp <- grid.locator(unit = unit)) & n != 0) {
        n <- n - 1
        ans$x <- c(ans$x, as.numeric(temp$x))
        ans$y <- c(ans$y, as.numeric(temp$y))
        if(length(plot.list)>0)
            lpoints(ans$x, ans$y, ...)
    }
    #check last case
    if(!is.null(temp)){
        ans$x <- c(ans$x, as.numeric(temp$x))
        ans$y <- c(ans$y, as.numeric(temp$y))
        if(length(plot.list)>0)
            lpoints(ans$x, ans$y, ...)
    }

    #rescale if required

############
#could be make this 
#function only
############

    if (!is.null(scale.correction)) 
        ans <- scale.correction(ans)

    #unfocus and return results
    trellis.unfocus()
    ans
}



#####################
#####################
##getLatLon
#####################
#####################

getLatLon <- function(..., map = NULL, object = trellis.last.object(),  
            scale.correction = function(x) {
                                temp <- XY2LatLon(map, x$x, x$y)
                                as.list(as.data.frame(temp))
                            })
{

    if(is.null(map))
        map <- getMapArg(object)

#need new error catcher

#    if(missing(map))
#        stop("need 'map' as reference. (See ?getLatLon)", 
#             call. = FALSE)

    #pass to getXY 
    getXY(..., scale.correction = scale.correction)
}




######################
######################
##screenLatticePlot
######################
######################

#needs a lot of tidying but it is promising


screenLatticePlot <- function(object = trellis.last.object(),...){
#temp

this.panel <- object$panel
xlim <- object$x.limits
ylim <- object$y.limits

x.ref <- (max(xlim) - min(xlim))/20
y.ref <- (max(ylim) - min(ylim))/20

x1 <- (x.ref *18) + min(xlim)
y1 <- (y.ref *18) + min(ylim) 

x2 <- (x.ref *15) + min(xlim)
y2 <- (y.ref *18) + min(ylim) 

x3 <- (x.ref *12) + min(xlim)
y3 <- (y.ref *18) + min(ylim) 

x4 <- (x.ref *9) + min(xlim)
y4 <- (y.ref *18) + min(ylim) 


my.x <- c(NA, NA)
my.y <- c(NA, NA)



test <- FALSE

while(!test){

    test2 <- FALSE

    temp <- object$x.limits
    x.ref <- (max(temp) - min(temp))/20
    x1 <- (x.ref *18) + min(temp)
    x2 <- (x.ref *15) + min(temp)
    x3 <- (x.ref *12) + min(temp)
    x4 <- (x.ref *9) + min(temp)

    temp <- object$y.limits
    y.ref <- (max(temp) - min(temp))/20
    y1 <- (y.ref *18) + min(temp)
    y2 <- y1
    y3 <- y1
    y4 <- y1

    polygon <- list(x=c(x.ref, x.ref, -x.ref, -x.ref),
                    y=c(y.ref, -y.ref, -y.ref, y.ref))

    control.panel <- function(...){
        temp <- listUpdate(list(...), list(x=x1, y=y1, polygon=polygon, 
                                       loa.scale=list(fit="absolute")))
        do.call(loaPolygon, temp)
        llines(x=c(x1+x.ref, x1-x.ref), y=c(y1+y.ref, y1-y.ref), col="green")

        temp <- listUpdate(list(...), list(x=x2, y=y2, polygon=polygon, 
                                           loa.scale=list(fit="absolute")))
        do.call(loaPolygon, temp)
        temp <- listUpdate(list(...), list(x=x3, y=y3, polygon=polygon, 
                                           loa.scale=list(fit="absolute")))
        do.call(loaPolygon, temp)

        temp <- listUpdate(list(...), list(x=x4, y=y4, polygon=polygon, 
                                           loa.scale=list(fit="absolute")))
        do.call(loaPolygon, temp)

    }

    reaction.panel <- function(x=x, y=y, ...){
         x <- na.omit(my.x)
         y <- na.omit(my.y)

         if(length(x)==1 && length(y)==1)
             lpoints(x=x,y=y, cex=4, pch=3, col="black")

         if(length(x)>1 && length(y)>1){
             x0 <- c(max(x), max(x), min(x), min(x))
             y0 <- c(max(y), min(y), min(y), max(y))
             lpolygon(x=x0, y=y0)
         } 

    } 

    object$panel <- function(...){
        this.panel(...)
        control.panel(...)
        reaction.panel(...)
    }

    plot(object)

    ans <- getXY(n=1)

    if(!test){
        if(length(ans$x)<1 || length(ans$y)<1){
            output <- "cancelled"
            test <- TRUE
        }
    }
    if(!test){
         if(ans$x > x1-x.ref & ans$x < x1+x.ref & ans$y > y1-y.ref & ans$y < y1+y.ref){
                 test <- TRUE
                 output <- "button" 
         }
    }

    if(!test){
         if(ans$x > x2-x.ref & ans$x < x2+x.ref & ans$y > y2-y.ref & ans$y < y2+y.ref){
                 test <- TRUE
                 output <- "button2" 
         }
    }


    if(!test){
         if(ans$x > x3-x.ref & ans$x < x3+x.ref & ans$y > y3-y.ref & ans$y < y3+y.ref){
                 #test <- TRUE
                 #output <- "button3"
                 #zoom
                if(length(na.omit(my.x))>1 && length(na.omit(my.y))>1){
                    x0 <- c(max(my.x), max(my.x), min(my.x), min(my.x))
                    y0 <- c(max(my.y), min(my.y), min(my.y), max(my.y))
                    object$x.limits <- x0
                    object$y.limits <- y0
                    object$panel.args.common$xlim <- x0
                    object$panel.args.common$ylim <- x0
                    my.x <- c(NA,NA)
                    my.y <- c(NA,NA)
                    test2 <- TRUE
                }  
         }
    }

    if(!test){
         if(ans$x > x4-x.ref & ans$x < x4+x.ref & ans$y > y4-y.ref & ans$y < y4+y.ref){
                 #test <- TRUE
                 #output <- "button4"
                 #reset
                 object$x.limits <- xlim
                 object$y.limits <- ylim
                 object$panel.args.common$xlim <- xlim
                 object$panel.args.common$ylim <- ylim
                 my.x <- c(NA,NA)
                 my.y <- c(NA,NA)
                 test2 <- TRUE             
         }
    }

    if(!test & !test2){
       my.x <- c(ans$x[1], my.x)[1:2]
       my.y <- c(ans$y[1], my.y)[1:2]
    }



}

print(output)

#end
object$panel <- this.panel
object

}


