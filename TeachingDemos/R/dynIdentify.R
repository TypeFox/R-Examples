dynIdentify <- function(x,y,labels=seq_along(x),
                        corners = cbind( c(-1,0,1,-1,1,-1,0,1),
                                         c(1,1,1,0,0,-1,-1,-1) ),
                        ...) {

    lx <- x  # label positions
    ly <- y
    lx[ is.na(labels) ] <- NA
    ly[ is.na(labels) ] <- NA

    llx <- lx  # line end positions
    lly <- ly

    replot <- function() {
        plot(x,y,...)
        segments(x,y, llx,lly)
        text(lx,ly,labels)
    }

    replot()

#    tmp <- cnvrt.coords(x,y, input='usr')$dev
    dx <- grconvertX(x, to='ndc')
    dy <- grconvertY(y, to='ndc')  # device coordinates of points

    widths <- strwidth(labels)/2
    heights <- strheight(labels)/2

    ci <- 0  # current label

    mouse.down <- function(buttons, x, y){

        if( any(buttons==2) ){
            out <- list( labels=list(x=lx, y=ly),
                         lineends=list(x=llx, y=lly) )
            return(out)
        }
#        tmp <- cnvrt.coords(lx,ly, input='usr')$dev
        i <- which.min( (grconvertX(lx,to='ndc')-x)^2 + (grconvertY(ly,to='ndc')-y)^2 )
        ci <<- i
        NULL
    }

    mouse.up <- function(buttons, x, y){
#        tmp <- cnvrt.coords(x,y, input='dev')$usr
        cx <- grconvertX(x, from='ndc')
        cy <- grconvertY(y, from='ndc')

        tmpx <- cx + corners[,1]*widths[ci]
        tmpy <- cy + corners[,2]*heights[ci]
#        tmp <- cnvrt.coords(tmpx,tmpy, input='usr')$dev
        i <- which.min( (dx[ci] - grconvertX(tmpx,to='ndc'))^2 +
                       (dy[ci] - grconvertY(tmpy, to='ndc'))^2 )
#        tmp <- lx; tmp[ci] <- cx; lx <<- tmp
#        tmp <- ly; tmp[ci] <- cy; ly <<- tmp
#        tmp <- llx; tmp[ci] <- tmpx[i]; llx <<- tmp
#        tmp <- lly; tmp[ci] <- tmpy[i]; lly <<- tmp
        lx[ci] <<- cx
        ly[ci] <<- cy
        llx[ci] <<- tmpx[i]
        lly[ci] <<- tmpy[i]
        replot()

        NULL
    }

    out <- getGraphicsEvent( prompt= "Click on points and drag label to position.  \nRight click to exit\n",
                     onMouseDown=mouse.down,
                            onMouseUp=mouse.up)

    invisible(out)
}



TkIdentify <- function(x,y,labels=seq_along(x), hscale=1.75, vscale=1.75,
                       corners = cbind( c(-1,0,1,-1,1,-1,0,1),
                                        c(1,1,1,0,0,-1,-1,-1) ),
                       ...) {
    if( !requireNamespace('tkrplot', quietly = TRUE) ) stop('This function depends on the tkrplot package being available')

    md <- FALSE

    lx <- x  # label positions
    ly <- y
    lx[ is.na(labels) ] <- NA
    ly[ is.na(labels) ] <- NA

    llx <- lx  # line end positions
    lly <- ly

    first <- TRUE
    dx <- dy <- dlx <- dly <- widths <- heights <- numeric(0)
    ul <- ur <- ut <- ub <- 0

    replot <- function() {

        plot(x,y,...)
        segments(x,y, llx,lly)
        text(lx,ly,labels)
        if(first) {
            first <<- FALSE
#            tmp <- cnvrt.coords(x,y, input='usr')$dev
            dx <<- grconvertX(x, to='ndc')
            dy <<- grconvertY(y, to='ndc')
            widths <<- strwidth(labels)/2
            heights <<- strheight(labels)/2
#            tmp <- cnvrt.coords(c(0,1),c(0,1), input='dev')$usr
            ul <<- grconvertX(0, from='ndc')
            ur <<- grconvertX(1, from='ndc')
            ub <<- grconvertY(0, from='ndc')
            ut <<- grconvertY(1, from='ndc')
        }
#        tmp <- cnvrt.coords(lx,ly, input='usr')$dev
        dlx <<- grconvertX(lx, to='ndc')
        dly <<- grconvertY(ly, to='ndc')
    }

    tt <- tcltk::tktoplevel()
    tcltk::tkwm.title(tt, "TkIdentify")

    img <- tkrplot::tkrplot(tt, replot, vscale=vscale, hscale=hscale)
    tcltk::tkpack(img, side='top')

    tcltk::tkpack( tcltk::tkbutton(tt, text='Quit', command=function() tcltk::tkdestroy(tt)),
           side='right')

    corners <- cbind( c(-1,0,1,-1,1,-1,0,1), c(1,1,1,0,0,-1,-1,-1) )

    ci <- 0  # current label
    cx <- cy <- 0

    iw <- as.numeric(tcltk::tcl('image','width',tcltk::tkcget(img,'-image')))
    ih <- as.numeric(tcltk::tcl('image','height',tcltk::tkcget(img,'-image')))

    mouse.move <- function(x,y) {
        if(md){
            tx <- (as.numeric(x)-1)/iw
            ty <- 1-(as.numeric(y)-1)/ih
            cx <<- tx*ur + (1-tx)*ul
            cy <<- ty*ut + (1-ty)*ub

#            tmp <- lx; tmp[ci] <- cx; lx <<- tmp
#            tmp <- ly; tmp[ci] <- cy; ly <<- tmp
            lx[ci] <<- cx
            ly[ci] <<- cy

            tmpx <- cx + corners[,1]*widths[ci]
            tmpy <- cy + corners[,2]*heights[ci]

            tmpxx <- (tmpx - ul)/(ur-ul)
            tmpyy <- (tmpy - ub)/(ut-ub)

            i <- which.min( (dx[ci] - tmpxx)^2 +
                            (dy[ci] - tmpyy)^2 )

#            tmp <- lx; tmp[ci] <- cx; lx <<- tmp
#            tmp <- ly; tmp[ci] <- cy; ly <<- tmp
#            tmp <- llx; tmp[ci] <- tmpx[i]; llx <<- tmp
#            tmp <- lly; tmp[ci] <- tmpy[i]; lly <<- tmp
            lx[ci] <<- cx
            ly[ci] <<- cy
            llx[ci] <<- tmpx[i]
            lly[ci] <<- tmpy[i]

            tkrplot::tkrreplot(img)
        }
    }

    mouse.down <- function(x,y){
        tx <- (as.numeric(x)-1)/iw
        ty <- 1-(as.numeric(y)-1)/ih

        ci <<- which.min( (tx - dlx)^2 + (ty - dly)^2 )
        md <<- TRUE
        mouse.move(x,y)
    }

    mouse.up <- function(x,y){
        md <<- FALSE
    }

    tcltk::tkbind(img, '<Motion>', mouse.move)
    tcltk::tkbind(img, '<ButtonPress-1>', mouse.down)
    tcltk::tkbind(img, '<ButtonRelease-1>', mouse.up)

    tcltk::tkwait.window(tt)
    out <- list( labels=list(x=lx, y=ly),
                 lineends=list(x=llx, y=lly) )
    invisible(out)
}

### old version, possibilities for the Tk version
##   dynIdentify <- function(x,y,labels=seq_along(x), ...) {
##       plot(x,y,...)
##
##       tmp <- cnvrt.coords(x,y, input='usr')$dev
##       dx <- tmp$x
##       dy <- tmp$y  # device coordinates of points
##
##   print(dx)
##   print(dy)
##
##       lx <- rep(NA, length(x) )  # label positions
##       ly <- rep(NA, length(y) )
##       llx <- lx  # line end positions
##       lly <- ly
##
##       widths <- strwidth(labels)/2
##       heights <- strheight(labels)/2
##
##       corners <- cbind( c(-1,0,1,-1,1,-1,0,1), c(1,1,1,0,0,-1,-1,-1) )
##
##       md <- FALSE # mouse button down
##
##       cx <- 0 # current
##       cy <- 0
##       ci <- 0
##
##       replot <- function() {
##           if(!md){return()}
##           plot(x,y,...)
##           segments(x,y, llx,lly)
##           text(lx,ly,labels)
##       }
##
##       mouse.move <- function(buttons, x, y){
##
##           tmp <- cnvrt.coords(x,y, input='dev')$usr
##           cx <<- tmp$x
##           cy <<- tmp$y
##           if(md){
##
##               tmpx <- cx + corners[,1]*widths[ci]
##               tmpy <- cy + corners[,2]*heights[ci]
##               tmp <- cnvrt.coords(tmpx,tmpy, input='usr')$dev
##               i <- which.min( (dx[ci] - tmp$x)^2 +
##                               (dy[ci] - tmp$y)^2 )
##               tmp <- lx; tmp[ci] <- cx; lx <<- tmp
##               tmp <- ly; tmp[ci] <- cy; ly <<- tmp
##               tmp <- llx; tmp[ci] <- tmpx[i]; llx <<- tmp
##               tmp <- lly; tmp[ci] <- tmpy[i]; lly <<- tmp
##               replot()
##           }
##           NULL
##       }
##
##       mouse.down <- function(buttons, x, y){
##
##           if( any(buttons==2) ){
##               out <- list( labels=list(x=lx, y=ly),
##                            lineends=list(x=llx, y=lly) )
##               return(out)
##           }
##           i <- which.min( (dx-x)^2 + (dy-y)^2 )
##           ci <<- i
##           md <<- TRUE
##           mouse.move(buttons, x, y)
##           NULL
##       }
##
##       mouse.up <- function(buttons, x, y){
##
##           tmp <- dx; tmp[ci] <- NA; dx <<- tmp
##           tmp <- dy; tmp[ci] <- NA; dy <<- tmp
##
##           if(all(is.na(dx))) {
##               out <- list( labels=list(x=lx, y=ly),
##                            lineends=list(x=llx, y=lly) )
##               return(out)
##           }
##
##           md <<- FALSE
##           ci <<- 0
##           NULL
##       }
##
##       out <- getGraphicsEvent( prompt= "Click on points and drag label to position.\nRight click to exit\n",
##                        onMouseDown=mouse.down,
##                        onMouseMove=mouse.move,
##                        onMouseUp=mouse.up)
##
##       invisible(out)
##   }
