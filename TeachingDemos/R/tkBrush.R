tkBrush <- function(mat,hscale=1.75,vscale=1.75,wait=TRUE,...){

  if( !requireNamespace('tkrplot', quietly = TRUE) ) stop('This function depends on the tkrplot package being available')

  first <- TRUE
  bp <- FALSE

  cols <- character(0)

  colhist <- function(x,...){
    tmp <- hist(x,plot=F)
    br <- tmp$breaks
    w <- as.numeric(cut(x,br,include.lowest=TRUE))
    sy <- unlist(lapply(tmp$counts,function(x)seq(length=x)))
    my <- max(sy)
    sy <- sy/my
    my <- 1/my
    sy <- sy[order(order(x))]
    tmp.usr <- par('usr'); on.exit(par(usr=tmp.usr))
    par(usr=c(tmp.usr[1:2],0,1.5))
    rect(br[w], sy-my, br[w+1], sy, col=cols, border=NA)
    rect(br[-length(br)], 0, br[-1], tmp$counts*my)
    if(first){
#      tmp <- cnvrt.coords((br[w]+br[w+1])/2,sy-my/2,'usr')$tdev
      tmp <- list(
                  x=grconvertX((br[w]+br[w+1])/2, to='ndc'),
                  y=grconvertY( sy-my/2, to='ndc')
                  )
      dx <<- c(dx,tmp$x)
      dy <<- c(dy,tmp$y)
      di <<- c(di,seq(along=tmp$x))
    }
  }

  pcols <- rep('black',nrow(mat))
  tcols <- rep(NA,nrow(mat))

  ppch <- rep(1,nrow(mat))
  tpch <- rep(NA,nrow(mat))

  dx <- dy <- di <- numeric(0)

  rx <- ry <- 0.5
  rw <- rh <- 0.05

  epch<-tcltk::tclVar(16)
  ecol<-tcltk::tclVar('red')

  devlims <- c(0.05,0.95,0.05,0.95)

  replot <- function(){
    if(first){
      cols <<- pcols
      pairs(mat, #upper.panel=NULL,
            panel=function(x,y,...){
              points(x,y,...)
#              tmp <- cnvrt.coords(x,y,'usr')$tdev
              tmp <- list(
                          x=grconvertX(x,to='ndc'),
                          y=grconvertY(y,to='ndc') )
              dx <<- c(dx,tmp$x)
              dy <<- c(dy,tmp$y)
              di <<- c(di,seq(tmp$x))
            },
            diag.panel=colhist)
      first <<- FALSE
    } else {

      cols <<- ifelse(is.na(tcols),pcols,tcols)

      pairs(mat, #upper.panel=NULL,
            diag.panel=colhist,
            pch=ifelse(is.na(tpch),ppch,tpch),
            col=ifelse(is.na(tcols),pcols,tcols))
      par(fig=c(0,1,0,1),plt=c(0,1,0,1),usr=c(0,1,0,1),xpd=TRUE)

      rect(rx-rw,ry,rx,ry+rh,border='green')
     }
  }

  tt <- tcltk::tktoplevel()
  tcltk::tkwm.title(tt,"Tk Brush")

  img <- tkrplot::tkrplot(tt, replot, vscale=vscale, hscale=hscale)

  tcltk::tkpack(img,side='left')


  tcltk::tkpack( tcltk::tklabel(tt,text='pch:'),side='top')
  tcltk::tkpack(tcltk::tkentry(tt,textvariable=epch),side='top')


  tcltk::tkpack( tcltk::tklabel(tt,text='Color:'),side='top')
  tcltk::tkpack( tcltk::tkentry(tt,textvariable=ecol),side='top')


  tcltk::tkpack( tcltk::tkbutton(tt, text='Quit', command=function()tcltk::tkdestroy(tt)),
         side='bottom')

  iw <- as.numeric(tcltk::tcl('image','width',tcltk::tkcget(img,'-image')))
  ih <- as.numeric(tcltk::tcl('image','height',tcltk::tkcget(img,'-image')))

  mm <- function(x,y){
    tx <- (as.numeric(x)-1)/iw
    ty <- 1-(as.numeric(y)-1)/ih

    if(tx-rw < devlims[1]) tx <- devlims[1]+rw
    if(tx > devlims[2]) tx <- devlims[2]
    if(ty < devlims[3]) ty <- devlims[3]
    if(ty+rh > devlims[4]) ty <- devlims[4] - rh

    rx <<- tx
    ry <<- ty

    tmp <- di[ dx >= rx-rw & dx <= rx & dy >= ry & dy <= ry+rh ]

    tmpc <- rep(NA,nrow(mat))
    tmpcol <- as.character(tcltk::tclvalue(ecol))
    if( !( tmpcol %in% colors() ) ) tmpcol <- 'black'
    tmpc[tmp] <- tmpcol
    tcols <<- tmpc

    tmpp <- rep(NA,nrow(mat))
    tmppch <-  as.numeric(tcltk::tclvalue(epch))
    if(is.na(tmppch)) tmppch <- as.character(tcltk::tclvalue(epch))
    tmpp[tmp] <- tmppch
    tpch <<- tmpp

    if(bp){
      ppch <<- ifelse(is.na(tpch),ppch,tpch)
      pcols <<- ifelse(is.na(tcols),pcols,tcols)
    }
    tkrplot::tkrreplot(img)
  }

  mmm <- function(){
    tmp <- di[ dx >= rx-rw & dx <= rx & dy >= ry & dy <= ry+rh ]

    tmpc <- rep(NA,nrow(mat))
    tmpcol <- as.character(tcltk::tclvalue(ecol))
    if( !( tmpcol %in% colors() ) ) tmpcol <- 'black'
    tmpc[tmp] <- tmpcol
    tcols <<- tmpc

    tmpp <- rep(NA,nrow(mat))
    tmppch <-  as.numeric(tcltk::tclvalue(epch))
    if(is.na(tmppch)) tmppch <- as.character(tcltk::tclvalue(epch))
    tmpp[tmp] <- tmppch
    tpch <<- tmpp

    if(bp){
      ppch <<- ifelse(is.na(tpch),ppch,tpch)
      pcols <<- ifelse(is.na(tcols),pcols,tcols)
    }
    tkrplot::tkrreplot(img)
  }

  tcltk::tkbind(img, '<Motion>', mm)
  tcltk::tkbind(img, '<ButtonPress-1>', function() {bp<<-TRUE;mmm()})
  tcltk::tkbind(img, '<ButtonRelease-1>', function() bp<<-FALSE)
  tcltk::tkbind(tt, '<Key-Up>',function(){rh <<- rh+0.01;mmm()})
  tcltk::tkbind(tt, '<Key-Down>',function(){rh <<- rh-0.01;mmm()})
  tcltk::tkbind(tt, '<Key-Left>',function(){rw <<- rw+0.01;mmm()})
  tcltk::tkbind(tt, '<Key-Right>',function(){rw <<- rw-0.01;mmm()})

  if(wait){
    tcltk::tkwait.window(tt)
    return(list(col=pcols, pch=ppch))
  } else {
    return(invisible(NULL))
  }
}

