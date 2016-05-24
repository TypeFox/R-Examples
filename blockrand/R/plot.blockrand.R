"plotblockrand" <-
function(x, file='blockrand.pdf', top, middle, bottom, blockrand.text,
                           width=11,height=8.5,par.args,
                           id.col='id',stratum.col='stratum',
                           treat.col='treatment', cut.marks=FALSE,
                           top.ho,top.vo,middle.ho,middle.vo,
                           bottom.ho,bottom.vo, nrow=2, ncol=2,
                           ...) {

  pdf(file,width=width,height=height,...)
  on.exit(dev.off())

  if(!missing(par.args)){
    par(par.args)
  }
  par(mfrow=c(nrow,ncol),mar=c(0,0,0,0),xaxs='i',yaxs='i')
  cur.par <- par()

  if(missing(top)){
    if(!missing(blockrand.text)) {
      top <- blockrand.text$top
    } else {
      top <- character(0)
    }
  }

  if(is.list(top)){
    t.text <- top$text
    t.col <- if(is.null(top$col)) cur.par$col else top$col
    t.fnt <- if(is.null(top$font)) cur.par$font else top$font
  } else {
    t.text <- top
    t.col <- cur.par$col
    t.fnt <- cur.par$font
  }


  if(missing(middle)){
    if(!missing(blockrand.text)) {
      middle <- blockrand.text$middle
    } else {
      middle <- character(0)
    }
  }

  if(is.list(middle)){
    m.text <- middle$text
    m.col <- if(is.null(middle$col)) cur.par$col else middle$col
    m.fnt <- if(is.null(middle$font)) cur.par$font else middle$font
  } else {
    m.text <- middle
    m.col <- cur.par$col
    m.fnt <- cur.par$font
  }


  if(missing(bottom)){
    if(!missing(blockrand.text)) {
      bottom <- blockrand.text$bottom
    } else {
      bottom <- character(0)
    }
  }

  if(is.list(bottom)){
    b.text <- bottom$text
    b.col <- if(is.null(bottom$col)) cur.par$col else bottom$col
    b.fnt <- if(is.null(bottom$font)) cur.par$font else bottom$font
  } else {
    b.text <- bottom
    b.col <- cur.par$col
    b.fnt <- cur.par$font
  }


  if(missing(top.ho)){
    if(!missing(blockrand.text) && length(blockrand.text$top.ho)){
      top.ho <- blockrand.text$top.ho
    } else {
      top.ho <- 0
    }
  }

  if(missing(top.vo)){
    if(!missing(blockrand.text) && length(blockrand.text$top.vo)){
      top.vo <- blockrand.text$top.vo
    } else {
      top.vo <- 0
    }
  }

  if(missing(middle.ho)){
    if(!missing(blockrand.text) && length(blockrand.text$middle.ho)){
      middle.ho <- blockrand.text$middle.ho
    } else {
      middle.ho <- 0
    }
  }

  if(missing(middle.vo)){
    if(!missing(blockrand.text) && length(blockrand.text$middle.vo)){
      middle.vo <- blockrand.text$middle.vo
    } else {
      middle.vo <- 0
    }
  }

  if(missing(bottom.ho)){
    if(!missing(blockrand.text) && length(blockrand.text$bottom.ho)){
      bottom.ho <- blockrand.text$bottom.ho
    } else {
      bottom.ho <- 0
    }
  }

  if(missing(bottom.vo)){
    if(!missing(blockrand.text) && length(blockrand.text$bottom.vo)){
      bottom.vo <- blockrand.text$bottom.vo
    } else {
      bottom.vo <- 0
    }
  }





  id <- x[[id.col]]
  if (is.null(id)) id <- seq(length.out=nrow(x))
  id <- as.character(id)

  stratum <- x[[stratum.col]]
  if (is.null(stratum)) stratum <- rep('',nrow(x))
  stratum <- as.character(stratum)

  treat <- x[[treat.col]]
  if (is.null(treat)) stop("Cannot proceed without valid treatment column")
  treat <- as.character(treat)

  plot.new()
  sh <- strheight("ABCjglw")*1.5*1.3
  sw <- strwidth("W")*1.3
  par(new=T)

  for (i in seq(along=id) ) {
    t.text2 <- gsub('%ID%',id[i],t.text)
    t.text2 <- gsub('%TREAT%',treat[i], t.text2)
    t.text2 <- gsub('%STRAT%',stratum[i], t.text2)

    m.text2 <- gsub('%ID%',id[i], m.text)
    m.text2 <- gsub('%STRAT%',stratum[i], m.text2)

    plot.new()
#axis(1,at=seq(0,1,.1),line=-10)
#axis(2,at=seq(0,1,.1),line=-10)

    text(0.15+sw*top.ho, 0.85+sh*top.vo - seq(0,length=length(t.text2))*sh,
         t.text2, col=t.col, font=t.fnt,adj=c(0,0.5),cex=1.3)

    text(0.5+sw*middle.ho,0.30+sh*middle.vo - seq(0,length=length(m.text2))*sh,
         m.text2, col=m.col, font=m.fnt, adj=c(0.5,0.5),cex=1.3)

    text(0.9+sw*bottom.ho,0.1+sh*bottom.vo, b.text, col=b.col, font=b.fnt,
         adj=c(1,0.5))

    if(cut.marks){
      lines(c(0,0,0.1),c(0.9,1,1),col='darkgrey')
      lines(c(0.9,1,1),c(1,1,0.9),col='darkgrey')
      lines(c(1,1,0.9),c(0.1,0,0),col='darkgrey')
      lines(c(0.1,0,0),c(0,0,0.1),col='darkgrey')
    }
  }



  return(invisible(NULL))
}

