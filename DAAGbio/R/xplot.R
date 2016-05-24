`xplot` <-
function(data=coralRG$R, images=1:6, layout=coralRG$printer, mfrow=c(3,2),
           FUN=imgplot, device=NULL, title=NULL, width=7.5, height=10,
           paneltitles=c("1:R/G","2:G/R", "3:R/G","4:G/R", "5:R/G","6:G/R"),
           paneltitles.line=0.5,
           mar=c(3.6,3.6,1.6,0.6), oma=c(0.6,0.6,1.6,0.6), file=NULL){
    if(is.null(title)){title <- as.character(substitute(data))
                       title <- paste(title[2], title[3], sep=":")
                     }
    if(is.null(file))file <- title
    nch <- nchar(title)
    if(!is.null(device)){devnam <- deparse(substitute(device))
                         ext <- switch(devnam, ps="ps", pdf="pdf", png="png",
                                       jpeg="jpg", bitmap="bmp")
                         file <- paste(title,".", ext, sep="")
                         print(file)
                         device(file=file, width=width, height=height)
                       }
    oldpar <- par(mfrow=mfrow, mgp=c(1,0.25,0), oma=oma, mar=mar)
    on.exit(par(oldpar))
    for(i in images){
      FUN(data[,i], layout=layout)
      mtext(side=3,line=paneltitles.line,paneltitles[i],adj=0)
    }
    mtext(side=3, line=0.25, title, outer=TRUE)
    if(!is.null(device))dev.off()
  }

