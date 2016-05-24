`hardcopy` <-
function(width=3.75, height=3.75, color=FALSE, trellis=FALSE,
           device=c("","pdf","ps"), path=getwd(), file=NULL,
           format=c("nn-nn", "name"), split="\\.", pointsize=c(8,4),
           fonts=NULL,
           horiz=FALSE, ...){
    if(!trellis)pointsize <- pointsize[1]
    funtxt <- sys.call(1)
    nam <- strsplit(as.character(funtxt), "(", fixed=TRUE)[[1]][1]
    suffix <- switch(device, ps=".eps", pdf=".pdf")
    if(is.character(path) & nchar(path)>1 & substring(path, nchar(path))!="/")
      path <- paste(path, "/", sep="")
    if(is.null(file)) if(format[1]=="nn-nn"){
      if(!is.null(split))dotsplit <- strsplit(nam, split)[[1]] else
      dotsplit <- nam
      if(length(dotsplit)==1)dotsplit <- c("", dotsplit)
      nn2 <- paste(if(nchar(dotsplit[2])==1)"0" else "", dotsplit[2],
                   sep="")
      if(nchar(dotsplit[1])>0){
        numstart <- which(unlist(strsplit(dotsplit[1], "")) %in% paste(0:9))[1]
        nn1 <- substring(dotsplit[1], numstart)
        nn1 <- paste(if(nchar(nn1) == 1) "0" else "", nn1, "-", sep="")
      } else nn1 <- ""
      file <- paste(nn1, nn2, sep="")
    } else file <- nam
    if(nchar(file)>4 & substring(file, nchar(file)-nchar(suffix)+1)==suffix)
      suffix <- ""
    file <- paste(path, file, suffix, sep="")
    print(paste("Output will be directed to file:", file))
    dev.out <- device[1]
    dev.fun <- switch(dev.out, pdf=pdf, ps=postscript)
    if(trellis){
      if(device=="ps")
        trellis.device(file=file, device=dev.fun,
                       color = color, horiz=horiz, fonts=fonts,
                       width=width, height=height, ...) else
      trellis.device(file=file, device=dev.fun, fonts=fonts,
                     color = color, width=width, height=height, ...)

      trellis.par.set(list(fontsize=list(text=pointsize[1], points=pointsize[2])))
    } else
    if (dev.out!=""){
      print(c(width, height))
      if(device=="ps")
        dev.fun(file=file, paper="special",  horiz=horiz, fonts=fonts,
                width=width, height=height, pointsize=pointsize[1], ...) else
      dev.fun(file=file, paper="special", fonts=fonts,
              width=width, height=height, pointsize=pointsize[1], ...)
    }
    if(trellis)trellis.par.set(list(fontsize=list(text=pointsize[1],
                                      points=pointsize[2])))
  }

