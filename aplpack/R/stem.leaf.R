## ms <-
stem.leaf <- function(data, unit, m, Min, Max, 
     rule.line=c("Dixon", "Velleman", "Sturges"),
     style=c("Tukey", "bare"), trim.outliers=TRUE, depths=TRUE,
     reverse.negative.leaves=TRUE,na.rm=FALSE,printresult=TRUE){
  if(missing(data)){cat("Author:  Peter Wolf 05/2003,", 
                        "(modified slightly by J. Fox, 20 July 03)",
                        "03/2006 additional rounding to prevent misclasification",
                        "07/2008 counting of NA's, 04/2009 improvement of rounding",
                        "syntax: stem.leaf(data.set)\n",sep="\n")
                    return("Warning: no data set found by stem.leaf")
  }                  
  rule.line <- match.arg(rule.line)
  style <- match.arg(style)
  n.na  <- sum(is.na(data))
  if(0<n.na){
    data <- data[!is.na(data)]
    if(na.rm){ # data<-data[!is.na(data)]
      print("Warning: NA elements have been removed!!")
    }else{
  #   data[is.na(data)] <- mean(data,na.rm=TRUE)
  #   print("Warning: NA elements have been exchanged by the mean value!!")
    }  
  }


  if(0<length(h<-find("debug.cond")) && ".GlobalEnv" %in% h){
    debug.cond<-get("debug.cond",envir=.GlobalEnv)
  } else debug.cond<-""
  debug.show<-function(name){
    if(!exists("debug.cond")) return()
    if(debug.cond=="all"|| (name %in% debug.cond) ){
      cat(name,":\n"); obj<-eval(parse(text=name))
      if(is.vector(obj)){ print(obj) }
      return()
    }
  }

  ##################################################################
  #Description:                                                    #
  #   stem.leaf  produces a stem-and-leaf-display of a data set    #
  #                                                                #
  #Usage:                                                          #
  #   stem.leaf(data)                                              #
  #   stem.leaf(data,unit=100,m=5,Min=50,Max=1000,                 #
  #     rule.line=c("Dixon", "Velleman", "Sturges"),               #
  #     style=c("Tukey", "bare"), trim.outliers=TRUE, depths=TRUE, #
  #     reverse.negative.leaves=TRUE,na.rm=FALSE)                  #
  #                                                                #
  #Arguments:                                                      #
  #   data:      vector of input data                              #
  #   unit:      unit of leaves in: { ...,100,10,1,.1,.01,... }    #
  #   m:         1, 2 or 5 -- 10/m=number of possible leaf digits  #
  #   Min:       minimum of stem                                   #
  #   Max:       maximum of stem                                   #
  #   rule.line:   = "Dixon"    => number of lines <- 10*log(n,10) #
  #                = "Velleman" => number of lines <- 2*sqrt(n)    #
  #                = "Sturges"  => number of lines <- 1 + log(n,2) #
  #   style:       = "Tukey"    => Tukey-like stem ( m = 2, 5 )    #
  #   trim.outliers=TRUE        => outliers are printed absent     #
  #   depths       =TRUE        => depths info is printed          #
  #   reverse.negative.leaves=TRUE => neg.leaves are rev. sorted   #
  #Author:                                                         #
  #   Peter Wolf 05/2003 (modified slightly by J. Fox, 20 July 03) #
  #   rounding operation for comparing added 29 March 06           #
  #   07/2008 NA-values are counted if na.rm==FALSE                #
  ##################################################################

  n <- length(data <- sort(data))
  row.max <- floor(  c(Dixon   =10*log(n,10),
                       Velleman=2*sqrt(n),
                       Sturges =1+log(n,2)        ))[rule.line]

  stats <- boxplot(data,plot=FALSE)
  if(missing(Min)) Min <- if (trim.outliers) stats$stats[1,1] else min(data, na.rm=TRUE)
  if(missing(Max)) Max <- if (trim.outliers) stats$stats[5,1] else max(data, na.rm=TRUE)
  spannweite.red<-Max - Min

  zeilen.intervall.laenge <- spannweite.red / row.max
  if(missing(unit)){
         factor <- 10^ceiling(log(zeilen.intervall.laenge,10))
  } else factor <- 10^round(log(unit*10,10))
  debug.show("factor")

  z <- zeilen.intervall.laenge/factor  # z in (0.1 ,1]
  delta.tick <- c(.2,.2,.5,1)[sum(z>c(0,.1,.2,.5))]

  if(missing(m)) m <- round(1/delta.tick) else delta.tick <- 1/m
  debug.show("delta.tick"); debug.show("m")

  data.tr <- data/factor
  Min.tr  <- Min/factor
  Max.tr  <- Max/factor

  spannweite.red <- Max.tr - Min.tr
  sk.min <- floor(Min.tr)
  sk.max <- ceiling(Max.tr)
  skala  <- seq(sk.min,sk.max,by=delta.tick)
  if(sk.min<0) skala <- c(sk.min-delta.tick,skala)
  if(sk.max<0) skala <- skala[-length(skala)]
  debug.show("skala")


  lo.limit <- if (trim.outliers) skala[1] else -Inf
  lo.log   <- if(skala[1   ] <  0) data.tr <= lo.limit else data.tr <  lo.limit
  n.sk <- length(skala)
  hi.limit <- if (trim.outliers) skala[n.sk] + delta.tick else Inf
  hi.log   <- if(skala[n.sk] >= 0) data.tr >= hi.limit else data.tr >  hi.limit

  n.lower.extr.values <- sum(lo.log); n.upper.extr.values <- sum(hi.log)
  if(0<n.lower.extr.values){
    lower.line <- paste("LO:", paste(data[lo.log],collapse=" "))
  }
  if(0<n.upper.extr.values){
    upper.line <- paste("HI:", paste(data[hi.log],collapse=" "))
  }
  data.tr.red <-data.tr[(!lo.log)&(!hi.log)]


  stem  <- ifelse(data.tr.red<0, ceiling(data.tr.red), floor(data.tr.red) )
  # eps <- 1e-12; leaf <- floor(abs(data.tr.red*10-stem*10)+eps)
  leaf  <- floor(10*abs(signif(data.tr.red-stem,10)))
  debug.show("leaf"); debug.show("stem")

  class.of.data.tr<-unlist(c(
     sapply(signif(data.tr.red[data.tr.red< 0],10),
       function(x,sk)length(sk)-sum(-sk<=-x),signif(skala,10))
    ,sapply(signif(data.tr.red[data.tr.red>=0],10),
       function(x,sk)sum( sk<= x),signif(skala,10))
  ))
  debug.show("class.of.data.tr")
  class.of.data.tr  <- c(1:length(skala),class.of.data.tr)
  leaf.grouped      <- split(c(rep(-1,length(skala)),leaf),class.of.data.tr)
  leaf.grouped      <- lapply(leaf.grouped, function(x){ sort(x[-1]) })
  # debug.show("leaf.grouped")

  class.negative <- skala < 0
  class.neg.zero <- floor(skala) == -1

  if (reverse.negative.leaves){
          for (i in seq(class.negative))
              if (class.negative[i]) leaf.grouped[[i]] <- rev(leaf.grouped[[i]])
  }

  leaf.grouped.ch <- paste("|",unlist(lapply(leaf.grouped,paste,collapse="")))
  # debug.show("leaf.grouped")

  line.names <- skala
  line.names[class.negative] <- line.names[class.negative]+1
  line.names <- as.character(floor(line.names))
  line.names[class.neg.zero] <- "-0"


  if(style=="Tukey"){
    switch(as.character(m),
    "1"={},
    "2"={
          h<-round(2*(skala%%1)) #; line.names[h!=0] <- ""
          line.names<-paste(line.names,
                  ifelse(skala<0,c(".","*")[1+h],c("*",".")[1+h]),sep="")
        },
    "5"={
          h<-round(5*(skala%%1)); line.names[h>0 & h<4] <- ""
          line.names<-paste(line.names, ifelse(skala<0,
                           c(".","s","f","t","*")[1+h],
                           c("*","t","f","s",".")[1+h]), sep="")
        }
    )
  }
  ragged.left <- function(ch.lines){
    max.n <- max(n.lines<-nchar(ch.lines))
    h     <- paste(rep(" ",max.n),collapse="")
    ch.lines <- paste( substring(h,1,1+max.n-n.lines), ch.lines)
    ch.lines
  }

  line.names <- ragged.left(line.names)


  n.class <- unlist(lapply(leaf.grouped,length))
  select <- (cumsum(n.class) > 0) & rev((cumsum(rev(n.class)) > 0))
  depth     <-     cumsum(n.class)          + n.lower.extr.values
  depth.rev <- rev(cumsum(rev(n.class))     + n.upper.extr.values)
  debug.show("depth")

  uplow <- depth>=depth.rev
  pos.median <- which(uplow)[1] + (-1:0)
  h <- abs(depth[pos.median]-depth.rev[pos.median])
  pos.median <- pos.median[1]+(h[1]>h[2])
  debug.show("pos.median")

  depth[uplow] <- depth.rev[uplow]
  depth <- paste(depth,"")
  depth[pos.median] <- paste("(",n.class[pos.median],")",sep="")
  depth[n.class==0] <- " "
  depth <- if (depths) ragged.left(depth) else ""


  info <- c(  paste("1 | 2: represents",1.2*factor),
           #  paste("    m:",m     ),
              paste(" leaf unit:",factor/10),
              paste("            n:",n     ))


  stem <- paste(depth, line.names, leaf.grouped.ch)
  stem <- if((m!=5)||sum(select)>4) stem[select] else stem
  result <- list(display=stem)
  if(exists("lower.line")) result<-c(lower=lower.line,result)
  if(exists("upper.line")) result<-c(result,upper=upper.line)
  if(0<n.na&&!na.rm) result<-c(result,NAs=paste("NA's:",n.na,collapse=" "))
  result <- c(list( info=info), result)
  if(printresult){ for(i in seq(result)) cat(result[[i]],sep="\n") }
  result <- c(result, list(depths=depth, stem=line.names, leaves=leaf.grouped.ch))
  invisible(result)

}

stem.leaf.backback <- function(x,y, unit, m, Min, Max, rule.line = c("Dixon", "Velleman", 
    "Sturges"), style = c("Tukey", "bare"), trim.outliers = TRUE, 
    depths = TRUE, reverse.negative.leaves = TRUE, na.rm = FALSE,
    printresult=TRUE, show.no.depths = FALSE, add.more.blanks = 0,
    back.to.back = TRUE){
  x.name <- paste(deparse(substitute(x)),collapse="")
  if(missing(y)){ y <- x; y.name <- x.name } 
  else y.name <- paste(deparse(substitute(y)),collapse="")
  n.na.x <- sum(is.na(x)); n.na.y <- sum(is.na(y))

  ## ms <-
  stem.leaf <- function(data, unit, m, Min, Max, 
       rule.line=c("Dixon", "Velleman", "Sturges"),
       style=c("Tukey", "bare"), trim.outliers=TRUE, depths=TRUE,
       reverse.negative.leaves=TRUE,na.rm=FALSE,printresult=TRUE){
    if(missing(data)){cat("Author:  Peter Wolf 05/2003,", 
                          "(modified slightly by J. Fox, 20 July 03)",
                          "03/2006 additional rounding to prevent misclasification",
                          "07/2008 counting of NA's, 04/2009 improvement of rounding",
                          "syntax: stem.leaf(data.set)\n",sep="\n")
                      return("Warning: no data set found by stem.leaf")
    }                  
    rule.line <- match.arg(rule.line)
    style <- match.arg(style)
    n.na  <- sum(is.na(data))
    if(0<n.na){
      data <- data[!is.na(data)]
      if(na.rm){ # data<-data[!is.na(data)]
        print("Warning: NA elements have been removed!!")
      }else{
    #   data[is.na(data)] <- mean(data,na.rm=TRUE)
    #   print("Warning: NA elements have been exchanged by the mean value!!")
      }  
    }


    if(0<length(h<-find("debug.cond")) && ".GlobalEnv" %in% h){
      debug.cond<-get("debug.cond",envir=.GlobalEnv)
    } else debug.cond<-""
    debug.show<-function(name){
      if(!exists("debug.cond")) return()
      if(debug.cond=="all"|| (name %in% debug.cond) ){
        cat(name,":\n"); obj<-eval(parse(text=name))
        if(is.vector(obj)){ print(obj) }
        return()
      }
    }

    ##################################################################
    #Description:                                                    #
    #   stem.leaf  produces a stem-and-leaf-display of a data set    #
    #                                                                #
    #Usage:                                                          #
    #   stem.leaf(data)                                              #
    #   stem.leaf(data,unit=100,m=5,Min=50,Max=1000,                 #
    #     rule.line=c("Dixon", "Velleman", "Sturges"),               #
    #     style=c("Tukey", "bare"), trim.outliers=TRUE, depths=TRUE, #
    #     reverse.negative.leaves=TRUE,na.rm=FALSE)                  #
    #                                                                #
    #Arguments:                                                      #
    #   data:      vector of input data                              #
    #   unit:      unit of leaves in: { ...,100,10,1,.1,.01,... }    #
    #   m:         1, 2 or 5 -- 10/m=number of possible leaf digits  #
    #   Min:       minimum of stem                                   #
    #   Max:       maximum of stem                                   #
    #   rule.line:   = "Dixon"    => number of lines <- 10*log(n,10) #
    #                = "Velleman" => number of lines <- 2*sqrt(n)    #
    #                = "Sturges"  => number of lines <- 1 + log(n,2) #
    #   style:       = "Tukey"    => Tukey-like stem ( m = 2, 5 )    #
    #   trim.outliers=TRUE        => outliers are printed absent     #
    #   depths       =TRUE        => depths info is printed          #
    #   reverse.negative.leaves=TRUE => neg.leaves are rev. sorted   #
    #Author:                                                         #
    #   Peter Wolf 05/2003 (modified slightly by J. Fox, 20 July 03) #
    #   rounding operation for comparing added 29 March 06           #
    #   07/2008 NA-values are counted if na.rm==FALSE                #
    ##################################################################

    n <- length(data <- sort(data))
    row.max <- floor(  c(Dixon   =10*log(n,10),
                         Velleman=2*sqrt(n),
                         Sturges =1+log(n,2)        ))[rule.line]

    stats <- boxplot(data,plot=FALSE)
    if(missing(Min)) Min <- if (trim.outliers) stats$stats[1,1] else min(data, na.rm=TRUE)
    if(missing(Max)) Max <- if (trim.outliers) stats$stats[5,1] else max(data, na.rm=TRUE)
    spannweite.red<-Max - Min

    zeilen.intervall.laenge <- spannweite.red / row.max
    if(missing(unit)){
           factor <- 10^ceiling(log(zeilen.intervall.laenge,10))
    } else factor <- 10^round(log(unit*10,10))
    debug.show("factor")

    z <- zeilen.intervall.laenge/factor  # z in (0.1 ,1]
    delta.tick <- c(.2,.2,.5,1)[sum(z>c(0,.1,.2,.5))]

    if(missing(m)) m <- round(1/delta.tick) else delta.tick <- 1/m
    debug.show("delta.tick"); debug.show("m")

    data.tr <- data/factor
    Min.tr  <- Min/factor
    Max.tr  <- Max/factor

    spannweite.red <- Max.tr - Min.tr
    sk.min <- floor(Min.tr)
    sk.max <- ceiling(Max.tr)
    skala  <- seq(sk.min,sk.max,by=delta.tick)
    if(sk.min<0) skala <- c(sk.min-delta.tick,skala)
    if(sk.max<0) skala <- skala[-length(skala)]
    debug.show("skala")


    lo.limit <- if (trim.outliers) skala[1] else -Inf
    lo.log   <- if(skala[1   ] <  0) data.tr <= lo.limit else data.tr <  lo.limit
    n.sk <- length(skala)
    hi.limit <- if (trim.outliers) skala[n.sk] + delta.tick else Inf
    hi.log   <- if(skala[n.sk] >= 0) data.tr >= hi.limit else data.tr >  hi.limit

    n.lower.extr.values <- sum(lo.log); n.upper.extr.values <- sum(hi.log)
    if(0<n.lower.extr.values){
      lower.line <- paste("LO:", paste(data[lo.log],collapse=" "))
    }
    if(0<n.upper.extr.values){
      upper.line <- paste("HI:", paste(data[hi.log],collapse=" "))
    }
    data.tr.red <-data.tr[(!lo.log)&(!hi.log)]


    stem  <- ifelse(data.tr.red<0, ceiling(data.tr.red), floor(data.tr.red) )
    # eps <- 1e-12; leaf <- floor(abs(data.tr.red*10-stem*10)+eps)
    leaf  <- floor(10*abs(signif(data.tr.red-stem,10)))
    debug.show("leaf"); debug.show("stem")

    class.of.data.tr<-unlist(c(
       sapply(signif(data.tr.red[data.tr.red< 0],10),
         function(x,sk)length(sk)-sum(-sk<=-x),signif(skala,10))
      ,sapply(signif(data.tr.red[data.tr.red>=0],10),
         function(x,sk)sum( sk<= x),signif(skala,10))
    ))
    debug.show("class.of.data.tr")
    class.of.data.tr  <- c(1:length(skala),class.of.data.tr)
    leaf.grouped      <- split(c(rep(-1,length(skala)),leaf),class.of.data.tr)
    leaf.grouped      <- lapply(leaf.grouped, function(x){ sort(x[-1]) })
    # debug.show("leaf.grouped")

    class.negative <- skala < 0
    class.neg.zero <- floor(skala) == -1

    if (reverse.negative.leaves){
            for (i in seq(class.negative))
                if (class.negative[i]) leaf.grouped[[i]] <- rev(leaf.grouped[[i]])
    }

    leaf.grouped.ch <- paste("|",unlist(lapply(leaf.grouped,paste,collapse="")))
    # debug.show("leaf.grouped")

    line.names <- skala
    line.names[class.negative] <- line.names[class.negative]+1
    line.names <- as.character(floor(line.names))
    line.names[class.neg.zero] <- "-0"


    if(style=="Tukey"){
      switch(as.character(m),
      "1"={},
      "2"={
            h<-round(2*(skala%%1)) #; line.names[h!=0] <- ""
            line.names<-paste(line.names,
                    ifelse(skala<0,c(".","*")[1+h],c("*",".")[1+h]),sep="")
          },
      "5"={
            h<-round(5*(skala%%1)); line.names[h>0 & h<4] <- ""
            line.names<-paste(line.names, ifelse(skala<0,
                             c(".","s","f","t","*")[1+h],
                             c("*","t","f","s",".")[1+h]), sep="")
          }
      )
    }
    ragged.left <- function(ch.lines){
      max.n <- max(n.lines<-nchar(ch.lines))
      h     <- paste(rep(" ",max.n),collapse="")
      ch.lines <- paste( substring(h,1,1+max.n-n.lines), ch.lines)
      ch.lines
    }

    line.names <- ragged.left(line.names)


    n.class <- unlist(lapply(leaf.grouped,length))
    select <- (cumsum(n.class) > 0) & rev((cumsum(rev(n.class)) > 0))
    depth     <-     cumsum(n.class)          + n.lower.extr.values
    depth.rev <- rev(cumsum(rev(n.class))     + n.upper.extr.values)
    debug.show("depth")

    uplow <- depth>=depth.rev
    pos.median <- which(uplow)[1] + (-1:0)
    h <- abs(depth[pos.median]-depth.rev[pos.median])
    pos.median <- pos.median[1]+(h[1]>h[2])
    debug.show("pos.median")

    depth[uplow] <- depth.rev[uplow]
    depth <- paste(depth,"")
    depth[pos.median] <- paste("(",n.class[pos.median],")",sep="")
    depth[n.class==0] <- " "
    depth <- if (depths) ragged.left(depth) else ""


    info <- c(  paste("1 | 2: represents",1.2*factor),
             #  paste("    m:",m     ),
                paste(" leaf unit:",factor/10),
                paste("            n:",n     ))


    stem <- paste(depth, line.names, leaf.grouped.ch)
    stem <- if((m!=5)||sum(select)>4) stem[select] else stem
    result <- list(display=stem)
    if(exists("lower.line")) result<-c(lower=lower.line,result)
    if(exists("upper.line")) result<-c(result,upper=upper.line)
    if(0<n.na&&!na.rm) result<-c(result,NAs=paste("NA's:",n.na,collapse=" "))
    result <- c(list( info=info), result)
    if(printresult){ for(i in seq(result)) cat(result[[i]],sep="\n") }
    result <- c(result, list(depths=depth, stem=line.names, leaves=leaf.grouped.ch))
    invisible(result)

  }

  data <- c(x,y)
  rule.line <- match.arg(rule.line)
  style <- match.arg(style)
  n.na  <- sum(is.na(data))
  if(0<n.na){
    data <- data[!is.na(data)]
    if(na.rm){ # data<-data[!is.na(data)]
      print("Warning: NA elements have been removed!!")
    }else{
  #   data[is.na(data)] <- mean(data,na.rm=TRUE)
  #   print("Warning: NA elements have been exchanged by the mean value!!")
    }  
  }


  if(0<length(h<-find("debug.cond")) && ".GlobalEnv" %in% h){
    debug.cond<-get("debug.cond",envir=.GlobalEnv)
  } else debug.cond<-""
  debug.show<-function(name){
    if(!exists("debug.cond")) return()
    if(debug.cond=="all"|| (name %in% debug.cond) ){
      cat(name,":\n"); obj<-eval(parse(text=name))
      if(is.vector(obj)){ print(obj) }
      return()
    }
  }

  n <- length(data <- sort(data))
  row.max <- floor(  c(Dixon   =10*log(n,10),
                       Velleman=2*sqrt(n),
                       Sturges =1+log(n,2)        ))[rule.line]

  stats <- boxplot(data,plot=FALSE)
  if(missing(Min)) Min <- if (trim.outliers) stats$stats[1,1] else min(data, na.rm=TRUE)
  if(missing(Max)) Max <- if (trim.outliers) stats$stats[5,1] else max(data, na.rm=TRUE)
  spannweite.red<-Max - Min

  zeilen.intervall.laenge <- spannweite.red / row.max
  if(missing(unit)){
         factor <- 10^ceiling(log(zeilen.intervall.laenge,10))
  } else factor <- 10^round(log(unit*10,10))
  debug.show("factor")

  z <- zeilen.intervall.laenge/factor  # z in (0.1 ,1]
  delta.tick <- c(.2,.2,.5,1)[sum(z>c(0,.1,.2,.5))]

  if(missing(m)) m <- round(1/delta.tick) else delta.tick <- 1/m
  debug.show("delta.tick"); debug.show("m")


  sl.xy <- stem.leaf(c(x,y),unit=unit,m=m,Min=Min, Max=Max, rule.line=rule.line,
                     style=style,trim.outliers=trim.outliers,depths=depths,
                     reverse.negative.leaves = reverse.negative.leaves, 
                     na.rm = na.rm, printresult=FALSE)
  h <- grep(" leaf unit: ",sl.xy$info,value=TRUE)
  unit <- as.numeric(sub(" leaf unit: ","",h))

    sl.x  <- stem.leaf(x,     unit=unit,m=m,Min=Min, Max=Max, rule.line=rule.line,
                       style=style,trim.outliers=trim.outliers,depths=depths,
                       reverse.negative.leaves = reverse.negative.leaves, na.rm = na.rm,
                       printresult=FALSE)
    sl.y  <- stem.leaf(y,     unit=unit,m=m,Min=Min, Max=Max, rule.line=rule.line,
                       style=style,trim.outliers=trim.outliers,depths=depths,
                       reverse.negative.leaves = reverse.negative.leaves, na.rm = na.rm,
                       printresult=FALSE)

  x.stem   <- gsub(" ","",sl.x$stem);          y.stem   <- gsub(" ","",sl.y$stem)
  x.leaves <- substring(sl.x$leaves,3);        y.leaves <- substring(sl.y$leaves,3)
  x.depths <- substring(sl.x$depths,3);        y.depths <- substring(sl.y$depths,3)
  x.digits <- grep("[0-9]",x.stem,value=TRUE); y.digits <- grep("[0-9]",y.stem,value=TRUE)
  h <- match(y.digits, x.digits); h <- h[!is.na(h)][1]
  # take the first ->[1] only otherwise an error occurs, see mail from John Fox 10/2013:
  x.pos <- which(x.stem==x.digits[h])[1]; y.pos <- which(y.stem==x.digits[h])[1]
  LZ <- rep(" ",d <- abs(y.pos - x.pos))
  if(x.pos < y.pos) { # x vorn verlängern 
    x.stem <- c(y.stem[1:d],x.stem); x.leaves <- c(LZ,x.leaves); x.depths <- c(LZ,x.depths)
  }
  if(y.pos < x.pos) { # y vorn verlängern 
    y.stem <- c(x.stem[1:d],y.stem); y.leaves <- c(LZ,y.leaves); y.depths <- c(LZ,y.depths)
  }
  x.l <- length(x.stem); y.l <- length(y.stem)
  LZ <- rep(" ",d <- abs(y.l-x.l))
  if(x.l < y.l) { # x hinten verlängern 
    x.stem <- c(x.stem,y.stem[-(1:x.l)]); x.leaves <- c(x.leaves,LZ)
    x.depths <- c(LZ,x.depths,LZ)
  }
  if(y.l < x.l) { # y hinten verlängern 
    y.stem <- c(y.stem,x.stem[-(1:y.l)]); y.leaves <- c(y.leaves,LZ)
    y.depths <- c(LZ,y.depths,LZ)
  }

  expand.text <- function(x,N=0,O=0,S=0,W=0,fill.right=TRUE,sep=" "){
    # bringt Elemente auf gleiche Laenge und ergaenzt Leerzeichen in den Himmelsrichtungen
    if(0<O) x <- paste(x,paste(rep(sep,O),collapse=""),sep="")
    if(0<W) x <- paste(  paste(rep(sep,W),collapse=""),x,sep="")
    if(0<S) x <- c(x,rep(sep,S)); if(0<N) x <- c(rep(sep,N),x)
    maxch <- max(nchar(x));  LZ <- paste(rep(sep,maxch),collapse="")
    ch <- substring(LZ,1,maxch-nchar(x))
    x <- if(fill.right) paste(x,ch,sep="") else paste(ch,x,sep="")
    return(x)
  } # test: cbind(expand.text("asdf",3,2,1,4,fill.right=FALSE,sep=" "))

  expand.a.to.b <- function(a,b,fill.right=TRUE,fill.tail=TRUE){
    # expandiert Vektor a auf die Dimensionen von Vektor b
    d <- length(b) - length(a)
    if(0 < d) a <- if(fill.tail) c(a, rep(" ",d)) else c(rep(" ",d), a)
    n.LZ <- max(nchar(b),nchar(a)); d <- n.LZ - nchar(a)
    LZ <- paste(rep(" ",n.LZ),collapse="")
    a <- ifelse( d <= 0, a, if(fill.right) paste(a,substring(LZ,1,d),sep="") 
                            else paste(substring(LZ,1,d),a,sep=""))
    return(a)
  }

  rotate.string <- function(x){
      x <- sapply(x,function(y) { h <- nchar(y); paste(substring(y,h:1,h:1),collapse="") })
  }

  vecpaste <- function(..., sep=" ", widths=NULL, fill.right = TRUE){  ## ???
    # vecpaste pastet Vektoren zusammen
    xyz <- list(...); n <- length(xyz)
    LZ <- paste(rep(" ",200),collapse="")
    if(0 < length(widths)) sep <- "" else widths <- rep(0,n)
    if(is.numeric(sep)) sep <- substring(LZ,1,sep)
    if(length(sep)==1)  sep <- rep(sep,n+1)
    if(length(fill.right)==1)  fill.right <- rep(fill.right,n)
    result <- ""; w <- NULL
    for(i in 1:n){ #print(fill.right[i]); cat("n",n,"i",i)
      h <- as.character(xyz[[i]]); h <- expand.text(h,fill.right=fill.right[i])
      if(nchar(h[1]) < widths[i]) 
        if(fill.right[i]){
          h <- paste(h,substring(LZ,1,widths[i]-nchar(h[1])),sep="")
        }else{
          h <- paste(substring(LZ,1,widths[i]-nchar(h[1])),h,sep="")
        }
      result <- paste(result,h,sep=sep[i])
      w <- c(w, nchar(sep[i]), nchar(h[1]))
    }
    if(0 == length(widths)) widths <- c(w,nchar(sep[n+1])) else widths <- w
    result <- paste(result,"",sep=sep[n+1])
    return(list(result=result,widths=widths))
  }

  line.to.textvec <- function(vec,width,sep=" "){
    # wandelt Textzeile in Textvektor vorgegebener maximaler Zeichenanzahl um
    vec <- paste(vec,collapse=" "); if(width==Inf) return(vec)
    d <- width - nchar(vec); LZ <- paste(rep(" ",width),collapse="")
    if( 0 < d ) return(paste(vec,substring(LZ,1,d),sep=""))
    vec.sep <- unlist(strsplit(vec,sep)); result <- NULL; imax <- 200
    while(0 < length(vec.sep)){        if( (imax<-imax-1) < 0 ) break
      h <- sum( cumsum(nchar(vec.sep)+1) <= width ) # +1 fuer Trennzeichen
      if( h==0 ) { 
        h <- 1; cat("Warning: word '",vec.sep[1],"' too long ...",sep="") 
        vec.sep <- c(substring(vec.sep[1],1,width),substring(vec.sep[1],1+width),vec.sep[-1])
        h <- 1; cat("  and has been split to '",vec.sep[1],"' and '",vec.sep[2],"'",sep="") 
      }
      new.line <- paste(vec.sep[1:h],collapse=sep)
      d <- width - nchar(new.line)
      if( d < 0 ) new.line <- substring(new.line,1,width)
      if( 0 < d ) new.line <- paste(new.line,substring(LZ,1,d),sep="")
      result <- c(result, new.line); vec.sep <- vec.sep[-(1:h)]
    }
    return(result)
  }

  x.leaves <- expand.text(x.leaves,fill.right=TRUE,O=add.more.blanks)
  x.leaves <- rotate.string(x.leaves)
  x.stem   <- expand.text(x.stem,fill.right=FALSE)
  y.leaves <- expand.text(y.leaves,fill.right=TRUE,O=add.more.blanks)
  # Zweige egalisieren
  x.leaves <- expand.a.to.b(x.leaves,y.leaves,fill.right=FALSE)
  y.leaves <- expand.a.to.b(y.leaves,x.leaves)

  if(back.to.back){ 
    if(show.no.depths){
      result <- vecpaste("",x.leaves,x.stem,y.leaves,"",sep=c("  ","","| "," |","","  ")) 
    } else {
      result <- vecpaste(x.depths,x.leaves,x.stem,y.leaves,y.depths,
                         sep=c("  ","  ","| "," |","  ","  "))
    }
    end.of.x.leaves <- sum(result$widths[1:4])
    end.of.x.attributes <- sum(result$widths[1:4])
    start.of.y.attributes <- sum(result$widths[1:7]) # falls x.name sehr lang
    space.x.to.y   <- sum(result$widths[5:7])
    end.of.left.stem <- 0 
    extr.width <- end.of.x.leaves 
  } else { # parallel displays
    if(show.no.depths){
      result <- vecpaste(x.stem,rotate.string(x.leaves),y.leaves,sep=c("  "," |","  |","  "))
      end.of.x.leaves <- sum(result$widths[1:4])
      end.of.x.attributes <- sum(result$widths[1:4])
      start.of.y.attributes <- sum(result$widths[1:4]) # falls x.name sehr lang
      space.x.to.y   <- sum(result$widths[5])
    } else {
      result <- vecpaste(x.stem,rotate.string(x.leaves),x.depths,y.leaves,y.depths,
                         sep=c("  "," |","  "," |","  ","  "))
      end.of.x.leaves <- sum(result$widths[1:6]) 
      end.of.x.attributes <- sum(result$widths[1:6]) 
      start.of.y.attributes <- sum(result$widths[1:7]) # falls x.name sehr lang
      space.x.to.y   <- sum(result$widths[7]) #- 5
    }
    end.of.left.stem <- sum(result$widths[1:3]) 
    extr.width <- end.of.x.leaves-end.of.left.stem
  }
  result <- result$result
  LZ <- paste(rep(" ",max(nchar(result))),collapse="")

  # compose info line
  info.line <- vecpaste(sl.x$info[1],sl.x$info[2],sep=c("  ",","," "))$result
  end.of.attr <- max(end.of.left.stem,4)
  # compose name line
  if( nchar(x.name) < end.of.x.leaves) {
     name.line <- vecpaste("    ",x.name," ",y.name,
          widths=c(end.of.attr, end.of.x.leaves-end.of.attr, space.x.to.y, end.of.x.leaves-4),
          fill.right=c(TRUE,!back.to.back,TRUE,TRUE))$result
  } else {
     name.line <- c(vecpaste(x.name,sep=c(end.of.x.leaves - nchar(x.name),0))$result,
                    vecpaste(y.name,sep=c(start.of.y.attributes,0))$result)
  }
  # compose length line
  n.line <- vecpaste("n:  ",length(x)," ",length(y),
          widths=c(end.of.attr, end.of.x.leaves-end.of.attr, space.x.to.y, end.of.x.leaves-4),
          fill.right=c(TRUE,!back.to.back,TRUE,TRUE))$result
  # compose NA line
  na.line <- vecpaste("NAs:",n.na.x," ",n.na.y,
          widths=c(end.of.attr, end.of.x.leaves-end.of.attr, space.x.to.y, end.of.x.leaves-4),
          fill.right=c(TRUE,!back.to.back,TRUE,TRUE))$result
  if(0==length(grep("[1-9]",na.line))) na.line <- "    "

  # compose info of HI values
  ## split HI info of x in smaller parts
  upper.x.res <- line.to.textvec(sl.x$upper, extr.width)
  ## split HI info of y in smaller parts
  upper.y.res <- line.to.textvec(sl.y$upper, extr.width)
  ## unify x and y HI values
  upper.x.res <- expand.a.to.b(upper.x.res, upper.y.res)
  upper.y.res <- expand.a.to.b(upper.y.res, upper.x.res)
  ## compose HI line
  if(back.to.back){
    upper.line <- vecpaste(upper.x.res,upper.y.res,sep=c(0,space.x.to.y,0))$result
  } else {
    upper.line <- vecpaste(" ",upper.x.res," ",upper.y.res,
       widths=c(end.of.attr,end.of.x.leaves-end.of.attr,space.x.to.y,end.of.x.leaves-4),
       fill.right=c(TRUE,TRUE,TRUE,TRUE))$result
  }    
  # compose info of LO values
  ## width=space by depths, leaves and 2 outer and 2 between
  ## split LO info of x in smaller parts
  lower.x.res <- line.to.textvec(sl.x$lower, extr.width)
  ## split LO info of y in smaller parts
  lower.y.res <- line.to.textvec(sl.y$lower, extr.width)
  ## unify x and y LO values
  lower.x.res <- expand.a.to.b(lower.x.res, lower.y.res)
  lower.y.res <- expand.a.to.b(lower.y.res, lower.x.res)
  ## compose LO line
  if(back.to.back){
    lower.line <- vecpaste(lower.x.res,lower.y.res,sep=c(0,space.x.to.y,0))$result
  } else {
    lower.line <- vecpaste(" ",lower.x.res," ",lower.y.res,
       widths=c(end.of.attr,end.of.x.leaves-end.of.attr,space.x.to.y,end.of.x.leaves-4),
       fill.right=c(TRUE,TRUE,TRUE,TRUE))$result
  }

  # komponiere Ergebnisvektor
  line   <- paste(rep("_",max(nchar(result))),collapse="")
  result <- c(line,info.line,name.line,lower.line,line,result,line,upper.line,n.line,
              if(substring(na.line,1,1)=="N") na.line,line)
  result <- result[grep("[^ ]",result)]
  # zeige Output
  if(printresult){ cat(result,sep="\n") }
  invisible(list(info=sl.x$info, display=result,lower.x=sl.x$lower,upper.x=sl.x$upper,
                 lower.y=sl.y$lower,upper.y=sl.y$upper,
                 x.depths=x.depths,y.depths=y.depths,stem=x.stem,x.leaves,y.leaves))

}

