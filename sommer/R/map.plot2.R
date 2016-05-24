map.plot2 <- function(data, trait=NULL, trait.scale="same", col.chr=NULL, col.trait=NULL, type="hist", cex=0.4, lwd=1, cex.axis=0.4, cex.trait=0.8, jump=5 ){
  ## transparent function
  ## data needs to have 2 columns; LG and Position
  ## trait needs to indicate the name to plot in the chromosome
  ## the trait can be expressed as "dot", "line" or "hist"
  ## the trait scale can be "same" or "ind", which is same for all or individual
  ## cex is only the cex for the ruler of cM
  ## cex.axis is for the titles
  transp <- function (col, alpha = 0.5) {
    res <- apply(col2rgb(col), 2, function(c) rgb(c[1]/255, c[2]/255, c[3]/255, alpha))
    return(res)
  }
  #####
  draw.arc <- function (x = 1, y = NULL, radius = 1, angle1 = deg1 * pi/180, 
            angle2 = deg2 * pi/180, deg1 = 0, deg2 = 45, n = 0.05, col = NA, 
            lwd = NA, ...) {
    getYmult<-function () {
      if (dev.cur() == 1) {
        warning("No graphics device open.")
        ymult <- 1
      }
      else {
        xyasp <- par("pin")
        xycr <- diff(par("usr"))[c(1, 3)]
        ymult <- xyasp[1]/xyasp[2] * xycr[2]/xycr[1]
      }
      return(ymult)
    }
    if (all(is.na(col))) 
      col <- par("col")
    if (all(is.na(lwd))) 
      lwd <- par("lwd")
    xylim <- par("usr")
    ymult <- getYmult()
    devunits <- dev.size("px")
    draw.arc.0 <- function(x, y, radius, angle1, angle2, n, col, 
                           lwd, ...) {
      delta.angle <- (angle2 - angle1)
      if (n != as.integer(n)) 
        n <- as.integer(1 + delta.angle/n)
      delta.angle <- delta.angle/n
      angleS <- angle1 + seq(0, length = n) * delta.angle
      angleE <- c(angleS[-1], angle2)
      if (n > 1) {
        half.lwd.user <- (lwd/2) * (xylim[2] - xylim[1])/devunits[1]
        adj.angle = delta.angle * half.lwd.user/(2 * (radius + 
                                                        half.lwd.user))
        angleS[2:n] = angleS[2:n] - adj.angle
        angleE[1:(n - 1)] = angleE[1:(n - 1)] + adj.angle
      }
      p1x <- x + radius * cos(angleS)
      p1y <- y + radius * sin(angleS) * ymult
      p2x <- x + radius * cos(angleE)
      p2y <- y + radius * sin(angleE) * ymult
      segments(p1x, p1y, p2x, p2y, col = col, lwd = lwd, ...)
    }
    xy <- xy.coords(x, y)
    x <- xy$x
    y <- xy$y
    a1 <- pmin(angle1, angle2)
    a2 <- pmax(angle1, angle2)
    angle1 <- a1
    angle2 <- a2
    args <- data.frame(x, y, radius, angle1, angle2, n, col, 
                       lwd, stringsAsFactors = FALSE)
    for (i in 1:nrow(args)) do.call("draw.arc.0", c(args[i, ], 
                                                    ...))
    invisible(args)
  }
  ## this function takes a dataframe with 2 basic columns:
  ## LG containint the linkage group
  ## Position, the position in cM
  len <- numeric()
  for(i in 1:max(unique(data$LG))){
    len[i] <- max(data[which(data$LG == i),"Position"])
  }
  coree <- which(len == max(len))
  coree1 <- data[which(data$LG == coree),]
  linesss <- coree1$Position / max(coree1$Position)
  # if one trait wants to be plotted or not
  if(!is.null(trait)){
    cols <- 1#(dim(data)[2] - 2)
  }else{cols <- 0}
  extra <- length(len) + cols *length(len)
  fact <- 1/ extra
  fact2 <- fact + (fact*cols) # real separation between chromosomes
  ## colors for traits
  if(!is.null(col.trait)){
    col.trait <- col.trait
  }else{col.trait <- c(1:6,1:6)}
  # colors for chromosomes
  
  
  plot.new()
  for(j in 1:max(unique(data$LG))){
    
    prov <- data[which(data$LG == j),] # extract the jth LG
    ## add zeros to the beggining and end of the LG so the curves look good
    dddd <- prov[1,]
    dddd2 <- prov[dim(prov)[1],]
    dddd[1,which(names(dddd) != "LG")] <- 0
    dddd2[1,which(names(dddd) != "LG" & names(dddd) != "Position")] <- 0
    prov <- rbind(dddd,prov,dddd2)
    ##-----------------------------------------------------------
    ## CHROMOSOME CHUNK
    chr <- prov$Position / max(coree1$Position)
    ### -------------------------------------------------
    ### depending if is a genetic map or physical map we decide how the ruler will be
    #if(max(prov$Position) < 1000){mark <- 10}
    #if(max(prov$Position) > 1000 & max(prov$Position) < 10000){mark <- 100}
    #if(max(prov$Position) > 10000 & max(prov$Position) < 100000){mark <- 1000}
    #if(max(prov$Position) > 100000 & max(prov$Position) < 1000000){mark <- 10000}
    #if(max(prov$Position) < 1000000){mark <- 100000}
    
    ruler <- 1 - (c(seq(0,max(prov$Position), by=jump), round(max(prov$Position),0 )) / max(coree1$Position) ); ruler2 <- c(seq(0,max(prov$Position), by=jump),round(max(prov$Position),0 ))
    if(!is.null(trait)){
      sss <- (fact2*j)-fact# + j*cols 
    }else{sss <- (fact2*j)}# + j*cols}
    
    
    ## heatmap fr density
    dd2 <- density(chr, n=length(chr))$y # regular density
    dd <- sort(density(chr, n=length(chr))$y, decreasing=T)
    ##
    if(!is.null(col.chr)){
      hc <- colorRampPalette(c(col.chr[1], col.chr[2]))( length(dd) )
    }else{hc <- gray.colors(n=length(dd), start = 0, end = 0.6, gamma = 2.2, alpha = NULL)}
    ###### for each chromosome
    for(k in 1:length(dd2)){
      # for each putative point which is the closest position
      ooo <- which(dd == dd2[k])
      lines(y=c(1-chr[k],1-chr[k]), x=c(sss,sss-(fact/3)), lwd=lwd, col=hc[ooo])
    }
    lines(y=c(1,1-max(chr)), x=c(sss,sss), lwd=3)
    lines(y=c(1,1-max(chr)), x=c(sss-(fact/3),sss-(fact/3)), lwd=3)
    text(x=sss-(fact/1.6), y=ruler, labels=ruler2, cex=cex) # cex of the cM ruler
    axis(3,at=(sss-(fact/3)) , labels=paste("LG",j, sep=""), cex.axis=cex.axis, font=2) # cex of the name of LGs
    #axis(2,at=0.275, labels="Density")
    ## --------------------------------------------------------
    draw.arc((sss + sss-(fact/3))/2, 1- max(chr), (sss - (sss-(fact/3)))/2, deg1=180, deg2=360, col="black", lwd=2, lend=1)
    draw.arc((sss + sss-(fact/3))/2, 1, (sss - (sss-(fact/3)))/2, deg1=0, deg2=180, col="black", lwd=2, lend=1)
    ## ----------------
    ## ----------------
    if(!is.null(trait)){
      ## ---------------------------
      ## IF TRAIT IS A NUMERIC TRAIT
      ## ---------------------------
      if(is.numeric(prov[,trait])){
        w1 <- which(names(data) == trait) # column of trait to plot dots
        
        if(trait.scale == "same"){
          bobo <- max(data[,trait], na.rm=TRUE)
        }else{bobo <- max(prov[,trait], na.rm=TRUE)}
        
        
        dotss <- fact * ( prov[,trait]/ bobo )
        dotss2 <- sss + (fact/8) + dotss
        ####### for ablines in the trait selected
        sections <- (bobo - min(prov[,trait], na.rm=TRUE))/5
        sections2 <- seq(min(prov[,trait], na.rm=TRUE), bobo, by=sections)
        sections3 <- fact * ( sections2/ bobo )
        sections4 <- sss + (fact/8) + sections3
        ####### if trait is true
        # ablines for different values
        for(d in 1:length(sections4)){
          lines( x=c(sections4[d], sections4[d]), y=c(1,1-max(chr)), col="black", lty=3,lwd=0.5)
          text(x=sections4[d], y=1, labels=round(sections2[d],1), cex=0.4, srt=270)
        }
      
        # plot the trait values
        if(type=="dot"){
          points(y=1-chr, x=dotss2, pch=20, cex=cex.trait, col=transp(col.trait[j],0.6))
        }
        if(type == "line"){
          polygon(y=1-chr, x=dotss2, pch=20, cex=cex.trait, col=transp(col.trait[j],0.4))
          lines(y=1-chr, x=dotss2, pch=20, cex=cex.trait, col=transp(col.trait[j],0.6))
          # density()
          
        }
        if(type == "hist"){
          for(l in 1:length(dotss2)){
            # y is the position in the chromosome
            # x is how long is the line
            lines( x=c(sss + (fact/8),dotss2[l]), y=c(1-chr[l],1-chr[l]),lwd=cex.trait, col=transp(col.trait[j],0.8))
          }
        }
        axis(3,at=sss + (fact/2), labels=trait, cex.axis=cex.axis) # cex of the scale of the trait
      }
      ## ---------------------------
      ## IF TRAIT IS A FACTOR TRAIT
      ## ---------------------------
      if(is.factor(prov[,trait])){
        riel <- sss + (fact/2)
        lines( x=c(riel, riel), y=c(1,1-max(chr)), col=transp(col.trait[j],0.8))
        ww2 <- which(!is.na(prov[,trait]))
        ##
        if(length(ww2) > 0){
          yy <- 1-chr[ww2]
          points(y=yy, x=rep(riel,length(yy)), pch=15, cex=0.9, col=transp(col.trait[j],0.6))
          text(x=riel+(fact/2), y=yy, labels=prov[ww2,trait], cex=0.3)
        }
        ##
        axis(3,at=riel, labels=trait, cex.axis=cex) # cex of the scale of the trait (letters)
      }
      ## ---------------------------
      ## ---------------------------
      
      ## ---------------------------
      ## IF TRAIT IS A character TRAIT
      ## ---------------------------
      if(is.character(prov[,trait])){
        riel <- sss + (fact/1.8)
        #lines( x=c(riel, riel), y=c(1,1-max(chr)), col=transp(col.trait[j],0.8))
        ww2 <- which(!is.na(prov[,trait]))
        ##
        if(length(ww2) > 0){
          yy <- 1-chr[ww2]
          #points(y=yy, x=rep(riel,length(yy)), pch=15, cex=0.9, col=transp(col.trait[j],0.6))
          text(x=riel, y=yy, labels=prov[ww2,trait], cex=0.17)
        }
        ##
        axis(3,at=riel, labels=trait, cex.axis=cex.axis)
      }
      ## ---------------------------
      ## ---------------------------
    }
    ################
    
  }
}

