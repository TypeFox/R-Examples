plotNeighbourLD <- function(LD,gpData,dense=FALSE,nMarker=TRUE,centr=NULL,file=NULL,fileFormat="pdf",...){
    oldPar <- par()
    if (class(gpData) == "gpData"){
       gpData.unit <- gpData$info$map.unit
       map <- gpData$map
       class(map) <- "data.frame"
    }
    else map.unit <- "unit"
    chr <- unique(map$chr)
    chr <- chr[!is.na(chr)]
    map <- map[!is.na(map$chr), ]
    if(class(map$chr) != "factor") bord <- NULL else if(par()$bg == "transparent") bord <- "white" else bord <- "transparent"

    # centromere positions of maize
    if(!is.null(centr)) if(centr == "maize") centr <- c(133,90,95,104.6,105.5,50,55.3,46.5,68.8,59.9)

    # norm pos
    if (!is.null(centr)) map$pos <- map$pos - centr[map$chr]

    # output in files
    if(!is.null(file)){
      if(substr(file, nchar(file)-nchar(fileFormat)+1, nchar(file)) != fileFormat | nchar(file) < 5)
        file <- paste(file, ".", fileFormat, sep="")
      if(fileFormat == "pdf") pdf(file)
      else if (fileFormat == "png") png(file)
      else stop("not supported file format choosen!")
    }
    
    # initialize map
    layout(matrix(1:2,ncol=2),widths=c(0.82,0.18))
    
    # make an empty plot 
    if(!is.null(centr)) {
    plot(map, type = "n", xaxt = "n", xlim = c(0.5, length(chr) +0.5), border=bord,
    ylim = c( max(map$pos,na.rm = TRUE) * 1.1, min(map$pos, na.rm = TRUE)),axes=FALSE, ...)
    }
    else{
    plot(map, type = "n", xaxt = "n", xlim = c(0.5, length(chr) + 0.5), border=bord,
         ylim = c( max(map$pos,na.rm = TRUE) * 1.1, min(map$pos, na.rm = TRUE)), ...)
    }   
 
   # x-axis     
    axis(side = 1, at = seq(along = chr), labels = chr)
   # y-axis
    if(!is.null(centr)){
        box()
        axis(side=2,at=-seq(-round(max(map$pos, na.rm = TRUE),-2),round(max(map$pos, na.rm = TRUE),-2),by=25),labels=abs(-seq(-round(max(map$pos, na.rm = TRUE),-2),round(max(map$pos, na.rm = TRUE),-2),by=25)),las=1)
    }

    # loop over chromosomes
    for (i in seq(along=chr)){
        mNam <- rownames(map[map$chr == chr[i],])
        LDm <- matrix(NA, ncol=10, nrow=length(mNam))
        if(class(LD) == "LDdf"){
            for(j in 1:10)
                LDm[1:(nrow(LDm)-j), j] <- LD[[i]]$r2[paste(LD[[i]]$marker1, LD[[i]]$marker2, sep="") %in% paste(mNam[j:(length(mNam)-1)], mNam[(1+j):length(mNam)], sep="")]
        } else if(class(LD) == "LDmat"){
            for(j in 1:10)
                LDm[1:(nrow(LDm)-j), j] <- diag(LD$LD[[i]][j:(nrow(LDm)-1),(1+j):nrow(LDm)])
        } else stop(paste("Wrong class of object", substitute(LD)))
    
        n <- sum(map$chr==i,na.rm=TRUE) 
        start <- min(map$pos[map$chr==chr[i]],na.rm=TRUE)
        end <- max(map$pos[map$chr==chr[i]],na.rm=TRUE)
      
        # display.brewer.pal(7, "Reds")s
        cols <- c("#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#99000D")

    if(dense){  # calculating averaged LD
        # map positions
        mapPos <- smooth(map$pos[map$chr == chr[i]])
        # matrices with 10 rows
        avLD <- rowMeans(LDm[1:(nrow(LDm)-11),])
        avPos <- matrix(NA, ncol=10, nrow=nrow(map[map$chr == i,]))
        for(j in 1:10)
            avPos[1:(nrow(LDm)-j), j] <- mapPos[1:(length(mapPos)-j)]
        avPos <- rowMeans(avPos[1:(nrow(avPos)-11),])

    }  # end of if,  smooth LD calculation
    else{  # using LD directly
        # map posittions
        mapPos <- map$pos[map$chr==chr[i]]
        avLD <- LDm[-nrow(LDm),1]
        avPos <- mapPos[-length(mapPos)]
    }
    LDPos <- data.frame(avLD,avPos)
    LDPos <- LDPos[!duplicated(LDPos$avPos), ]
    #LDPos <- LDPos[LDPos$avLD>=0.2, ]
    # visualisation of map
    image(seq(i-0.35,i+0.35,length=20), LDPos$avPos,matrix(rep(LDPos$avLD,20),nrow=20,byrow=TRUE),col=cols,add=TRUE,zlim=c(0,1))
    if(!is.null(centr)){
        # centromere
        polygon(x=c(i-0.4,i-0.1,i-0.1,i-0.4,i-0.4),y=c(-10,-1,1,10,-10),col="white",border="white")
        polygon(x=c(i+0.4,i+0.1,i+0.1,i+0.4,i+0.4),y=c(-10,-1,1,10,-10),col="white",border="white")
    }
    if(nMarker) text(i,max(map$pos)*1.05,sum(map$chr==chr[i],na.rm=TRUE))
 
  } # end chromosome loop
    
    # add legend
        par(mar = c(5, 1, 4, 3.8) + 0.1)
      image(seq(-0.4,0.4,length=20),seq(from=0,to=1,length=6),matrix(rep(seq(from=0,to=1,length=6),20),nrow=20,byrow=TRUE),main=expression(paste(r^{2})),col=cols,axes=FALSE,xlab="")
      axis(side=4,at=round(seq(from=0,to=1,length=6),4),las=1)
    # close graphic device
    if(!is.null(file)) dev.off() else {
      oldPar$cin <- oldPar$cra <- oldPar$csi <- oldPar$cxy <- oldPar$din <- NULL
      par(oldPar)
    }

}
