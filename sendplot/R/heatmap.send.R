#
# heatmap.send wrapper to default heatmap function in R 
#



heatmap.send <- function (x,
                          Rowv = NULL,
                          Colv = if (symm) "Rowv" else NULL, 
                          distfun = dist,
                          hclustfun = hclust,
                          reorderfun = function(d,w) reorder(d, w),
                          add.expr,
                          symm = FALSE,
                          revC = identical(Colv,"Rowv"),
                          scale = c("row", "column", "none"),
                          na.rm = TRUE, 
                          margins = c(5, 5),
                          ColSideColors,
                          RowSideColors,
                          MainColor = heat.colors(12),
                          cexRow = 0.2 +  1/log10(nr),
                          cexCol = 0.2 + 1/log10(nc),
                          labRow = NULL, 
                          labCol = NULL,
                          main = NULL,
                          xlab = NULL,
                          ylab = NULL,
                          keep.dendro = FALSE, 
                          verbose = getOption("verbose"),

                          x.labels=NA,
                          y.labels=NA,
                          xy.labels=NA,
                          x.links=NA,
                          y.links=NA,
                          xy.links=NA,
                          asLinks=NA,
                          x.images=NA,
                          y.images=NA,
                          xy.images=NA,
                          spot.radius=5,
                          source.plot=NA,
                          image.size="800x1100",
                          fname.root="test",dir="./", header="v3",
                          window.size = "800x1100", # in px

                          ...
                          ) {


   # code from heatmap function
  
    scale <- if (symm && missing(scale)) 
        "none"
    else match.arg(scale)
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("'x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) 
        stop("'x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2)
        stop("'margins' must be a numeric vector of length 2")
    doRdend <- !identical(Rowv, NA)
    doCdend <- !identical(Colv, NA)
    if (is.null(Rowv)) 
        Rowv <- rowMeans(x, na.rm = na.rm)
    if (is.null(Colv)) 
        Colv <- colMeans(x, na.rm = na.rm)
    if (doRdend) {
        if (inherits(Rowv, "dendrogram")) 
            ddr <- Rowv
        else {
            hcr <- hclustfun(distfun(x))
            ddr <- as.dendrogram(hcr)
            if (!is.logical(Rowv) || Rowv) 
                ddr <- reorderfun(ddr, Rowv)
        }
        if (nr != length(rowInd <- order.dendrogram(ddr))) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else rowInd <- 1:nr
    if (doCdend) {
        if (inherits(Colv, "dendrogram")) 
            ddc <- Colv
        else if (identical(Colv, "Rowv")) {
            if (nr != nc) 
                stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
            ddc <- ddr
        }
        else {
            hcc <- hclustfun(distfun(if (symm) 
                x
            else t(x)))
            ddc <- as.dendrogram(hcc)
            if (!is.logical(Colv) || Colv) 
                ddc <- reorderfun(ddc, Colv)
        }
        if (nc != length(colInd <- order.dendrogram(ddc))) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else colInd <- 1:nc
    x <- x[rowInd, colInd]
    labRow <- if (is.null(labRow)) 
        if (is.null(rownames(x))) 
            (1:nr)[rowInd]
        else rownames(x)
    else labRow[rowInd]
    labCol <- if (is.null(labCol)) 
        if (is.null(colnames(x))) 
            (1:nc)[colInd]
        else colnames(x)
    else labCol[colInd]
    if (scale == "row") {
        x <- sweep(x, 1, rowMeans(x, na.rm = na.rm))
        sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        x <- sweep(x, 2, colMeans(x, na.rm = na.rm))
        sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    lmat <- rbind(c(NA, 3), 2:1)
    lwid <- c(if (doRdend) 1 else 0.05, 4)
    lhei <- c((if (doCdend) 1 else 0.05) + if (!is.null(main)) 0.2 else 0, 
        4)

 
    if (!missing(ColSideColors)) {
        if (!is.character(ColSideColors) || length(ColSideColors) != 
            nc) 
            stop("'ColSideColors' must be a character vector of length ncol(x)")
        lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
        lhei <- c(lhei[1], 0.2, lhei[2])
    }
    if (!missing(RowSideColors)) {
        if (!is.character(RowSideColors) || length(RowSideColors) != 
            nr) 
            stop("'RowSideColors' must be a character vector of length nrow(x)")
        lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 
            1), lmat[, 2] + 1)
        lwid <- c(lwid[1], 0.2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
    if (!symm || scale != "none") 
        x <- t(x)
    if (revC) {
        iy <- nr:1
        ddr <- rev(ddr)
        x <- x[, iy]
    }
    else iy <- 1:nr

    # now move image and plotting calls to match structures of sendplot
    mai.mat=NA
    mai.prc=FALSE

    
   # in all cases the first plot must be our image
   plot1 = "image(x=1:nc, y=1:nr, z=z, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = '', ylab = '', col=MainColor);axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,cex.axis = cexCol);"
    if (!is.null(xlab)) plot1 = paste(plot1,"mtext(xlab, side = 1, line = margins[1] - 1.25);", sep="")
    plot1 = paste(plot1, "axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,cex.axis = cexRow);", sep="")
     if (!is.null(ylab)) plot1 = paste(plot1,"mtext(ylab, side = 4, line = margins[2] - 1.25);",sep="")
     plot.extras=list(); plot.extras$plot1=NA;plot.extras$plot2=NA;plot.extras$plot3=NA
    if (!missing(add.expr)) plot.extras$plot1=add.expr
    
    # case 1: no color band

    #
    # adjust all future so that ddr are represented row = 2 col = 3 
    # 
    
    # if row dendrogram included
    if(doRdend){
      plot2 = "plot(ddr, horiz = TRUE, axes = FALSE, yaxs = 'i', leaflab = 'none')"
    }else{
      plot2 = "frame()"
    }
    # if column dendrogram include
    if (doCdend){
      plot3 = "plot(ddc, axes = FALSE, xaxs = 'i', leaflab = 'none');title(main=main, cex=1.5)"
    }else{
      plot3 = "frame();title(main=main, cex=1.5)"
    }

    plot4 = NA
    plot5 = NA
  
    # case2: one color
    if(max(lmat, na.rm=TRUE)==4){
      # if col color
      if((dim(lmat)[1])==3 & (dim(lmat)[2])==2){
        #adjust matrix layout so that image = 1, row ddr = 2,and col.ddr=3
        lmat[,2] = c(3,4,1)
        lmat[3,1] = 2        
        plot4 = "image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)"
        plot.extras$plot4=NA
      }# end if col color
      # if row color
      if((dim(lmat)[1])==2 & (dim(lmat)[2])==3){
        #adjust matrix layout so that image = 1, row ddr = 2,and col.ddr=3
        lmat[2,] = c(2,4,1)
        lmat[1,3] = 3
        plot4 = "image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)"
        plot.extras$plot4=NA
      }# end if row color
    }# end case 2: one color
        
    # case3: both color
    if(max(lmat, na.rm=TRUE)==5){
        #adjust matrix layout so that image = 1, row ddr = 2,and col.ddr=3
        lmat[,3] = c(3,5,1)
        lmat[3,] = c(2,4,1)
        plot5 = "image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)"
        plot.extras$plot5=NA
        plot4 = "image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)"
        plot.extras$plot4=NA
    }# end case 3: both color

    # initialize some variables for sendplot function call
    plt.calls = c(plot1, plot2, plot3)
    if(!is.na(plot4)) plt.calls = c(plt.calls, plot4)
    if(!is.na(plot5)) plt.calls = c(plt.calls, plot5)
    sendX = 1:nc
    sendY = 1:nr
    z = x
    

    # now default layout makes equally spaced regions.
    # we want to change this option for more  appealing
    mat=matrix(c(rep(c(rep(1,10)),17)),ncol=10,byrow=TRUE)
    if(max(lmat, na.rm=TRUE)==3) {
      mat[1:2,] = 3
      mat[,1:2] = 2
      mat[1:2,1:2] = 0

      if(is.na(mai.mat)){
        mai.mat = matrix(0,nrow=3,ncol=4)
        mai.mat[1,]= c(.75,  .1,  .1,  .75)
        mai.mat[2,]= c(.75,  0,   .1,   0)
        mai.mat[3,]= c( 0,   .1,  .1,   .75)
      }        
    }
   if(max(lmat, na.rm=TRUE)==4){
      # row
      if((dim(lmat)[1])==2 & (dim(lmat)[2])==3){
        mat[1:2,] = 3
        mat[,1:2] = 2
        mat[3:(dim(mat)[1]),3] = 4
        mat[1:2,1:3] = 0

        if(is.na(mai.mat)){
          mai.mat = matrix(0,nrow=4,ncol=4)
          mai.mat[1,]= c(.75,   .05,   .05,   .75)
          mai.mat[2,]= c(.75,   0,     .05,    0)
          mai.mat[3,]= c( 0,    .05,   .1,     .75)
          mai.mat[4,]= c(.75,   .05,   .05,    0)
        }                
      }      
      # col
      if((dim(lmat)[1])==3 & (dim(lmat)[2])==2){
        mat[1:2,] = 3
        mat[,1:2] = 2
        mat[3,3:(dim(mat)[2])] = 4
        mat[1:3,1:2] = 0

        
        if(is.na(mai.mat)){
          mai.mat = matrix(0,nrow=4,ncol=4)
          mai.mat[1,]= c(.75,   .05,   .025,   .75)
          mai.mat[2,]= c(.75,   0,     .025,     0)
          mai.mat[3,]= c(0,     .05,    .1,     .75)
          mai.mat[4,]= c(.025,   .05,  .05,    .75)
        }   
      }               
    }
    if(max(lmat, na.rm=TRUE)==5){
      mat[1:2,] = 3
      mat[,1:2] = 2
      mat[3:(dim(mat)[1]),3] = 4
      mat[3,3:(dim(mat)[2])] = 5
      mat[1:3,1:3] = 0

      if(is.na(mai.mat)){
        mai.mat = matrix(0,nrow=5,ncol=4)
        mai.mat[1,]= c(.75,   .05,   .05,   .75)
        mai.mat[2,]= c(.75,   0,     .05,     0)
        mai.mat[3,]= c( 0,    .05,    .1,     .75)
        mai.mat[4,]= c(.75,   .05,   .05,     0)
        mai.mat[5,]= c(.025,  .05,   .05,    .75)
      }   
      
    }

    Splot = initSplot(mat=mat, plot.calls= plt.calls, mai.mat=mai.mat, Iflag=c(TRUE,rep(FALSE, (length(plt.calls)-1))), figTypes=c("image",rep(FALSE, (length(plt.calls)-1))),source.plot="png",  image.size=image.size, returnVl=TRUE, saveFlag=FALSE)
    
    # check dimensions and fix index of x.labels,y.labels, xy.labels
    br = FALSE
    if(is.null(dim(x.labels))){
      if(length(x.labels) != 1) br = TRUE
    }
    if(!is.null(dim(x.labels))){
      if( (dim(x.labels)[1] == length(colInd))){
        xnm = names(x.labels)
        x.labels = as.data.frame(x.labels[colInd,])
        names(x.labels) = xnm
      }else{
        br = TRUE
      }
    }
    if(br){
      cat(paste("x.labels dimension is not correct. \n x.labels should be a matrix with number of rows equal to ",nc ," \n continuing with x.labels =NA \n", sep=""))
      x.labels = NA
    }
    br = FALSE
    if(is.null(dim(y.labels))){
      if(length(y.labels) != 1) br = TRUE
    }
    if(!is.null(dim(y.labels))){
      if( (dim(y.labels)[1] == length(rowInd))){
        y.labels = y.labels[rowInd,]
      }else{
        br = TRUE
      }
    }
    if(br){
      cat(paste("y.labels dimension is not correct. \n y.labels should be a matrix with number of rows equal to ",nr ," \n continuing with y.labels =NA \n", sep=""))
      y.labels = NA
    }

    br = FALSE
    if(!is.null(dim(xy.labels))){
      if(dim(xy.labels)[1]!=nr) {
        br = TRUE
      }
      if(dim(xy.labels)[2]!=nc) {
        
         br = TRUE
      }
      if(br){
        cat(paste("xy.labels dimension is not correct. \n xy.labels should be a matrix with \n number of rows equal to ",nr ," \n and the number of columns equal to ", nc, "\n continuing with y.labels =NA \n", sep=""))
        xy.labels=NA
      } 
      if(!is.null(dim(xy.labels))){
        newDF = list()
        newDF$xy.lbl = xy.labels
        xy.labels = newDF
      }    
    }
    br = FALSE
    if(length(xy.labels)==1){
      if(!is.na(xy.labels[1])){
        if(dim(xy.labels[[1]])[1]==nr) {
          xy.labels[[1]] = xy.labels[[1]][rowInd,]
        }else{
          br = TRUE
        }        
        if(dim(xy.labels[[1]])[2]==nc) {
          xy.labels[[1]] = xy.labels[[1]][,colInd]
        }else{
          br = TRUE
        }
        if(br){
          cat(paste("xy.labels dimension is not correct. \n xy.labels should be a matrix with \n number of rows equal to ",nr ," \n and the number of columns equal to ", nc, "\n continuing with y.labels =NA \n", sep=""))
          xy.labels=NA
        } 
      }
    }       
    br = FALSE
    if(length(xy.labels)>1){
      for(k in 1:length(xy.labels)){
        if(dim(xy.labels[[k]])[1]==nr) xy.labels[[k]] = xy.labels[[k]][rowInd,]
        if(dim(xy.labels[[k]])[1]!=nr) br = TRUE
        if(dim(xy.labels[[k]])[2]==nc) xy.labels[[k]] = xy.labels[[k]][,colInd]
        if(dim(xy.labels[[k]])[2]!=nc) br = TRUE
      }
    }
    if(br){
      cat(paste("an xy.labels dimension is not correct.\n number of rows should be equal to ",nr ," \n and the number of columns equal to ", nc, "\n continuing with xy.labels =NA \n", sep=""))
      xy.labels = NA
    }

    # check dimensions and fix index of x.links,y.links, xy.links, asLinks
    br = FALSE
    if(is.null(dim(x.links))){
      if(length(x.links) != 1) br = TRUE
    }
    if(!is.null(dim(x.links))){
      if( (dim(x.links)[1] == length(colInd))){
        xlnm = names(x.links)
        x.links = as.data.frame(x.links[colInd,])
        names(x.links)=xlnm
      }else{
        br = TRUE
      }
    }
    if(br){
      cat(paste("x.links dimension is not correct. \n x.links should be a matrix with number of rows equal to ",nc ," \n continuing with x.links =NA \n", sep=""))
      x.links = NA
    }
    br = FALSE
    if(is.null(dim(y.links))){
      if(length(y.links) != 1) br = TRUE
    }
    if(!is.null(dim(y.links))){
      if( (dim(y.links)[1] == length(rowInd))){
        y.links = y.links[rowInd,]
      }else{
        br = TRUE
      }
    }
    if(br){
      cat(paste("y.links dimension is not correct. \n y.links should be a matrix with number of rows equal to ",nr ," \n continuing with y.links =NA \n", sep=""))
      y.links = NA
    }

    br = FALSE
    if(!is.null(dim(xy.links))){
      if(dim(xy.links)[1]!=nr) {
        br = TRUE
      }
      if(dim(xy.links)[2]!=nc) {
        
         br = TRUE
      }
      if(br){
        cat(paste("xy.links dimension is not correct. \n xy.links should be a matrix with \n number of rows equal to ",nr ," \n and the number of columns equal to ", nc, "\n continuing with y.links =NA \n", sep=""))
        xy.links=NA
      } 
      if(!is.null(dim(xy.links))){
        newDF = list()
        newDF$xy.lbl = xy.links
        xy.links = newDF
      }    
    }
    br = FALSE
    if(length(xy.links)==1){
      if(!is.na(xy.links[1])){
        if(dim(xy.links[[1]])[1]==nr) {
          xy.links[[1]] = xy.links[[1]][rowInd,]
        }else{
          br = TRUE
        }        
        if(dim(xy.links[[1]])[2]==nc) {
          xy.links[[1]] = xy.links[[1]][,colInd]
        }else{
          br = TRUE
        }
        if(br){
          cat(paste("xy.links dimension is not correct. \n xy.links should be a matrix with \n number of rows equal to ",nr ," \n and the number of columns equal to ", nc, "\n continuing with xy.links =NA \n", sep=""))
          xy.links=NA
        } 
      }
    }       
    br = FALSE
    if(length(xy.links)>1){
      for(k in 1:length(xy.links)){
        if(dim(xy.links[[k]])[1]==nr) xy.links[[k]] = xy.links[[k]][rowInd,]
        if(dim(xy.links[[k]])[1]!=nr) br = TRUE
        if(dim(xy.links[[k]])[2]==nc) xy.links[[k]] = xy.links[[k]][,colInd]
        if(dim(xy.links[[k]])[2]!=nc) br = TRUE
      }
    }
    if(br){
      cat(paste("an xy.links dimension is not correct.\n number of rows should be equal to ",nr ," \n and the number of columns equal to ", nc, "\n continuing with yx.links =NA \n", sep=""))
      xy.links = NA
    }
    
    br = FALSE
    if(is.null(dim(asLinks))){
      if(length(asLinks) == (length(rowInd)*length(colInd))) asLinks = matrix(asLinks, ncol=length(colInd))
    }    
    if(!is.null(dim(asLinks))){
      if( (dim(asLinks)[1] == length(rowInd))){
        asLinks = asLinks[rowInd,]
      }else{
        br = TRUE
      }
      if( (dim(asLinks)[2] == length(colInd))){
        asLinks = asLinks[,colInd]
      }else{
        br = TRUE
      }
    }
    if(is.null(dim(asLinks))){
      if(length(asLinks) == 1){
        if(!is.na(asLinks)) asLinks = rep(asLinks, (length(rowInd)*length(colInd)))
      }
      if(length(asLinks) == length(rowInd)) asLinks = rep(asLinks, length(colInd))
      if(length(asLinks) == length(colInd)) asLinks = rep(asLinks, length(rowInd))
      
      if(length(asLinks) != (length(rowInd)*length(colInd))) br = TRUE
    }
    if(br){
      if(!is.na(asLinks[1])) cat(paste("asLinks is not of correct dimension or length \n continuing with asLinks =NA \n", sep=""))
      asLinks = NA
  
    }



   # check dimensions and fix index of x.images,y.images, xy.images
    br = FALSE
    if(is.null(dim(x.images))){
      if(length(x.images) != 1) br = TRUE
    }
    if(!is.null(dim(x.images))){
      if( (dim(x.images)[1] == length(colInd))){
        xinm = names(x.images)
        x.images = as.data.frame(x.images[colInd,])
        names(x.imagse) = xinm
      }else{
        br = TRUE
      }
    }
    if(br){
      cat(paste("x.images dimension is not correct. \n x.images should be a matrix with number of rows equal to ",nc ," \n continuing with x.images =NA \n", sep=""))
      x.images = NA
    }
    br = FALSE
    if(is.null(dim(y.images))){
      if(length(y.images) != 1) br = TRUE
    }
    if(!is.null(dim(y.images))){
      if( (dim(y.images)[1] == length(rowInd))){
        y.images = y.images[rowInd,]
      }else{
        br = TRUE
      }
    }
    if(br){
      cat(paste("y.images dimension is not correct. \n y.images should be a matrix with number of rows equal to ",nr ," \n continuing with y.images =NA \n", sep=""))
      y.images = NA
    }

    br = FALSE
    if(!is.null(dim(xy.images))){
      if(dim(xy.images)[1]!=nr) {
        br = TRUE
      }
      if(dim(xy.images)[2]!=nc) {
        
         br = TRUE
      }
      if(br){
        cat(paste("xy.images dimension is not correct. \n xy.images should be a matrix with \n number of rows equal to ",nr ," \n and the number of columns equal to ", nc, "\n continuing with y.images =NA \n", sep=""))
        xy.images=NA
      } 
      if(!is.null(dim(xy.images))){
        newDF = list()
        newDF$xy.lbl = xy.images
        xy.images = newDF
      }    
    }
    br = FALSE
    if(length(xy.images)==1){
      if(!is.na(xy.images[1])){
        if(dim(xy.images[[1]])[1]==nr) {
          xy.images[[1]] = xy.images[[1]][rowInd,]
        }else{
          br = TRUE
        }        
        if(dim(xy.images[[1]])[2]==nc) {
          xy.images[[1]] = xy.images[[1]][,colInd]
        }else{
          br = TRUE
        }
        if(br){
          cat(paste("xy.images dimension is not correct. \n xy.images should be a matrix with \n number of rows equal to ",nr ," \n and the number of columns equal to ", nc, "\n continuing with xy.images =NA \n", sep=""))
          xy.images=NA
        } 
      }
    }       
    br = FALSE
    if(length(xy.images)>1){
      for(k in 1:length(xy.images)){
        if(dim(xy.images[[k]])[1]==nr) xy.images[[k]] = xy.images[[k]][rowInd,]
        if(dim(xy.images[[k]])[1]!=nr) br = TRUE
        if(dim(xy.images[[k]])[2]==nc) xy.images[[k]] = xy.images[[k]][,colInd]
        if(dim(xy.images[[k]])[2]!=nc) br = TRUE
      }
    }
    if(br){
      cat(paste("an xy.images dimension is not correct.\n number of rows should be equal to ",nr ," \n and the number of columns equal to ", nc, "\n continuing with yx.images =NA \n", sep=""))
      xy.images = NA
    }







    

   # save(list=ls(), file="Heatmap.Helper.RData", compress=T) 


    # set environments so can access local workspace
    environment(addBounding) <- environment()
    environment(addDefault) <- environment()
    environment(automapPts) <- environment()
    environment(getBounds) <- environment()
    #environment(getPlotsBounds) <- environment()
    environment(initSplot) <- environment()
    environment(makeCharacter) <- environment()
    environment(makeImageDF) <- environment()
    environment(makeImap) <- environment()
    environment(makePolyDF) <- environment()
    environment(makeRectDF) <- environment()
    environment(makeScatterDF) <- environment()
    environment(makeSplot) <- environment()
    environment(mapMethod) <- environment()
    environment(removeImap) <- environment()
    environment(writeDefault1) <- environment()
    environment(writeDefault2) <- environment()
    environment(writeCircle.1) <- environment()
    environment(writeRect.1) <- environment()
    environment(writePoly.1) <- environment()
    environment(writeCircle.2) <- environment()
    environment(writeRect.2) <- environment()
    environment(writePoly.2) <- environment()
    environment(writeToHTML1) <- environment()
    environment(writeToHTML2) <- environment()
    
    # make boxes interactive not midpoint
    x.pos = c( (sendX-0.5), (sendX[length(sendX)]+0.5) )
    y.pos = c( (sendY-0.5), (sendY[length(sendY)]+0.5) )
    
    Splot = makeImap(Splot, figure=1, xy.type="image.box", x.pos=x.pos, y.pos=y.pos, spot.radius=spot.radius, x.labels=x.labels, y.labels=y.labels, xy.labels=xy.labels, x.links=x.links, y.links=y.links, xy.links=xy.links, asLinks=asLinks, x.images=x.images, y.images=y.images, xy.images=xy.images, fname.root=fname.root, dir=dir, returnVl=TRUE, ...)
    
    Splot = makeSplot(Splot, fname.root=fname.root, dir=dir, window.size=window.size, returnObj=TRUE)

    
    
  }
