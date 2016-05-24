


heatmap.send.legacy <- function (x,
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
                          cexRow = 0.2 +  1/log10(nr),
                          cexCol = 0.2 + 1/log10(nc),
                          labRow = NULL, 
                          labCol = NULL,
                          main = NULL,
                          xlab = NULL,
                          ylab = NULL,
                          keep.dendro = FALSE, 
                          verbose = getOption("verbose"),
                          mai.mat=NA, mai.prc=FALSE,
                          z.value="value",
                          x.lbls=NA,
                          y.lbls=NA,
                          xy.lbls=NA,
                          x.links=NA, y.links=NA,
                          xy.links=NA,asLinks=NA,
                          bound.pt = FALSE, source.plot=NA,
                          resize="800x1100",
                          ps.paper="letter",ps.width=8,ps.height=11,
                          fname.root="test",dir="./", header="v2",
                          paint=FALSE, img.prog = NA,
                          up.left=c(288,203),low.right=c(620,940),
                          spot.radius=5, automap=FALSE, automap.method="mode"
                          ) {

  
  cat("NOTE:  heatmap.send.legacy function is deprecated\n      Please see heatmap.send \n\n\n")

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

   # in all cases the first plot must be our image
   plot1 = "image(x=1:nc, y=1:nr, z=z, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = '', ylab = '');axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,cex.axis = cexCol);"
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
    tempVar = x
    

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


    # check dimensions and fix index of x.lbls,y.lbls, xy.lbls
    br = FALSE
    if(is.null(dim(x.lbls))){
      if(length(x.lbls) != 1) br = TRUE
    }
    if(!is.null(dim(x.lbls))){
      if( (dim(x.lbls)[1] == length(colInd))){
        x.lbls = x.lbls[colInd,]
      }else{
        br = TRUE
      }
    }
    if(br){
      cat(paste("x.lbls dimension is not correct. \n x.lbls should be a matrix with number of rows equal to ",nc ," \n continuing with x.lbls =NA \n", sep=""))
      x.lbls = NA
    }
    br = FALSE
    if(is.null(dim(y.lbls))){
      if(length(y.lbls) != 1) br = TRUE
    }
    if(!is.null(dim(y.lbls))){
      if( (dim(y.lbls)[1] == length(rowInd))){
        y.lbls = y.lbls[rowInd,]
      }else{
        br = TRUE
      }
    }
    if(br){
      cat(paste("y.lbls dimension is not correct. \n y.lbls should be a matrix with number of rows equal to ",nr ," \n continuing with y.lbls =NA \n", sep=""))
      y.lbls = NA
    }

    br = FALSE
    if(!is.null(dim(xy.lbls))){
      if(dim(xy.lbls)[1]!=nr) {
        br = TRUE
      }
      if(dim(xy.lbls)[2]!=nc) {
        
         br = TRUE
      }
      if(br){
        cat(paste("xy.lbls dimension is not correct. \n xy.lbls should be a matrix with \n number of rows equal to ",nr ," \n and the number of columns equal to ", nc, "\n continuing with y.lbls =NA \n", sep=""))
        xy.lbls=NA
      } 
      if(!is.null(dim(xy.lbls))){
        newDF = list()
        newDF$xy.lbl = xy.lbls
        xy.lbls = newDF
      }    
    }
    br = FALSE
    if(length(xy.lbls)==1){
      if(!is.na(xy.lbls[1])){
        if(dim(xy.lbls[[1]])[1]==nr) {
          xy.lbls[[1]] = xy.lbls[[1]][rowInd,]
        }else{
          br = TRUE
        }        
        if(dim(xy.lbls[[1]])[2]==nc) {
          xy.lbls[[1]] = xy.lbls[[1]][,colInd]
        }else{
          br = TRUE
        }
        if(br){
          cat(paste("xy.lbls dimension is not correct. \n xy.lbls should be a matrix with \n number of rows equal to ",nr ," \n and the number of columns equal to ", nc, "\n continuing with y.lbls =NA \n", sep=""))
          xy.lbls=NA
        } 
      }
    }       
    br = FALSE
    if(length(xy.lbls)>1){
      for(k in 1:length(xy.lbls)){
        if(dim(xy.lbls[[k]])[1]==nr) xy.lbls[[k]] = xy.lbls[[k]][rowInd,]
        if(dim(xy.lbls[[k]])[1]!=nr) br = TRUE
        if(dim(xy.lbls[[k]])[2]==nc) xy.lbls[[k]] = xy.lbls[[k]][,colInd]
        if(dim(xy.lbls[[k]])[2]!=nc) br = TRUE
      }
    }
    if(br){
      cat(paste("an xy.lbls dimension is not correct.\n number of rows should be equal to ",nr ," \n and the number of columns equal to ", nc, "\n continuing with xy.lbls =NA \n", sep=""))
      xy.lbls = NA
    }

    # check dimensions and fix index of x.links,y.links, xy.links, asLinks
    br = FALSE
    if(is.null(dim(x.links))){
      if(length(x.links) != 1) br = TRUE
    }
    if(!is.null(dim(x.links))){
      if( (dim(x.links)[1] == length(colInd))){
        x.links = x.links[colInd,]
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
      cat(paste("asLinks is not of correct dimension or length \n continuing with asLinks =NA \n", sep=""))
      asLinks = NA
  
    }
    
    # load required library
    require("sendplot")

    # set environment to use variables in memory
    environment(sendplot) <- environment()

    # call sendplot 
    sendplot(mat=mat, x=sendX, y=sendY,plot.calls = plt.calls, z=tempVar, type="image",mai.mat=mai.mat, mai.prc=mai.prc, plt.extras = plot.extras, x.lbls=x.lbls, y.lbls=y.lbls, xy.lbls=xy.lbls,x.links=x.links, y.links=y.links, xy.links=xy.links, asLinks=asLinks,bound.pt=bound.pt, source.plot=source.plot, z.value=z.value, resize=resize, ps.paper=ps.paper, ps.width=ps.width,ps.height=ps.height,fname.root=fname.root,dir=dir, header= header,paint=paint,img.prog=img.prog,up.left=up.left,low.right=low.right,spot.radius=spot.radius, automap=automap, automap.method=automap.method)

       
}
