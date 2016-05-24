############################################
##  R package pairheatmap
##  compare heatmaps
############################################


# para
pairheatmap <- function(   data1, data2,
                    scale="none", dendrogram="both",
                    matDist=0.5, matrixBorderCol = "grey",
                    colorStyle="s1",
                    rowGroup=rowGroup,
                    orderRowGroup=NULL,
                    rowGroupColor=FALSE,
                    rowGroupColor.choice,
                    groupBorder = "line",
                    groupBorder.selectList=list(),
                    groupBorder.lwd=3,
                    groupBorder.col="green",
                    rowNameColor="blue",
                    colNameColor="blue",
                    rowNameFontSize=7,
                    colNameFontSize=7,
                    rowNameGroupColor=NULL,
                    #colNameGroupColor=FALSE,
                    clusterMethod="complete", clusterMembers=NULL,
                    clusterRow=TRUE, clusterCol=TRUE, clusterColTogether=FALSE,
                    legend.pos="middle", legend.percent=0.5, legend.fontsize=7
                    )
{ 

## parameter checking
    if (missing(data2)) stop("data2 is required!")
    if (nrow(data1)==0 || ncol(data1)==0) stop("data1 must NOT be 0 row/col")
    if (nrow(data2)==0 || ncol(data2)==0) stop("data2 must NOT be 0 row/col")
    if (nrow(data1) != nrow(data2)) stop("data1 must have same number of row as data2!")
    
    if (!(scale %in% c("row", "col", "none", "rowsep", "colsep")))
        stop("scale can only be none, row, col!")
    if (!(dendrogram %in% c("row", "col", "both")))
        stop("dendrogram can only be row, col, both!")
    if (!(colorStyle %in% c("s1", "s2", "s3", "s4")))
        stop("colorStyle can only be s1, s2, s3, s4!")

    if (missing(rowGroup)) rowGroup <- rep(1, nrow(data1))
    #if (missing(colGroup)) colGroup <- rep(1, ncol(data1))
    
    if (length(rowGroup) != nrow(data1)) stop("rowGroup has different length from the row of data1!")
    #if (length(colGroup) != min(ncol(data1), ncol(data2))) stop("colGroup has different length from the column of data1 or data2!")

    if (rowGroupColor && !missing(rowGroupColor.choice ) && (length(rowGroupColor.choice)!= length(unique(rowGroup)))) stop("rowGroupColor.choice has different length from rowGroup!")
    
    if (  (!is.null(orderRowGroup)) && (length(orderRowGroup)!=length(unique(rowGroup))) )
        stop("orderRowGroup is not NULL and has group number different from those in rowGroup!")

    if (any(!(orderRowGroup %in% rowGroup))) stop("orderRowGroup must match the groups in rowGroup!")

    if (!is.null(rowNameGroupColor) &&(length(rowNameGroupColor) != length(unique(rowGroup)))) stop("rowNameGroupColor must have same length as the group number in rowGroup!")
    if (!(legend.pos %in% c("top", "middle", "bottom")))
        stop("legend.pos can only be top, middle, bottom!")
    if (legend.percent > 1 || legend.percent < 0)
        stop("legend.percent can only between 0 and 1!")

    matSepRow <- pairheatmapENV$matSepRow
    

## scale data
    resultScale <- scaleData(data1, data2, scale=scale)
    data1 <- resultScale$data1
    data2 <- resultScale$data2


## setup configure
    configure.list <- getConfigure(data1, data2, matSepRow, matDist
                            ,rowGroupColor
                            )
    barwidthv <- configure.list$barwidthv
    barwidth <- configure.list$barwidth
    barheight <- configure.list$barheight

    x <-  configure.list$x
    y <-  configure.list$y

    data1.x <- x[1:ncol(data1)]
    data1.y <- y[1:nrow(data1)]
    barSep.x <- x[ncol(data1)+1]

    data2.x <- x[(ncol(data1)+2):(ncol(data1)+1+ncol(data2))]

    totalwidth <- ncol(data1) + ncol(data2) + matDist*matSepRow
    totalheight <- nrow(data1)

## cluster the matrix here
    mlist <- rankMatrix(data1, data2, rowGroup,
                        clusterMethod, clusterMembers,
                        orderRowGroup, clusterRow, clusterCol, clusterColTogether)

    rowInd <- mlist$rowInd
    colInd.m1 <- mlist$colInd.m1
    colInd.m2 <- mlist$colInd.m2
    rowGroup <- rowGroup[rowInd]
    
    data1 <- data1[rowInd, colInd.m1]
    data2 <- data2[rowInd, colInd.m2]

## get color
    colorSegment <- colorMat(style=colorStyle)

    datav <- unlist(cbind(data1, data2))
    colorMatrix <- matchColor(data1, colorSegment, min(datav), max(datav))
    colorMatrix2 <- matchColor(data2, colorSegment, min(datav), max(datav))

## draw figure
    # set para for grid
    deviceName <- names(dev.cur())
    if (!(deviceName %in% c("pdf", "bmp", "jpeg", "png", "tiff", "postscript"))) grid.newpage()

    rname <- rownames(data1)
    rind <- which.max(nchar(rname))
    layout.colnameWidth <- stringWidth(rname[rind])

    treeHeight <- unit(1.5, "cm")
    colNameHeight <- unit(10, "mm")
    dataHeight <- unit(1, "npc") -  treeHeight - colNameHeight

    rowNameWidth <- unit(20, "mm")
    legendWidth <- unit(15, "mm")

    dataWidth <- unit(1, "npc") - rowNameWidth - legendWidth - treeHeight

    layout.height <- unit.c(treeHeight, dataHeight, colNameHeight )
    layout.width <- unit.c(treeHeight, dataWidth, rowNameWidth, legendWidth)

    layout.rowNo <- 3
    layout.colNo <- 4

    cheat.layout <- grid.layout(layout.rowNo, layout.colNo,
                                    heights=layout.height,
                                    widths=layout.width)
    pushViewport(viewport(layout = cheat.layout

                ))

    ## draw dendrogram - vertical (2, 1)
    if (dendrogram == "row"  || dendrogram == "both" )
    {
        pushViewport(viewport(layout.pos.row=2, layout.pos.col=1))
        draw_dendrogram.mod(data1, clusterMethod, horizontal = FALSE, rowGroup, xgroup1=data1.y)
        upViewport()
    }

    #  draw dendrogram - hroizon(1, 2)
    if (dendrogram == "col"  || dendrogram == "both" )
    {
        pushViewport(viewport(layout.pos.row=1, layout.pos.col=2))

        if (clusterColTogether)
        {
            draw_dendrogram.mod(t(data1), clusterMethod, horizontal = TRUE, rowGroup, xgroup1=data1.x, xgroup2=data2.x)
        }
        else
        {
            draw_dendrogram.mod(t(data1), clusterMethod, horizontal = TRUE, rowGroup, xgroup1=data1.x, xgroup2=data2.x, t(data2))
        }

        upViewport()
    }

    ## data matrix (2,2)
    pushViewport(viewport(layout.pos.row=2, layout.pos.col=2))

    totalcol <- ncol(data1) + ncol(data2) + matSepRow
    totalrow <- nrow(data1)

    for (i in 1:totalcol)
    {

        if (i <= ncol(data1))
        {
            grid.rect(x = x[i], y = y[1:totalrow],
                      width=barwidthv[i], height=barheight,
                              just=c("left", "bottom")
                              ,gp=gpar(fill=colorMatrix[,i], col=matrixBorderCol))
        }
        else if ( (i <= ncol(data1) + matSepRow) && (i > ncol(data1)) )
        {
            grid.rect(x = x[i], y = y[1:totalrow],
                      width=barwidthv[i], height=barheight,
                              just=c("left", "bottom")
                              ,gp=gpar(fill="white", col="white"))
        }
        else
        {
            color.ind <- ncol(data2) + i - totalcol
            grid.rect(x = x[i], y = y[1:totalrow],
                      width=barwidthv[i], height=barheight,
                              just=c("left", "bottom")
                              ,gp=gpar(fill=colorMatrix2[,color.ind], col=matrixBorderCol))
        }

    }

    ## global param
    annoWidth <- convertWidth(unit(5, "mm"), "npc", valueOnly=TRUE)
    annoHeight <- convertHeight(unit(5, "mm"), "npc", valueOnly=TRUE)

    real.1.width <- 1 - annoWidth
    real.1.height <- 1 - annoHeight

    if (!missing(rowGroup))
    {
        ## order the group from top to bottom
        rowGroup <- rowGroup
        row.group <- rle(rowGroup)$lengths
        cum.group <- cumsum(row.group)
        ybar <- y[cum.group+1]  # note 0 based y value, 1 based index value
        ybar <- c(0, ybar[-length(ybar)], y[length(y)])
        yheight <- diff(ybar)

        ## add anntation color bar
        if (rowGroupColor)
        {
            if (!missing(rowGroupColor.choice)) rowBarCol <- rowGroupColor.choice
            else rowBarCol <- colorMatrix[cum.group[1:length(row.group)],1]
            
            for (i in 1:length(row.group))
            {
              grid.rect(x = x[length(x)], y = ybar[i],
                          width=barwidthv[length(barwidthv)], height=barheight*row.group[i],
                          just=c("left", "bottom")
                                  ,gp=gpar(fill=rowBarCol[i],
                                           col="grey"))
             }

        }
    }

    ## add group border
    if (rowGroupColor)
      xborder <- x[length(x)]
    else xborder <- 1

    ## add group
    ybar.nozero <- ybar[-c(1, length(ybar))]

    if (groupBorder=="line")
    {
       for (i in ybar.nozero)
       {
          grid.lines(x=c(0, xborder), y=c(i,i),
                    gp=gpar(lwd=groupBorder.lwd, col=groupBorder.col)
                    )
       }

    }
    else if (groupBorder=="rect")
    {

       for (i in 1:(length(ybar)-1) )
       {
          grid.rect(x = 0, y = ybar[i],
                        width=xborder, height=yheight[i],
                        just=c("left", "bottom")
                                ,gp=gpar(lwd=groupBorder.lwd, col=groupBorder.col))

       }
    }
    
    ## Group selection
    #group 1
    if (length(groupBorder.selectList)>0)
    {

        for (i in 1:length(groupBorder.selectList$xgroup.start))
        {
                  grid.rect(x = x[groupBorder.selectList$xgroup.start[i]],
                            y = y[groupBorder.selectList$ygroup.start[i]],
                            width=x[groupBorder.selectList$xgroup.end[i]+1]-x[groupBorder.selectList$xgroup.start[i]],
                            height=y[groupBorder.selectList$ygroup.end[i]+1]-y[groupBorder.selectList$ygroup.start[i]],
                            just=c("left", "bottom")
                                        ,gp=gpar(lwd=groupBorder.lwd, col=groupBorder.col))

        #group 2
                  grid.rect(x = x[groupBorder.selectList$xgroup.start[i] + ncol(data1) + 1],
                            y = y[groupBorder.selectList$ygroup.start[i]],
                            width=x[groupBorder.selectList$xgroup.end[i]+1]-x[groupBorder.selectList$xgroup.start[i]],
                            height=y[groupBorder.selectList$ygroup.end[i]+1]-y[groupBorder.selectList$ygroup.start[i]],
                            just=c("left", "bottom")
                                        ,gp=gpar(lwd=groupBorder.lwd, col=groupBorder.col))
        }
    }


    upViewport()

    ## draw rowname (2,3)
    pushViewport(viewport(layout.pos.row=2, layout.pos.col=3))

    # row color
    if (!is.null(rowNameGroupColor)) rowname.col <- rep(rowNameGroupColor, row.group) #rowname.col <- rep(colorMatrix[cum.group,1], row.group)
    else rowname.col <- rep(rowNameColor, length(rownames(data1)))

    drawRownames(txt=rownames(data1),
                 x=rep(0, length(y))+0.05, y=y+1/totalheight/2,
                 fontsize=rowNameFontSize,
                 col=rowname.col
                 )

    upViewport()

    ## draw colname (3,2)
    pushViewport(viewport(layout.pos.row=3, layout.pos.col=2))

    row.x <- c(data1.x, data2.x) + 1/totalwidth/2
    row.y <- rep(1, length(row.x))*0.95

    col.name <- c(colnames(data1), colnames(data2))

    # col color
    #if (colNameGroupColor) colname.col <- rep(rep(colorMatrix[1, cum.group.y], col.group),2)
    #else
    colname.col <- rep(colNameColor, length(col.name))

    drawColnames(txt=col.name, x=row.x, y=row.y, fontsize=colNameFontSize, col=colname.col)

    upViewport()

    ## draw legend (2, 4)
    pushViewport(viewport(layout.pos.row=2, layout.pos.col=4))


    tmin <- min(as.vector(data1), as.vector(data2))
    tmax <- max(as.vector(data1), as.vector(data2))

    legendLabel <- round(seq(tmin, tmax, length.out= pairheatmapENV$totalLegendLabel),0)

## TODO
#   - check should do here
#  if (length(colorSegment) > 5) 

    drawLegend(colorSegment, legendLabel, legend.pos, legend.percent, legend.fontsize)
    upViewport()

    popViewport(0)

    return(invisible(list(rowInd=rowInd, colInd1=colInd.m1, colInd2=colInd.m2)))
}


                    