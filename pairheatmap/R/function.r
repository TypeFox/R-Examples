############################################
##  R package pairheatmap
##  function: cfunction.r
############################################

## scale data
scaleData <- function(data1, data2, scale="none")
{
    if (scale=="row")
    {
        ctwo <- cbind(data1, data2)
        crowMean <- apply(ctwo, 1, mean, na.rm=TRUE)
        crowSD <- apply(ctwo, 1, sd, na.rm=TRUE)
        cdone <- (ctwo-crowMean)/crowSD
        data1 <- cdone[,1:ncol(data1)]
        data2 <- cdone[,(ncol(data1)+1):ncol(cdone)]
    }
    else if (scale=="rowsep")
    {
        crowMean <- apply(data1, 1, mean, na.rm=TRUE)
        crowSD <- apply(data1, 1, sd, na.rm=TRUE)
        data1 <- (data1-crowMean)/crowSD

        crowMean <- apply(data2, 1, mean, na.rm=TRUE)
        crowSD <- apply(data2, 1, sd, na.rm=TRUE)
        data2 <- (data2-crowMean)/crowSD
    }
    else if (scale=="col")
    {
        crowMean <- apply(data1, 2, mean, na.rm=TRUE)
        crowSD <- apply(data1, 2, sd, na.rm=TRUE)
        data1 <- (data1-crowMean)/crowSD
        
        crowMean <- apply(data2, 2, mean, na.rm=TRUE)
        crowSD <- apply(data2, 2, sd, na.rm=TRUE)
        data2 <- (data2-crowMean)/crowSD
    }

    
    return(list(data1=data1, data2=data2))
}

## generate colorMatrix
colorMat <- function(style="s1")
{
    colorSegment <- switch(style,
                           s1 = colorRampPalette(c("blue", "red"))(100),
                           s2 = colorRampPalette(c("green", "red"))(100),
                           s3 = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF", "#E0F3F8", "#91BFDB", "#4575B4")))(100),
                           s4 = colorRampPalette(c("white", "black"))(100)
                           )
    return(colorSegment)
}

## match data with color, return a colorMatrix with that data format
matchColor <- function(dataMatrix, colorSegment, dmin, dmax)
{
    # split data
    datav <- unlist(dataMatrix)
    data.boundary <- seq(dmin, dmax, length.out=length(colorSegment))
    data.color <- colorSegment[as.numeric(cut(datav, data.boundary, include.lowest =TRUE))]

    colorMatrix <- matrix(data.color, ncol=ncol(dataMatrix))
    return(colorMatrix)
}

getConfigure <- function(data1=data1, data2=data2, 
                    matSepRow=matSepRow, matDist=matDist, 
                    rowGroupColor=TRUE, colGroupColor=FALSE)
{
    totalcol <- ncol(data1) + ncol(data2) + matSepRow  
    totalrow <- nrow(data1)
    datacol <-  ncol(data1) + ncol(data2) + matSepRow
    
    totalwidth <- ncol(data1) + ncol(data2) + matDist*matSepRow # 1: paramter as line
    totalheight <- nrow(data1)
    
    annoWidth <- convertX(unit(5, "mm"), "npc", valueOnly=TRUE)
    
    if (rowGroupColor)
      real.1.width <- 1 - annoWidth
    else
      real.1.width <- 1 
    
    if (colGroupColor)  
      real.1.height <- 1 - annoWidth
    else
      real.1.height <- 1
      
    barheight <- unit(real.1.height/totalheight, "native")

    barwidth <- unit(real.1.width/totalwidth, "native")
    sepwidth <- unit(real.1.width/totalwidth*matDist*matSepRow, "native")
    
    if (rowGroupColor)
      barwidthv <- c(rep(barwidth, ncol(data1)), sepwidth, 
                     rep(barwidth, ncol(data2)), annoWidth 
                     )
    else
      barwidthv <- c(rep(barwidth, ncol(data1)), sepwidth, 
                     rep(barwidth, ncol(data2))                
                     )
    x <- c(0, cumsum(barwidthv))
    x <- round(x[-length(x)],3)
    
    y <- (0:(totalrow))*real.1.height/totalheight
    
    return(list(barwidthv=barwidthv, barheight=barheight, x=x, y=y, barwidth=barwidth))
}

## sub function
distMatrix <- function(hc, xgroup1)
{
	h = hc$height / max(hc$height)/1.1
	m = hc$merge
	o = hc$order
	n = length(o)

	m[m > 0] = n + m[m > 0]
	m[m < 0] = abs(m[m < 0])

	dist = matrix(0, nrow = 2 * n - 1, ncol = 2, dimnames = list(NULL, c("x", "y")))

	## calcuate x
	if (is.null(xgroup1))
      dist[1:n, 1] = 1 / n / 2 + (1 / n) * (match(1:n, o) - 1)
  else
      dist[1:n, 1] = diff(xgroup1)[1]/ 2 + diff(xgroup1)[1] * (match(1:n, o) - 1)

	for(i in 1:nrow(m))
  {
		dist[n + i, 1] = (dist[m[i, 1], 1] + dist[m[i, 2], 1]) / 2
		dist[n + i, 2] = h[i]
	}
	return(list(dist=dist, h=h, m=m))
}

draw_connection = function(x1, x2, y1, y2, y)
{
		grid.lines(x = c(x1, x1), y = c(y1, y))
		grid.lines(x = c(x2, x2), y = c(y2, y))
		grid.lines(x = c(x1, x2), y = c(y, y))
}

## draw dengdram within xgroup2 specify range
drawDistanceDen <- function(hc2, xgroup2)
{
        xtmp <- xgroup2 - xgroup2[1]
        distList <- distMatrix(hc2, xtmp)
        dist2 <- distList$dist
        h2 <- distList$h
        m2 <- distList$m
        dist2[1:nrow(dist2), 1] <- dist2[1:nrow(dist2), 1] + xgroup2[1]

    		for(i in 1:nrow(m2))
        {
            draw_connection(dist2[m2[i, 1], 1], dist2[m2[i, 2], 1], dist2[m2[i, 1], 2], dist2[m2[i, 2], 2], h2[i])
        }

}

## method and some codes are from R package: pheatmap
draw_dendrogram.mod <- function(data1, clusterMethod, horizontal = TRUE, rowGroup, xgroup1=NULL, xgroup2=NULL, data2=NULL)
{
  hc <- hclust(dist(data1), clusterMethod)
  if (!is.null(data2)) hc2 <- hclust(dist(data2), clusterMethod)
  distList <- distMatrix(hc, xgroup1)
  dist <- distList$dist
  h <- distList$h
  m <- distList$m

	if(horizontal)
  {

    for(i in 1:nrow(m))
    {
      draw_connection(dist[m[i, 1], 1], dist[m[i, 2], 1], dist[m[i, 1], 2], dist[m[i, 2], 2], h[i])
    }

    if (!is.null(data2) && !is.null(xgroup2))
    {
        xtmp <- xgroup2 - xgroup2[1]
        distList <- distMatrix(hc2, xtmp)
        dist2 <- distList$dist
        h2 <- distList$h
        m2 <- distList$m
        dist2[1:nrow(dist2), 1] <- dist2[1:nrow(dist2), 1] + xgroup2[1]
        
    		for(i in 1:nrow(m2))
        {
            draw_connection(dist2[m2[i, 1], 1], dist2[m2[i, 2], 1], dist2[m2[i, 1], 2], dist2[m2[i, 2], 2], h2[i])
        }

    }

	}
	else
  {  ## edit code
  		gr = rectGrob()

      pushViewport(viewport(height = unit(1, "grobwidth", gr), width = unit(1, "grobheight", gr), angle = 90))
      if (length(unique(rowGroup))==1 )
      {
      		for(i in 1:nrow(m))
          {
      			draw_connection(dist[m[i, 1], 1], dist[m[i, 2], 1], dist[m[i, 1], 2], dist[m[i, 2], 2], h[i])
      		}
      }
      else
      {
          dengroup.unique <- unique(rowGroup)
          
          for (i in 1:length(dengroup.unique))
          {
              select.ind <- which(rowGroup==dengroup.unique[i])
              if (length(select.ind) > 2)
              {
                  den.xgroup <- xgroup1[select.ind]
                  hcgroup <- hclust(dist(data1[select.ind,]), clusterMethod)
                  drawDistanceDen(hcgroup, den.xgroup)
              }
          }
        
      }
  		upViewport()
	}
}


## draw legend
drawLegend <- function(colorSegment, legendLabel, legend.pos = "top",
                       legend.percent=0.5,
                       legend.fontsize=6
                       )
{
    ## total legend length is 0.5 npc
    x <- unit(1, "mm")+unit(1, "picas")
    y <- 1:length(colorSegment)*legend.percent/length(colorSegment)
    height <- diff(y)[1]
    
    # calculate middle
    yfinal <- switch(legend.pos,
            top=y+(1-legend.percent),
            middle = y + 0.5 - y[round(length(y)/2)],
            bottom = y
            )
    legGrob <- rectGrob(name="leg.rect", x=x, y=yfinal,
                           width=unit(1, "picas"),
                           height=height, 
                           just=c("centre", "bottom"),
                           gp=gpar(fill=colorSegment, col=colorSegment))    
    xtxt <- x + unit(1, "picas") + unit(2, "mm") 
    ytxt <- seq(min(yfinal), max(yfinal), length.out=pairheatmapENV$totalLegendLabel)
    
    txtGrob <- textGrob(name="leg.txt", label=legendLabel, 
                x=xtxt, y=ytxt, gp=gpar(fontsize=legend.fontsize))
    grid.draw(legGrob)
    grid.draw(txtGrob)    
}


drawGroupBar <- function(x, y, height=height, width=width)
{
    sideg <- rectGrob(name="sideGrob", x = x, y = y,
                  width=width, height=height,
                          just=c("left", "bottom")
                          ,gp=gpar(fill="pink", col="blue"))
    grid.draw(sideg)
    
}


drawRownames <- function(txt, x, y, col, ...)
{   ## TODO, add gp for label
    #labelGrob <- textGrob(name="labelGrob", txt, x=x, y=y,                          
    #grid.draw(labelGrob)
    
    for (i in 1:length(txt))
    {
        grid.text(txt[i], x=x[i], y=y[i],
                  just=c("left", "centre"), gp=gpar(col=col[i],...)) 
    }
    
}

drawColnames <- function(txt, x, y, col,...)
{   ## TODO, add gp for label
    #labelGrob <- textGrob(name="labelGrob", txt, x=x, y=y,                          
    #grid.draw(labelGrob)

    for (i in 1:length(txt))
    {
        grid.text(txt[i], x=x[i]
                  , y=y[i],
                  just=c("left", "centre"), 
                  rot=270, 
                  gp=gpar(col=col[i],...)) 
    }
    
}