makeElipseCoords <- function(x0 = 0, y0 = 0, b = 1, a = 1, alpha = 0, pct.range = c(0,1), len = 50){
  rad.range <- 2 * pi * pct.range
  theta <- seq(rad.range[1], rad.range[2], length.out=len)
  x <- x0 + a * cos(theta) * cos(alpha) - b * sin(theta) * sin(alpha)
  y <- y0 + a * cos(theta) * sin(alpha) + b * sin(theta) * cos(alpha)
  return(tab2df(cbind(x,y)))
}


plotClock <- function(hour, minute, x0 = 0, y0 = 0, r = 1){  
  
  circleXY <- makeElipseCoords(x0 = x0, y0 = y0, b = 1.1*r, a = 1.1*r, alpha = 0, 
                               pct.range = c(0,1), len = 50)
  quarHourTickMarksXY <- makeElipseCoords(x0 = x0, y0 = y0, b = 1.05*r, a = 1.05*r, alpha = (pi/2), 
                               pct.range = c((12*4-1)/(12*4),0), len = 12*4)
  hourLabelsXY <- makeElipseCoords(x0 = x0, y0 = y0, b = .9*r, a = .9*r, alpha = (pi/2), 
                               pct.range = c(11/12,0), len = 12)

  polygon(circleXY)
  text(hourLabelsXY[,1],hourLabelsXY[,2],seq(1,12), cex=.5)
  text(quarHourTickMarksXY[,1],quarHourTickMarksXY[,2],".")

  minuteV <- minute/60
  minuteVXY <- makeElipseCoords(x0 = x0, y0 = y0, b = r, a = r, alpha = 0, 
  				               pct.range =  (.25 - rep(minuteV,2)), len = 1)
  segments(x0,y0,minuteVXY$x[1],minuteVXY$y[1])

  hourV <- hour/12
  hourVXY <- makeElipseCoords(x0 = x0, y0 = y0, b = .7*r, a =.7*r, alpha = 0, 
  				               pct.range = (.25 - rep(hourV,2)), len = 1)
  segments(x0,y0,hourVXY$x,hourVXY$y)  

}

## function to grab par's usr param for use in subsequence plots
usr2lims <- function(adj=.04){
  
  par.usr <- par('usr')
  xlims <- par.usr[c(1,2)]
  ylims <- par.usr[c(3,4)]

  adj <- adj - adj^2 *2 # this simplifies the math below
  
  xlims <- xlims + c(adj,-adj) *  diff(xlims)
  ylims <- ylims + c(adj,-adj) *  diff(ylims)

  return(list(x=xlims, y=ylims))
}




### a better pie function with origin positions ###
pies <- function(x, show.labels = FALSE, show.slice.labels = FALSE, color.table = NULL, 
		radii = rep(2,length(x)), x0=NULL, y0=NULL, 
		edges = 200,  clockwise = FALSE, 
                init.angle = if (clockwise) 90 else 0, density = NULL, angle = 45, 
                border = NULL, lty = NULL,  
                other.color='gray', na.color='white', ...) 
{
  
  if(!par()$new){
    plot(x0, y0, pch='', ...)
    par(new=TRUE)
  }
  if(class(x)!='list')
    stop("x must be a list")
  
  if(length(x) != length(x0) | length(x0) != length(y0))
    stop(paste("x0 and y0 lengths (",length(x0),',',length(y0),") must match length of x (",length(x),")", sep=''))
  
  if(length(radii) < length(x))
    radii <- rep(radii, length.out=length(x))
  
  ## calculate the char size to pie radius conversions
  cx <- .25 * par('cxy')[1]
  cy <- .19 * par('cxy')[2]
  # old -> * (par('csi')/par('pin')[2]) * diff(ylim) * .2 # inches to coords scaling
  
  radii <- radii  
  xlim <- usr2lims()$x
  ylim <- usr2lims()$y
  y2x.asp <- diff(xlim)/diff(ylim)

  pie.labels <- names(x)

  if (is.null(color.table)) {
    unique.labels <- unique(unlist(lapply(x,names)))
    color.table <- rainbow(length(unique.labels))
    names(color.table) <- unique.labels
  }
  
  ## loop through the list of pie tables
  for(j in seq(along=x)){
    X <- x[[j]]
    data.labels <- names(X)

    if(j != 1)
      par(new=TRUE)
    if (!is.numeric(X) || any(is.na(X) | X < 0)) 
        stop("'x' values must be non-missing positive.")


    if(length(X) == 0){
      warning(paste(names(x)[[j]], 'has zero length vector'))

    }else{
      
      ## generate a slice fraction vector
      X <- c(0, cumsum(X)/sum(X))
      names(X) <- data.labels  #re-label it
      
      dx <- diff(X)
      nx <- length(dx)
      plot.new()
      pin <- par("pin")
      if(all(xlim == c(-1, 1)) && all(ylim == c(-1,1))){
        if (pin[1] > pin[2]) 
          xlim <- (pin[1]/pin[2]) * xlim
        else ylim <- (pin[2]/pin[1]) * ylim
      }
      
      plot.window(xlim, ylim, "") #, asp = 1)

      ## change to gray all of the X names without colors in the color table
      nolgnd <- names(X)[!names(X) %in% names(color.table) ]
      color.table <- c(color.table, nv(rep(other.color,length(nolgnd)),nolgnd))
      col <- color.table[names(X)]
      col[names(col)=='NA'] <- na.color 

      if(length(border)> 1){
        if(length(border) != length(x))
          stop('length of border doesnt equal length of x')
        this.brdr <- border[j]
      }else{
        this.brdr <- border
      }
      
      lty <- rep(lty, length.out = nx)
      angle <- rep(angle, length.out = nx)
      density <- rep(density, length.out = nx)
      twopi <- if (clockwise) 
        -2 * pi
      else 2 * pi
      ## function to turn theta into xy coordinates
      t2xy <- function(t) {
        t2p <- twopi * t + init.angle * pi/180
        list(x = radii[j] *cx * cos(t2p), y =  radii[j] * cy * sin(t2p)) #y
      }
      
      ## loop through each slice of the pie
      for (i in 1:nx) {
        lab <- as.character(names(X)[i])
        nx <- as.character(names(X)[i])
        n <- max(2, floor(edges * dx[i]))
        P <- t2xy(seq.int(X[i], X[i + 1], length.out = n))
        polygon(c(x0[j]+ P$x, x0[j]), c(y0[j] +P$y, y0[j]), density = density[i], angle = angle[i], 
                border = this.brdr, col = col[lab], lty = lty[i], ...)
        P <- t2xy(mean(X[i + 0:1]))
        
        if (!is.na(lab) && nzchar(lab)) {
          if(show.slice.labels){
            lines(x0[j] +c(1, 1.05) * P$x, y0[j] +c(1, 1.05) * P$y)
            text(x0[j] +1.1 * P$x, y0[j] + 1.1 * P$y, lab, xpd = TRUE, 
                 adj = ifelse(P$x < 0, 1, 0))
          }
        }
      }
      if(show.labels)
        text(x0[j],y0[j] + radii[j]+.2, pie.labels[j])
      
      invisible(NULL)
    }    
  }
}



#annulus()  #donut or ring plot


vennMatrix <- function(l){
  ## much of the code in this function was inspired by parts of Yongmin Sun's doVennDiagram function
  if(is.null(names(l)))
    stop("The list 'l' must have named elements")
  l.all <- unique(do.call(c,lapply(l, as.character)))
  l.mat <- matrix(0, nrow = length(l.all), ncol = length(l))
  colnames(l.mat) <- names(l)
  rownames(l.mat) <- l.all

  for(i in 1:length(l.all)) 
    for(nm in names(l))
      l.mat[i,nm] <- l.all[i] %in% l[[nm]]
  return(l.mat)
}



textplot <- function(..., x=1, y=1){
  plot(x, y, pch='', bty='n',xaxt='n',yaxt='n', xlab='', ylab='')
  text(x, y, ...)
}




mvlabs <- function(df, n=nrow(df), x='x', y='y', l='lab', cols=colors()[grep("dark",colors())], ...){

  for(i in (1:n)+1){

    ## identify point to move
    idx <- identify(x=df[,x], y=df[,y], labels=df[,l],n=1)
    print(df[idx,])
    ## locate new location to move to
    locs <- locator(n=1)

    ## move the point  (refresh the plot?)
    df[idx, x] <- locs$x[1]
    df[idx, y] <- locs$y[1]
    #points(x=df[idx, x], y=df[idx, y], col=colors()[i])
    text(x=df[idx, x], y=df[idx, y], labels=as.character(df[idx,l]), col=cols[i], ...)

  }
  return(df)
}



labsegs <- function(x0, y0, x1, y1, buf=.3, ...){

  a <- x1 - x0
  b <- y1 - y0
  c0 <- sqrt(a^2 + b^2)
  theta <- atan(b/a)
  theta[a<0] <- theta[a<0] + pi
  
  c1 <- c0 - buf
  
  if(any(c1<0))
    stop('buffer size too large or annotations too close')
  
  a1 <- c1*cos(theta)
  b1 <- c1*sin(theta)

  x1 <- x0 + a1
  y1 <- y0 + b1

  segments(x0,y0,x1,y1,...)
}

  
hyperplot <- function(x, y=NULL, annout=1:length(x), name='hyperplot.imagemap', w=72*8, h=72*6, link='internal', browse=TRUE, cex=1, ...){
 
  ## generate output paths
  img.path <- paste(name,'.png',sep='')
  html.path <- paste(name,'.html',sep='')

  ## create image
  png(img.path, width=w, height=h)

  if(class(x)=='data.frame')
    stop('for data.frame input: x and y vectors should be in the annout table with x & y used to specify column names')
    	
  ## plot
  if(is.null(y)){  
    # x as an object with a native 'plot' method and xyc return data.frame
    xyc <- plot(x, ...)
    if(ncol(xyc)<2 | ncol(xyc)>3)
      stop('plot(x) must return a rownamed dataframe or matrix of 2 or 3 columns: x,y and optionally cex (in that order)')
    x <- xyc[,1]
    y <- xyc[,2]
    if(ncol(xyc)==3)
      cex <- xyc[,3]
    else
      cex <- 1
    idx <- rownames(xyc)
  }else{          

    if(all(is.character(c(x,y))) & all(sapply(list(x,y), length)==1)){ 
      # x and y as column names
      if(all(c(x,y)%in% names(annout)))
        idx <- rownames(annout)
      else
        stop('x and y must exist in annout')
      
      x <- annout[,x]
      y <- annout[,y]
      idx <- rownames(annout)
    }else{    
      # x and y as named numeric vectors and annout as index to x (and y)

      if(is.data.frame(annout) | is.matrix(annout)){
        annout.idx  <- row.names(annout)
      }else{ # is.vector
        annout.idx <- annout  
        link='none'
      }
      
      if(is.null(names(x))){ # annout as a character vector or named dataframe
        if(is.numeric(annout.idx)){
          idx <- as.character(1:length(x))  # annout index is converted below in df check
        }else{
          if(length(annout.idx) == length(x))
            idx <- annout.idx 
          else
            stop("length of name index supplied form annout must be same as length as x")
        }    
      }else{  # annout as index into names of x
        idx <- names(x)
        if(is.numeric(annout.idx))
          annout <- names(x)[annout.idx]
      }
    }
    
    if(!is.numeric(x) | !is.numeric(y) | length(x) != length(y))
      stop('x and y must be numeric vectors of equal length')
    
    plot(x, y, ...)

  }
   
  ## build x, y coordinates, names & point sizes 
  map <- data.frame(idx=idx, x=x, y=y, cex=cex, row.names=idx)
  
  ## grab figure & plot dimentions
  mai <- par('mai')*72 #margins
  pin <- par('pin')*72 #plot dim
  usr <- par('usr')
  
  usr.xd <-diff(usr[1:2])
  usr.yd <-diff(usr[3:4])

  pin.xd <-pin[1]
  pin.yd <-pin[2]

  ## save image
  dev.off()
  
  ## determine outlier points to annotate
  if(is.data.frame(annout))
    annout$nm <- rownames(annout)
  else
    annout <- data.frame(nm=as.character(annout))

  if(length(link) != 1 | !is.character(link))
    stop("'link' must be a character string specifying if points are linked 'internal'ly \
           or which column of annot to use as external hyperlinks")

  if(!link %in% c('none','internal',colnames(annout)))
    stop("'link' must be either 'none', 'internal' or a column name of annout")
  
  if('out' %in% names(annout))
    annout <- annout[annout$out, ]

  ## subset coords by only the desired outlier list  
  map <- subset(map, idx %in% annout$nm)
  
  if(nrow(map) < 1)
    stop('no points to map because none of the annotations matched?')
  
  if(link == 'none')
    map$href <- ''
  else    
    if(link == 'internal')
      map$href <- paste('#',sub(' ','_',map$idx), sep='')
    else
      map <- nerge(list(map=map, href=nv(as.character(annout[,link]), annout$nm)))

  for(clmn in names(annout)[sapply(annout, is.character)])
    if(substr(annout[1,clmn],1,4) %in% c('http','www.'))
      annout[,clmn] <- paste("<a target='_blank' href='",annout[,clmn],"'>",annout[,clmn],"</a>",sep='')
  
  ## translate user coordinate space into image output coordinates
  map$x <-  (map$x - usr[1])/usr.xd * pin.xd + mai[2]
  map$y <- (-map$y + usr[4])/usr.yd * pin.yd + mai[3]
  map$r <- map$cex * 3

  ## write HTML map
  sink(html.path)

  cat('<html>\n')
  cat(paste('<img src="',img.path,'" width="',w,'" height="',h,'" usemap="#',name,'"/>', sep=''),'\n')
  cat(paste('<map name="',name,'">',sep=''))    
  with(map, cat(paste('<area shape="circle" coords="',x,',',y,',',r,'" href="',href,'" title="',idx,'"/>\n',sep='')),'\n')
  cat('</map>')

  if(link == 'internal' & class(annout)=='data.frame'){
    cat('<br><h5>Annotations</h5>')
    cat('<table border=1>','\n')
    cat(paste('<tr><th></th>', paste('<th>', names(annout),'</th>',collapse=''), '</tr>\n'))
    for(i in 1:nrow(annout))
      cat(paste('<tr><td><a name="',sub(' ','_',annout$nm[i]),'"></a></td>',
                paste('<td>',annout[i,],'</td>',collapse=''), '</tr>\n',sep=''))
    cat('</table>','\n')
  }
  cat('</html>')
  cat(rep('<br>',60)) ## so hyperlink jumps don't run into the bottom of the browser window')
  sink()

  ## open browser 
  if(browse)
    browseURL(html.path)
}



.rect3venn <- function(xt){

    lim <- max(xt)
    nms <- names(dimnames(xt))

    t0 <- xt[1,1,1]
    t1 <- xt[2,1,1]
    t2 <- xt[1,2,1]
    t3 <- xt[1,1,2]
    t12 <- xt[2,2,1]
    t23 <- xt[1,2,2]
    t13 <- xt[2,1,2]
    t123 <- xt[2,2,2]

    s123 <- sqrt(t123)
    s12 <- (t12-t123)/s123
    s23 <- (t23-t123)/s123
    s13 <- (t13-t123)/s123

    s1 <- t1/(s13+s123+s12)
    s2 <- t2/(s123+s23)
    s3 <- t3/(s123+s23)

    s2 <- sqrt(t2)
    #s3 <- sqrt(t3)

    plot.new()
    plot.window(c(-lim,lim),c(-lim,lim),xaxs="i",yaxs="i")

    #par(cex=0.9)

    rect(-s13    ,0        ,s123+s12,s1  , col="pink")
    #rect(0       ,-s23     ,s2     ,s123, col="springgreen")
    rect(0       ,-s2+s123 ,s2      ,s123, col="springgreen")
    rect(s123-s3 ,-s23     ,s123    ,s123, col="lightblue")
    #rect(-s3+s123,-s3+s123 ,s123    ,s123, col="lightblue")

    rect(s123,0   ,s123+s12,s123, col="yellow")
    rect(-s13,0   ,0       ,s123, col="violet")
    rect(0   ,-s23,s123    ,0   , col="cyan")

    rect(0   ,0   ,s123    ,s123,col="white")

    text(s123/2           ,s123/2, t123    )
    text(-s23/2           ,s123/2, t13-t123)
    text(s123+(s12/2)     ,s123/2, t12-t123)
    text(s123/2           ,-s23/2, t23-t123)

    text((-s13+s123+s12)/2,s1-3      ,paste(nms[1],"n=",t1))

    text(s2/2             ,-s2+s123+3,paste(nms[2],"n=",t1))

    text(-(s3-s123)/2     ,-s23+3    ,paste(nms[3],"n=",t3))
}

heatmatrix <- function(x, values=TRUE, clp=c('bottom','top'), rlp=c('left','right'), xadj=.02, yadj=.3, ylab.cntr=FALSE, cex=1, cex.axis=1, ...){

  image(1L:ncol(x), 1L:nrow(x), t(x[nrow(x):1,]),  xaxt='n', yaxt='n', ... )

  if(!is.null(rownames(x))){
    clp2par <- nv(c(1,2),clp)
    clp2xadj <- nv(c(-1,1),clp) * xadj
    clp2adj <- nv(c(1,0),clp); if(ylab.cntr) clp2adj <- nv(rep(.5,2),clp)
    clp <- match.arg(clp)
    text(x=par("usr")[clp2par[clp]] +clp2xadj[clp], y=nrow(x):1,
adj=clp2adj[clp],
         labels = rownames(x), xpd = TRUE, cex=cex.axis)
  }
  if(!is.null(colnames(x))){
    rlp2par <- nv(c(3,4),rlp)
    rlp2yadj <- nv(c(-1,1),rlp) * yadj
    rlp <- match.arg(rlp)
    text(x=1:ncol(x), y=par("usr")[rlp2par[rlp]] + rlp2yadj[rlp], adj=.5,
         labels = colnames(x), xpd = TRUE, cex=cex.axis)
  }
  if(values)
     text(col(x),row(x), round(x[nrow(x):1,],2), cex=cex)

}

