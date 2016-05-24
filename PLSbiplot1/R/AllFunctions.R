
rm(list=ls(all=TRUE))  #clear R memory

# Date: 18 - October - 2014

# R CODES USED IN THIS DISSERTATION


#A.1  Additional R codes for running PCA and PLS biplots
# Additional functions for running PCA and PLS biplots.
##1. Axis calibration 
calibrate.axis = function(j, unscaled.X, means, sd, axes.rows, ax.which, 
                          ax.tickvec, ax.orthogxvec, ax.orthogyvec,
                          ax.oblique, p) 
{
  ax.num = ax.which[j]
  tick = ax.tickvec[j]
  ax.direction = axes.rows[j,]  
  r = ncol(axes.rows)
  ax.orthog = rbind(ax.orthogxvec, ax.orthogyvec)
  if (nrow(ax.orthog)<r) 
    ax.orthog = rbind(ax.orthog, 0)
  if (nrow(axes.rows)>1) 
    phi.vec = diag(1/diag(axes.rows %*% t(axes.rows))) %*% axes.rows %*% ax.orthog[,j]
  else 
    phi.vec = (1/(axes.rows%*%t(axes.rows))) %*% axes.rows %*% ax.orthog[,j]
  
  number.points = 100
  std.ax.tick.label = pretty(unscaled.X[, ax.num], n=tick)
  std.range = c(min(std.ax.tick.label), max(std.ax.tick.label))
  std.ax.tick.label.min = std.ax.tick.label - (std.range[2] - std.range[1])
  std.ax.tick.label.max = std.ax.tick.label + (std.range[2] - std.range[1])
  std.ax.tick.label = c(std.ax.tick.label, std.ax.tick.label.min, std.ax.tick.label.max)
  interval = (std.ax.tick.label - means[ax.num]) / sd[ax.num]
  axis.vals = seq(from=min(interval), to=max(interval), length=number.points)
  axis.vals = sort(unique(c(axis.vals, interval)))
  number.points = length(axis.vals)
  
  # axis.points is a matrix
  # -- col1, ... col(r): co-ordinates to trace the axis
  # -- col (r+1): standard values in terms of unscaled.X
  # -- col (r+2): indicator whether co-ordinate should be calibrated
  
  axis.points = matrix(0, nrow=number.points, ncol=r)
  for (i in 1:r)
    axis.points[, i] = ax.orthog[i,ax.num] + (axis.vals-phi.vec[ax.num]) * ax.direction[i]
  if (!is.null(ax.oblique)) 
    for (i in 1:r)
      axis.points[,i] = axis.vals*ax.direction[i] - ((ax.oblique[ax.num]-means[ax.num])/sd[ax.num]) * 
    ax.direction[i] + (((ax.oblique-means)/sd) %*% axes.rows[,i])/p
  axis.points = cbind(axis.points, axis.vals*sd[ax.num]+means[ax.num], 0)
  for (i in 1:number.points) 
    if (any(zapsmall(axis.points[i, r+1]-std.ax.tick.label) == 0)) 
      axis.points[i, r+2] = 1
  axis.points
}

##2. Biplot sample control  
biplot.sample.control = function(J, col = c(1:8,3:1), pch=3, cex=1, 
                                 label=F, label.cex=0.75,
                                 label.side="bottom", alpha=1)
{
  while (length(col)<J) col = c(col, col)
  col = as.vector(col[1:J])
  
  while (length(pch)<J) pch = c(pch, pch)
  pch = as.vector(pch[1:J])
  
  while (length(cex)<J) cex = c(cex, cex)
  cex = as.vector(cex[1:J])
  
  while (length(label)<J) label = c(label, label)
  label = as.vector(label[1:J])
  
  while (length(label.cex)<J) label.cex = c(label.cex, label.cex)
  label.cex = as.vector(label.cex[1:J])
  
  while (length(label.side)<J) label.side = c(label.side, label.side)
  label.side = as.vector(label.side[1:J])
  
  while (length(alpha)<J) alpha = c(alpha, alpha)
  alpha = as.vector(alpha[1:J])
  
  list(col=col, pch=pch, cex=cex, label=label, label.cex=label.cex, 
       label.side=label.side, alpha=alpha)
}

##3. Biplot axis control
biplot.ax.control = function(p, X.names, which=1:p, type = "prediction",
                             col="grey", lwd=1, lty=1, label="Orthog", 
                             label.col=col, label.cex=0.75, label.dist=0,
                             ticks=5, tick.col=col, tick.size=1,
                             tick.label=T, tick.label.col=tick.col, 
                             tick.label.cex=0.6, tick.label.side="left", 
                             tick.label.offset=0.5, tick.label.pos=1, 
                             predict.col=col, predict.lwd=lwd, 
                             predict.lty=lty, ax.names=X.names, 
                             rotate=NULL, orthogx=0, orthogy=0, 
                             oblique=NULL) 
{
  if (!all(is.numeric(which)))  
    which = match(which, X.names, nomatch=0)
  which = which[which <= p]
  which = which[which > 0]
  ax.num = length(which)
  
  if (type != "prediction" & type != "interpolation") 
    stop ("Incorrect type of biplot axes specified")
  
  while (length(col)<ax.num) col = c(col, col)
  col = as.vector(col[1:ax.num])
  
  while (length(lwd)<ax.num) lwd = c(lwd, lwd)
  lwd = as.vector(lwd[1:ax.num])
  
  while(length(lty)<ax.num) lty = c(lty, lty)
  lty = as.vector(lty[1:ax.num])
  
  if (label != "Orthog" & label != "Hor" & label != "Paral") 
    stop ("Incorrect specification of axis label direction")
  
  while (length(label.col)<ax.num) label.col = c(label.col, label.col)
  label.col = as.vector(label.col[1:ax.num])
  
  while (length(label.cex)<ax.num) label.cex = c(label.cex, label.cex)
  label.cex = as.vector(label.cex[1:ax.num])
  
  while (length(label.dist)<ax.num) label.dist = c(label.dist, label.dist)
  label.dist = as.vector(label.dist[1:ax.num])
  
  while (length(ticks)<ax.num) ticks = c(ticks, ticks)
  ticks = as.vector(ticks[1:ax.num])   
  
  while (length(tick.col)<ax.num) tick.col = c(tick.col, tick.col)
  tick.col = as.vector(tick.col[1:ax.num])
  
  while (length(tick.size)<ax.num) tick.size = c(tick.size, tick.size)
  tick.size = as.vector(tick.size[1:ax.num]) 
  
  while (length(tick.label)<ax.num) tick.label = c(tick.label, tick.label)
  tick.label = as.vector(tick.label[1:ax.num])
  
  while (length(tick.label.col)<ax.num) tick.label.col = c(tick.label.col, tick.label.col)
  tick.label.col = as.vector(tick.label.col[1:ax.num])
  
  while (length(tick.label.cex)<ax.num) tick.label.cex = c(tick.label.cex, tick.label.cex)
  tick.label.cex = as.vector(tick.label.cex[1:ax.num])
  
  while (length(tick.label.side)<ax.num) tick.label.side = c(tick.label.side, tick.label.side)
  tick.label.side = as.vector(tick.label.side[1:ax.num])
  
  while (length(tick.label.offset)<ax.num) tick.label.offset = c(tick.label.offset, tick.label.offset)
  tick.label.offset = as.vector(tick.label.offset[1:ax.num])
  
  while (length(tick.label.pos)<ax.num) tick.label.pos = c(tick.label.pos, tick.label.pos)
  tick.label.pos = as.vector(tick.label.pos[1:ax.num])
  
  while (length(predict.col)<ax.num) predict.col = c(predict.col, predict.col)
  predict.col = as.vector(predict.col[1:ax.num])
  
  while (length(predict.lwd)<ax.num) predict.lwd = c(predict.lwd, predict.lwd)
  predict.lwd = as.vector(predict.lwd[1:ax.num])
  
  while (length(predict.lty)<ax.num) predict.lty = c(predict.lty, predict.lty)
  predict.lty = as.vector(predict.lty[1:ax.num])
  
  while (length(ax.names)<p) ax.names = c(ax.names, "")
  ax.names = as.vector(ax.names[1:ax.num])
  
  if (!is.null(oblique)) 
    if (length(oblique) != p) 
      stop ("For oblique translations values must be specified for each variable")
  
  while (length(orthogx)<p) orthogx = c(orthogx, orthogx)
  orthogx = as.vector(orthogx[1:p])
  
  while (length(orthogy)<p) orthogy = c(orthogy, orthogy)
  orthogy = as.vector(orthogy[1:p])
  
  list(which=which, type=type, col=col, lwd=lwd, lty=lty, label=label, 
       label.col=label.col, label.cex=label.cex, label.dist=label.dist, 
       ticks=ticks, tick.col=tick.col, tick.size=tick.size,	   
       tick.label=tick.label, tick.label.col=tick.label.col, 
       tick.label.cex=tick.label.cex, tick.label.side=tick.label.side, 
       tick.label.offset=tick.label.offset, tick.label.pos=tick.label.pos, 
       predict.col=predict.col, predict.lty=predict.lty,
       predict.lwd=predict.lwd, names=ax.names, rotate=rotate,
       orthogx=orthogx, orthogy=orthogy, oblique=oblique)
}

##4. Draw biplot
draw.biplot = function(Z, G=matrix(1,nrow=nrow(Z),ncol=1), classes=1, 
                       Z.means=NULL, z.axes=NULL, z.bags=NULL, 
                       z.ellipse=NULL, Z.new=NULL, Z.density=NULL,
                       sample.style, mean.style=NULL, ax.style=NULL, 
                       bag.style=NULL, ellipse.style=NULL, 
                       new.sample.style=NULL, density.style=NULL, 
                       predict.samples=NULL, predict.means=NULL, 
                       Title=NULL, VV=NULL, exp.factor=1.2, ...)
{
  .samples.plot = function(Z, G, classes, sample.style) 
  {
    x.vals = Z[, 1]
    y.vals = Z[, 2]
    invals = x.vals < usr[2] & x.vals > usr[1] & y.vals < usr[4] & y.vals >  
      usr[3]
    Z = Z[invals, ]
    for (j in 1:length(classes)) 
    {
      class.num = classes[j]
      Z.class = Z[G[,class.num] == 1, , drop=FALSE]
      text.pos = match(sample.style$label.side[j], c("bottom","left","top","right"))
      
      if (sample.style$label[j]) 
        text(Z.class[,1], Z.class[,2], labels=dimnames(Z.class)[[1]], 
             cex=sample.style$label.cex[j], pos=text.pos)
      for (i in 1:nrow(Z.class)) 
        points(x=Z.class[i,1], y=Z.class[i,2], pch=sample.style$pch[j], 
               col=sample.style$col[j], cex=sample.style$cex[j])
    }    
  }
  # ---------------------------------------------------------
  
  .predict.func = function(p.point, coef, col, lty, lwd)
  {
    if (is.na(coef[2]))  #horizontal axis
      lines(c(p.point[1],coef[1]), rep(p.point[2],2), col=col, lwd=lwd, lty=lty)
    else 
      if (coef[2]==0) #vertical axis
        lines(rep(p.point[1],2), p.point[2:3], col=col, lwd=lwd, lty=lty)
    else  #axis with diagonal slope
    {
      intercept.projection = p.point[2] + p.point[1]/coef[2]
      project.on.x = (intercept.projection - coef[1]) / (coef[2] + 1/coef[2])
      project.on.y = coef[1] + coef[2]*project.on.x
      lines(c(p.point[1],project.on.x), c(p.point[2],project.on.y), col=col, lwd=lwd, lty=lty)
    }
  }
  .marker.label.cm = function(x, y, grad, marker.val, expand = 1, col, 
                              label.on.off, side, pos, offset, label.col, 
                              cex) 
  {
    uin = par("pin")/c(usr[2] - usr[1], usr[4] - usr[3])
    mm = 1/(uin[1] * 25.4)
    d = expand * mm
    
    if (grad=="v") 
    {  
      lines(rep(x,2), c(y-d,y+d), col=col)
      if (label.on.off==1) 
        text(x, y-d, label=marker.val, pos=pos, offset=offset, col=label.col, cex=cex)
    }
    if (grad=="h") 
    { 
      lines(c(x-d,x+d), rep(y,2), col=col)
      if (label.on.off==1) 
        if (side=="right") 
          text(x+d, y, label=marker.val, pos=pos, offset=offset, col=label.col, cex=cex)
      else 
        text(x-d, y, label=marker.val, pos=pos, offset=offset, col=label.col, cex=cex)
    }
    if (is.numeric(grad))
    {
      b = d * sqrt(1/(1 + grad * grad))
      a = b * grad
      lines(c(x - b, x + b), c(y - a, y + a), col=col)
      if (label.on.off==1) 
        if (side=="right") 
          text(x+b, y+a, label=marker.val, pos=pos, offset=offset, col=label.col, cex=cex)
      else
        text(x-b, y-a, label=marker.val, pos=pos, offset=offset, col=label.col, cex=cex)
    }
  }
  .marker.func = function(vec, coef, col, tick.size, side, pos, offset, label.col, cex)
  {
    x = vec[1]
    y = vec[2]
    marker.val = vec[3]
    label.on.off = vec[4]
    
    if (is.na(coef[2]))
      .marker.label.cm(x, y, grad="h", marker.val, expand=tick.size, 
                       col=col, label.on.off=label.on.off, side=side, 
                       pos=pos, offset=offset, label.col=label.col, cex=cex)
    else 
      if (coef[2]==0)
        .marker.label.cm(x, y, grad="v", marker.val, expand=tick.size, 
                         col=col, label.on.off=label.on.off, side=side, 
                         pos=pos, offset=offset, label.col=label.col, 
                         cex=cex)
    else
      .marker.label.cm(x, y, grad=-1/coef[2], marker.val, 
                       expand=tick.size, col=col, 
                       label.on.off=label.on.off, side=side, 
                       pos=pos, offset=offset, label.col=label.col, cex=cex)
  }
  .lin.axes.plot = function(z.axes, ax.style, predict.mat)
  {
    for (i in 1:length(ax.style$which))
    {
      ax.num = ax.style$which[i]
      marker.mat = z.axes[[ax.num]][z.axes[[ax.num]][, 4]==1, 1:3]
      marker.mat = marker.mat[rev(order(marker.mat[,3])),]
      x.vals = marker.mat[, 1]
      y.vals = marker.mat[, 2]
      #draw axis
      lin.coef = coefficients(lm(y.vals~x.vals))
      if (is.na(lin.coef[2])) 
        abline(v=x.vals, col=ax.style$col[i], lwd=ax.style$lwd[i], lty=ax.style$lty[i])
      else
        abline(coef=lin.coef, col=ax.style$col[i], lwd=ax.style$lwd[i], lty=ax.style$lty[i])
      #label axis; marker.mat ordered in decreasing marker values; label at 
      #positive side of axis
      if (ax.style$label == "Hor") 
      { 
        par(las=1)
        adjust = c(0.5, 1, 0.5, 0) 
      }
      if (ax.style$label == "Orthog") 
      { 
        par(las=2)
        adjust = c(1, 1, 0, 0) 
      }	
      if (ax.style$label == "Paral") 
      { 
        par(las=0)
        adjust = c(0.5, 0.5, 0.5, 0.5) 
      }
      h = nrow(marker.mat)
      #vertical axis
      if (is.na(lin.coef[2]))
      {
        if (y.vals[1] < y.vals[h]) #label bottom
          mtext(text=ax.style$names[i], side=1, line=ax.style$label.dist[i], adj=adjust[1], 
                at=x.vals[1], col=ax.style$label.col[i], cex=ax.style$label.cex[i])
        else  #label top
          mtext(text=ax.style$names[i], side=3, line=ax.style$label.dist[i], adj=adjust[3], 
                at=y.vals[1], col=ax.style$label.col[i], cex=ax.style$label.cex[i])
      }
      else
      {	
        y1.ster = lin.coef[2] * usr[1] + lin.coef[1]
        y2.ster = lin.coef[2] * usr[2] + lin.coef[1]
        x1.ster = (usr[3] - lin.coef[1]) / lin.coef[2]
        x2.ster = (usr[4] - lin.coef[1]) / lin.coef[2]
        
        #horizontal axis  
        if (lin.coef[2]==0)
        {
          if (x.vals[1] < x.vals[h]) #label left
            mtext(text=ax.style$names[i], side=2, line=ax.style$label.dist[i], adj=adjust[2],
                  at=y.vals[1], col=ax.style$label.col[i], cex=ax.style$label.cex[i])
          else  #label right
            mtext(text=ax.style$names[i], side=4, line=ax.style$label.dist[i], adj=adjust[4], 
                  at=y.vals[1], col=ax.style$label.col[i], cex=ax.style$label.cex[i])
        }
        #positive gradient  
        if (lin.coef[2]>0)
        {
          if (x.vals[1] < x.vals[h]) #label left or bottom
            if (y1.ster <= usr[4] & y1.ster >= usr[3]) #left
              mtext(text=ax.style$names[i], side=2, line=ax.style$label.dist[i], adj=adjust[2], 
                    at=y1.ster, col=ax.style$label.col[i], cex=ax.style$label.cex[i])
          else  #bottom
            mtext(text=ax.style$names[i], side=1, line=ax.style$label.dist[i], 
                  adj=adjust[1], at=x1.ster, col=ax.style$label.col[i], 
                  cex=ax.style$label.cex[i])
          else  #label right or top
            if (y2.ster <= usr[4] & y2.ster >= usr[3]) #right
              mtext(text=ax.style$names[i], side=4, line=ax.style$label.dist[i], 
                    adj=adjust[4], at=y2.ster, col=ax.style$label.col[i], 
                    cex=ax.style$label.cex[i])
          else  #top
            mtext(text=ax.style$names[i], side=3, line=ax.style$label.dist[i], 
                  adj=adjust[3], at=x2.ster, col=ax.style$label.col[i], 
                  cex=ax.style$label.cex[i])
        }
        #negative gradient
        if (lin.coef[2]<0)
        {
          if (x.vals[1] < x.vals[h]) #label left or top
            if (y1.ster <= usr[4] & y1.ster >= usr[3]) #left
              mtext(text=ax.style$names[i], side=2, line=ax.style$label.dist[i], adj=adjust[2], 
                    at=y1.ster, col=ax.style$label.col[i], cex=ax.style$label.cex[i])
          else  #top
            mtext(text=ax.style$names[i], side=3, line=ax.style$label.dist[i], adj=adjust[3], 
                  at=x2.ster, col=ax.style$label.col[i], cex=ax.style$label.cex[i])
          else  #label right or bottom
            if (y2.ster <= usr[4] & y2.ster >= usr[3]) #right
              mtext(text=ax.style$names[i], side=4, line=ax.style$label.dist[i], 
                    adj=adjust[4], at=y2.ster, col=ax.style$label.col[i], 
                    cex=ax.style$label.cex[i])
          else  #top
            mtext(text=ax.style$names[i], side=1, line=ax.style$label.dist[i], 
                  adj=adjust[1], at=x1.ster, col=ax.style$label.col[i], 
                  cex=ax.style$label.cex[i])
        }
      }	
      
      #axis tick marks		
      invals = x.vals<usr[2] & x.vals>usr[1] & y.vals<usr[4] & y.vals>usr[3] 
      std.markers = zapsmall(marker.mat[invals, 3])
      x.vals = x.vals[invals]
      y.vals = y.vals[invals]
      if (ax.style$tick.label[i]) 
        label.on.off = rep(1,sum(invals))
      else 
        rep(0,sum(invals))
      if (!ax.style$tick.label[i]) 
        label.on.off[c(1,length(label.on.off))] = 1
      apply(cbind(x.vals,y.vals,std.markers,label.on.off), 1, .marker.func, coef=lin.coef, 
            col=ax.style$tick.col[i], tick.size=ax.style$tick.size[i], 
            side=ax.style$tick.label.side[i], pos=ax.style$tick.label.pos[i], 
            offset=ax.style$tick.label.offset[i], label.col=ax.style$tick.label.col[i], 
            cex=ax.style$tick.label.cex[i])
      
      #predictions on axis
      #for horizontal axis, info is needed on the axis y-values
      if (!is.null(predict.mat))
        apply(cbind(predict.mat,y.vals[1]), 1, .predict.func, coef=lin.coef, 
              col=ax.style$predict.col[i], lty=ax.style$predict.lty[i], lwd=ax.style$predict.lwd[i])
    }
  }
  # ---------------------------------------------------------
  
  .bags.plot = function(z.bags, bag.style)
  {
    for (i in 1:length(z.bags))
    {
      mat = cbind(unlist(z.bags[[i]][1]), unlist(z.bags[[i]][2]))
      mat = rbind(mat, mat[1,])
      lines(mat, col=bag.style$col[i], lty=bag.style$lty[i], lwd=bag.style$lwd[i])
      if (bag.style$Tukey.median[i])
        points(unlist(z.bags[[i]][3]), col=bag.style$col[i], pch=bag.style$pch, cex=bag.style$cex)
    }
  }
  # ---------------------------------------------------------
  
  .ellipse.plot = function(z.ellipse, ellipse.style)
  {
    for (i in 1:length(z.ellipse))
    {
      #mat = cbind (unlist(z.bags[[i]][1]), unlist(z.bags[[i]][2]))
      #mat = rbind (mat, mat[1,])
      lines(z.ellipse[[i]], col=ellipse.style$col[i], lty=ellipse.style$lty[i], lwd=ellipse.style$lwd[i])
    }
  }
  # ---------------------------------------------------------
  
  .new.samples.plot = function (Z.new, new.sample.style)
  {
    points(Z.new[,1], Z.new[,2], pch=new.sample.style$pch, col=new.sample.style$col, 
           cex=new.sample.style$cex)
    pos.vec = rep(0,nrow(Z.new))
    pos.vec = match(new.sample.style$label.side, c("bottom","left","top","right"))
    if (any(new.sample.style$label)) 
      text(Z.new[new.sample.style$label,1], Z.new[new.sample.style$label,2],     
           labels=dimnames(Z.new)[[1]][new.sample.style$label],
           cex=new.sample.style$label.cex[new.sample.style$label], pos=pos.vec[new.sample.style$label])
  }
  # ---------------------------------------------------------
  
  .class.means.plot = function(Z.means, mean.style)
  {
    points(Z.means[,1], Z.means[,2], pch=mean.style$pch, col=mean.style$col, cex=mean.style$cex)
    pos.vec=rep(0,nrow(Z.means))
    pos.vec=match(mean.style$label.side, c("bottom","left","top","right"))
    if (any(mean.style$label)) 
      text(Z.means[mean.style$label,1], Z.means[mean.style$label,2], 
           labels=dimnames(Z.means)[[1]][mean.style$label], 
           cex=mean.style$label.cex[mean.style$label], pos=pos.vec[mean.style$label])
  }
  # ---------------------------------------------------------
  
  .density.plot = function(Z.density, density.style)
  {
    levels.rect = pretty(range(Z.density$z), n = density.style$cuts)
    col.use = colorRampPalette(density.style$col)
    col.use = col.use(length(levels.rect) - 1)
    image(Z.density, breaks=levels.rect, col=col.use, add=TRUE)
    if (density.style$contours) 
      contour(Z.density, levels=levels.rect, col=density.style$contour.col, add = TRUE)  
    list(levels.rect, col.use)	
  }
  
  .density.legend = function(levels.rect, col.use)
  {
    par(pty="m", mar=density.style$legend.mar)	
    plot(range(levels.rect), y=1:2, ylim=c(10,100), xaxs="i", yaxs="i", xlab="", ylab="", xaxt="n",
         yaxt="n", type="n", asp=1, frame.plot=FALSE)
    rect(xleft=levels.rect[-length(levels.rect)], ybottom=10, xright=levels.rect[-1], ytop=50, 
         col=col.use, border=FALSE)
    axis(side=1, at=pretty(levels.rect,n=8), labels=pretty(levels.rect,n=8), line=0, 
         cex.axis=density.style$cex, mgp=density.style$mgp, tcl=density.style$tcl, las=0)
  }
  
  # ---------------------------------------------------------
  old.par = par(no.readonly=TRUE)
  on.exit(par(old.par))
  par(pty = "s", ...)
  
  if (!is.null(Z.density))
    layout(mat=matrix(1:2,ncol=1), heights=density.style$layout.heights)
  
  #setting up plot	
  plot(Z[,1]*exp.factor, Z[,2]*exp.factor, xlim=range(Z[,1]*exp.factor), ylim=range(Z[,2]*exp.factor), 
       xaxt="n", yaxt="n", xlab="", ylab="", type="n", xaxs="i", yaxs="i", asp=1)
  usr = par("usr")
  
  if (!is.null(predict.samples)) 
    predict.mat = Z[predict.samples,,drop=F] 
  else 
    predict.mat = NULL
  if (!is.null(predict.means)) 
    predict.mat = rbind(predict.mat, Z.means[predict.means,,drop=F])
  
  #density plots
  if (!is.null(Z.density))
    density.out = .density.plot(Z.density, density.style)	
  else
    density.out = NULL	
  
  #plotting biplot axes
  if (!is.null(z.axes))
    .lin.axes.plot(z.axes, ax.style, predict.mat)
  
  #plotting of sample points
  if (length(classes)>0) 
    .samples.plot(Z, G, classes, sample.style)
  
  #plotting class means
  if (length(mean.style$which)>0)
    .class.means.plot(Z.means, mean.style)
  
  #plotting new samples
  if (!is.null(Z.new)) 
    .new.samples.plot(Z.new, new.sample.style)
  
  #plotting alpha bags	  
  if (length(z.bags)>0)
    .bags.plot(z.bags, bag.style)
  
  #plotting kappa ellipses	  
  if (length(z.ellipse)>0)
    .ellipse.plot(z.ellipse, ellipse.style)
  
  #if (!is.null(Title)) 
  #  main(Title)
  
  #density plots
  if (!is.null(density.out))
    .density.legend(density.out[[1]], density.out[[2]])	
  
  #create vectors
  arrows(VV[, 1]*0, VV[, 2]*0, VV[, 1], VV[, 2], length=0.09, col="green")        
}

##5. Draw triangles in area biplot
draw.biplot.with.triangles = function(Z, G=matrix(1,nrow=nrow(Z),ncol=1), 
                                      classes=1, Z.means=NULL, z.axes=NULL, 
                                      z.bags=NULL, z.ellipse=NULL, Z.new=NULL, 
                                      Z.density=NULL, sample.style, mean.style=NULL, 
                                      ax.style=NULL, bag.style=NULL, ellipse.style=NULL, 
                                      Title=NULL, new.sample.style=NULL, BB=NULL,
                                      density.style=NULL, QQ=NULL, predict.samples=NULL, 
                                      base.tri=NULL, predict.means=NULL, exp.factor=1.2,
                                      bi.value=NULL, VV=NULL, ...)
{
  .samples.plot = function(Z, G, classes, sample.style) 
  {
    x.vals = Z[, 1]
    y.vals = Z[, 2]
    invals = x.vals < usr[2] & x.vals > usr[1] & y.vals < usr[4] & y.vals >  
      usr[3]
    Z = Z[invals, ]
    for (j in 1:length(classes)) 
    {
      class.num = classes[j]
      Z.class = Z[G[,class.num] == 1, , drop=FALSE]
      text.pos = match(sample.style$label.side[j], c("bottom","left","top","right"))
      
      if (sample.style$label[j]) 
        text(Z.class[,1], Z.class[,2], labels=dimnames(Z.class)[[1]], 
             cex=sample.style$label.cex[j], pos=text.pos)
      for (i in 1:nrow(Z.class)) 
        points(x=Z.class[i,1], y=Z.class[i,2], pch=sample.style$pch[j], 
               col=sample.style$col[j], cex=sample.style$cex[j])
    }    
  }
  # ---------------------------------------------------------
  
  .predict.func = function(p.point, coef, col, lty, lwd)
  {
    if (is.na(coef[2]))  #horizontal axis
      lines(c(p.point[1],coef[1]), rep(p.point[2],2), col=col, lwd=lwd, lty=lty)
    else 
      if (coef[2]==0) #vertical axis
        lines(rep(p.point[1],2), p.point[2:3], col=col, lwd=lwd, lty=lty)
    else  #axis with diagonal slope
    {
      intercept.projection = p.point[2] + p.point[1]/coef[2]
      project.on.x = (intercept.projection - coef[1]) / (coef[2] + 1/coef[2])
      project.on.y = coef[1] + coef[2]*project.on.x
      lines(c(p.point[1],project.on.x), c(p.point[2],project.on.y), col=col, lwd=lwd, lty=lty)
    }
  }
  .marker.label.cm = function(x, y, grad, marker.val, expand = 1, col, 
                              label.on.off, side, pos, offset, label.col, 
                              cex) 
  {
    uin = par("pin")/c(usr[2] - usr[1], usr[4] - usr[3])
    mm = 1/(uin[1] * 25.4)
    d = expand * mm
    
    if (grad=="v") 
    {  
      lines(rep(x,2), c(y-d,y+d), col=col)
      if (label.on.off==1) 
        text(x, y-d, label=marker.val, pos=pos, offset=offset, col=label.col, cex=cex)
    }
    if (grad=="h") 
    { 
      lines(c(x-d,x+d), rep(y,2), col=col)
      if (label.on.off==1) 
        if (side=="right") 
          text(x+d, y, label=marker.val, pos=pos, offset=offset, col=label.col, cex=cex)
      else 
        text(x-d, y, label=marker.val, pos=pos, offset=offset, col=label.col, cex=cex)
    }
    if (is.numeric(grad))
    {
      b = d * sqrt(1/(1 + grad * grad))
      a = b * grad
      lines(c(x - b, x + b), c(y - a, y + a), col=col)
      if (label.on.off==1) 
        if (side=="right") 
          text(x+b, y+a, label=marker.val, pos=pos, offset=offset, col=label.col, cex=cex)
      else
        text(x-b, y-a, label=marker.val, pos=pos, offset=offset, col=label.col, cex=cex)
    }
  }
  .marker.func = function(vec, coef, col, tick.size, side, pos, offset, label.col, cex)
  {
    x = vec[1]
    y = vec[2]
    marker.val = vec[3]
    label.on.off = vec[4]
    
    if (is.na(coef[2]))
      .marker.label.cm(x, y, grad="h", marker.val, expand=tick.size, 
                       col=col, label.on.off=label.on.off, side=side, 
                       pos=pos, offset=offset, label.col=label.col, cex=cex)
    else 
      if (coef[2]==0)
        .marker.label.cm(x, y, grad="v", marker.val, expand=tick.size, 
                         col=col, label.on.off=label.on.off, side=side, 
                         pos=pos, offset=offset, label.col=label.col, 
                         cex=cex)
    else
      .marker.label.cm(x, y, grad=-1/coef[2], marker.val, 
                       expand=tick.size, col=col, 
                       label.on.off=label.on.off, side=side, 
                       pos=pos, offset=offset, label.col=label.col, cex=cex)
  }
  .lin.axes.plot = function(z.axes, ax.style, predict.mat)
  {
    for (i in 1:length(ax.style$which))
    {
      ax.num = ax.style$which[i]
      marker.mat = z.axes[[ax.num]][z.axes[[ax.num]][, 4]==1, 1:3]
      marker.mat = marker.mat[rev(order(marker.mat[,3])),]
      x.vals = marker.mat[, 1]
      y.vals = marker.mat[, 2]
      #draw axis
      lin.coef = coefficients(lm(y.vals~x.vals))
      if (is.na(lin.coef[2])) 
        abline(v=x.vals, col=ax.style$col[i], lwd=ax.style$lwd[i], lty=ax.style$lty[i])
      else
        abline(coef=lin.coef, col=ax.style$col[i], lwd=ax.style$lwd[i], lty=ax.style$lty[i])
      #label axis; marker.mat ordered in decreasing marker values; label at 
      #positive side of axis
      if (ax.style$label == "Hor") 
      { 
        par(las=1)
        adjust = c(0.5, 1, 0.5, 0) 
      }
      if (ax.style$label == "Orthog") 
      { 
        par(las=2)
        adjust = c(1, 1, 0, 0) 
      }	
      if (ax.style$label == "Paral") 
      { 
        par(las=0)
        adjust = c(0.5, 0.5, 0.5, 0.5) 
      }
      h = nrow(marker.mat)
      #vertical axis
      if (is.na(lin.coef[2]))
      {
        if (y.vals[1] < y.vals[h]) #label bottom
          mtext(text=ax.style$names[i], side=1, line=ax.style$label.dist[i], adj=adjust[1], 
                at=x.vals[1], col=ax.style$label.col[i], cex=ax.style$label.cex[i])
        else  #label top
          mtext(text=ax.style$names[i], side=3, line=ax.style$label.dist[i], adj=adjust[3], 
                at=y.vals[1], col=ax.style$label.col[i], cex=ax.style$label.cex[i])
      }
      else
      {	
        y1.ster = lin.coef[2] * usr[1] + lin.coef[1]
        y2.ster = lin.coef[2] * usr[2] + lin.coef[1]
        x1.ster = (usr[3] - lin.coef[1]) / lin.coef[2]
        x2.ster = (usr[4] - lin.coef[1]) / lin.coef[2]
        
        #horizontal axis  
        if (lin.coef[2]==0)
        {
          if (x.vals[1] < x.vals[h]) #label left
            mtext(text=ax.style$names[i], side=2, line=ax.style$label.dist[i], adj=adjust[2],
                  at=y.vals[1], col=ax.style$label.col[i], cex=ax.style$label.cex[i])
          else  #label right
            mtext(text=ax.style$names[i], side=4, line=ax.style$label.dist[i], adj=adjust[4], 
                  at=y.vals[1], col=ax.style$label.col[i], cex=ax.style$label.cex[i])
        }
        #positive gradient  
        if (lin.coef[2]>0)
        {
          if (x.vals[1] < x.vals[h]) #label left or bottom
            if (y1.ster <= usr[4] & y1.ster >= usr[3]) #left
              mtext(text=ax.style$names[i], side=2, line=ax.style$label.dist[i], adj=adjust[2], 
                    at=y1.ster, col=ax.style$label.col[i], cex=ax.style$label.cex[i])
          else  #bottom
            mtext(text=ax.style$names[i], side=1, line=ax.style$label.dist[i], 
                  adj=adjust[1], at=x1.ster, col=ax.style$label.col[i], 
                  cex=ax.style$label.cex[i])
          else  #label right or top
            if (y2.ster <= usr[4] & y2.ster >= usr[3]) #right
              mtext(text=ax.style$names[i], side=4, line=ax.style$label.dist[i], 
                    adj=adjust[4], at=y2.ster, col=ax.style$label.col[i], 
                    cex=ax.style$label.cex[i])
          else  #top
            mtext(text=ax.style$names[i], side=3, line=ax.style$label.dist[i], 
                  adj=adjust[3], at=x2.ster, col=ax.style$label.col[i], 
                  cex=ax.style$label.cex[i])
        }
        #negative gradient
        if (lin.coef[2]<0)
        {
          if (x.vals[1] < x.vals[h]) #label left or top
            if (y1.ster <= usr[4] & y1.ster >= usr[3]) #left
              mtext(text=ax.style$names[i], side=2, line=ax.style$label.dist[i], adj=adjust[2], 
                    at=y1.ster, col=ax.style$label.col[i], cex=ax.style$label.cex[i])
          else  #top
            mtext(text=ax.style$names[i], side=3, line=ax.style$label.dist[i], adj=adjust[3], 
                  at=x2.ster, col=ax.style$label.col[i], cex=ax.style$label.cex[i])
          else  #label right or bottom
            if (y2.ster <= usr[4] & y2.ster >= usr[3]) #right
              mtext(text=ax.style$names[i], side=4, line=ax.style$label.dist[i], 
                    adj=adjust[4], at=y2.ster, col=ax.style$label.col[i], 
                    cex=ax.style$label.cex[i])
          else  #top
            mtext(text=ax.style$names[i], side=1, line=ax.style$label.dist[i], 
                  adj=adjust[1], at=x1.ster, col=ax.style$label.col[i], 
                  cex=ax.style$label.cex[i])
        }
      }	
      
      #axis tick marks		
      invals = x.vals<usr[2] & x.vals>usr[1] & y.vals<usr[4] & y.vals>usr[3] 
      std.markers = zapsmall(marker.mat[invals, 3])
      x.vals = x.vals[invals]
      y.vals = y.vals[invals]
      if (ax.style$tick.label[i]) 
        label.on.off = rep(1,sum(invals))
      else 
        rep(0,sum(invals))
      if (!ax.style$tick.label[i]) 
        label.on.off[c(1,length(label.on.off))] = 1
      apply(cbind(x.vals,y.vals,std.markers,label.on.off), 1, .marker.func, coef=lin.coef, 
            col=ax.style$tick.col[i], tick.size=ax.style$tick.size[i], 
            side=ax.style$tick.label.side[i], pos=ax.style$tick.label.pos[i], 
            offset=ax.style$tick.label.offset[i], label.col=ax.style$tick.label.col[i], 
            cex=ax.style$tick.label.cex[i])
      
      #predictions on axis
      #for horizontal axis, info is needed on the axis y-values
      if (!is.null(predict.mat))
        apply(cbind(predict.mat,y.vals[1]), 1, .predict.func, coef=lin.coef, 
              col=ax.style$predict.col[i], lty=ax.style$predict.lty[i], lwd=ax.style$predict.lwd[i])
    }
  }
  # ---------------------------------------------------------
  
  .bags.plot = function(z.bags, bag.style)
  {
    for (i in 1:length(z.bags))
    {
      mat = cbind(unlist(z.bags[[i]][1]), unlist(z.bags[[i]][2]))
      mat = rbind(mat, mat[1,])
      lines(mat, col=bag.style$col[i], lty=bag.style$lty[i], lwd=bag.style$lwd[i])
      if (bag.style$Tukey.median[i])
        points(unlist(z.bags[[i]][3]), col=bag.style$col[i], pch=bag.style$pch, cex=bag.style$cex)
    }
  }
  # ---------------------------------------------------------
  
  .ellipse.plot = function(z.ellipse, ellipse.style)
  {
    for (i in 1:length(z.ellipse))
    {
      #mat = cbind (unlist(z.bags[[i]][1]), unlist(z.bags[[i]][2]))
      #mat = rbind (mat, mat[1,])
      lines(z.ellipse[[i]], col=ellipse.style$col[i], lty=ellipse.style$lty[i], lwd=ellipse.style$lwd[i])
    }
  }
  # ---------------------------------------------------------
  
  .new.samples.plot = function (Z.new, new.sample.style)
  {
    points(Z.new[,1], Z.new[,2], pch=new.sample.style$pch, col=new.sample.style$col, 
           cex=new.sample.style$cex)
    pos.vec = rep(0,nrow(Z.new))
    pos.vec = match(new.sample.style$label.side, c("bottom","left","top","right"))
    if (any(new.sample.style$label)) 
      text(Z.new[new.sample.style$label,1], Z.new[new.sample.style$label,2],     
           labels=dimnames(Z.new)[[1]][new.sample.style$label],
           cex=new.sample.style$label.cex[new.sample.style$label], pos=pos.vec[new.sample.style$label])
  }
  # ---------------------------------------------------------
  
  .class.means.plot = function(Z.means, mean.style)
  {
    points(Z.means[,1], Z.means[,2], pch=mean.style$pch, col=mean.style$col, cex=mean.style$cex)
    pos.vec=rep(0,nrow(Z.means))
    pos.vec=match(mean.style$label.side, c("bottom","left","top","right"))
    if (any(mean.style$label)) 
      text(Z.means[mean.style$label,1], Z.means[mean.style$label,2], 
           labels=dimnames(Z.means)[[1]][mean.style$label], 
           cex=mean.style$label.cex[mean.style$label], pos=pos.vec[mean.style$label])
  }
  # ---------------------------------------------------------
  
  .density.plot = function(Z.density, density.style)
  {
    levels.rect = pretty(range(Z.density$z), n = density.style$cuts)
    col.use = colorRampPalette(density.style$col)
    col.use = col.use(length(levels.rect) - 1)
    image(Z.density, breaks=levels.rect, col=col.use, add=TRUE)
    if (density.style$contours) 
      contour(Z.density, levels=levels.rect, col=density.style$contour.col, add = TRUE)  
    list(levels.rect, col.use)	
  }
  
  .density.legend = function(levels.rect, col.use)
  {
    par(pty="m", mar=density.style$legend.mar)	
    plot(range(levels.rect), y=1:2, ylim=c(10,100), xaxs="i", yaxs="i", xlab="", ylab="", xaxt="n",
         yaxt="n", type="n", asp=1, frame.plot=FALSE)
    rect(xleft=levels.rect[-length(levels.rect)], ybottom=10, xright=levels.rect[-1], ytop=50, 
         col=col.use, border=FALSE)
    axis(side=1, at=pretty(levels.rect,n=8), labels=pretty(levels.rect,n=8), line=0, 
         cex.axis=density.style$cex, mgp=density.style$mgp, tcl=density.style$tcl, las=0)
  }
  
  # ---------------------------------------------------------
  old.par = par(no.readonly=TRUE)
  on.exit(par(old.par))
  par(pty = "s", ...)
  
  if (!is.null(Z.density))
    layout(mat=matrix(1:2,ncol=1), heights=density.style$layout.heights)
  
  #setting up plot	
  plot(Z[,1]*exp.factor, Z[,2]*exp.factor, xlim=range(Z[,1]*exp.factor), ylim=range(Z[,2]*exp.factor), 
       xaxt="n", yaxt="n", xlab="", ylab="", type="n", xaxs="i", yaxs="i", asp=1)
  usr = par("usr")
  
  if (!is.null(predict.samples)) 
    predict.mat = Z[predict.samples,,drop=F] 
  else 
    predict.mat = NULL
  if (!is.null(predict.means)) 
    predict.mat = rbind(predict.mat, Z.means[predict.means,,drop=F])
  
  #density plots
  if (!is.null(Z.density))
    density.out = .density.plot(Z.density, density.style)	
  else
    density.out = NULL	
  
  #plotting biplot axes
  if (!is.null(z.axes))
    .lin.axes.plot(z.axes, ax.style, predict.mat)
  
  #plotting of sample points
  if (length(classes)>0) 
    .samples.plot(Z, G, classes, sample.style)
  
  #plotting class means
  if (length(mean.style$which)>0)
    .class.means.plot(Z.means, mean.style)
  
  #plotting new samples
  if (!is.null(Z.new)) 
    .new.samples.plot(Z.new, new.sample.style)
  
  #plotting alpha bags	  
  if (length(z.bags)>0)
    .bags.plot(z.bags, bag.style)
  
  #plotting kappa ellipses	  
  if (length(z.ellipse)>0)
    .ellipse.plot(z.ellipse, ellipse.style)
  
  #if (!is.null(Title)) 
  # main(Title)
  
  #density plots
  if (!is.null(density.out))
    .density.legend(density.out[[1]], density.out[[2]])	
  
  #create triangles
  axis.marker = QQ[base.tri,]
  mat = NULL 
  for (i in bi.value)   
  {
    mat = rbind(mat, rbind(0, axis.marker, BB[i,], NA))
  }   
  polygon(mat, density=15, col=rainbow(nrow(BB)))  
  
  #create vectors
  arrows(VV[, 1]*0, VV[, 2]*0, VV[, 1], VV[, 2], length=0.09, col="green")   
} 

# --------------------------------------------------------------------------------------------------------------------


#A.2  Chapter  2: PCA biplots
# Example to illustrate biplot. 
# Example to illustrate Principal Component Analysis (PCA) biplot.
# One data set: 
#         - olive oil   
# Users can choose their preferred tick marker length ("ax.tickvec") per 
# variable for the PCA biplots.
# RGui(32-bit).

#' Principal Component Analysis (PCA) 
#' 
#' Takes in a samples by variables data matrix and gives the PCA parameters.
#' @param D A samples by variables data matrix    
#' @param r The number of PCA components
#' @param ... Other arguments. Currently ignored
#' @return The PCA parameters of D 
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export
mod.PCA = function(D, r, ...)
{
  #r=initial number of components    
  #D=centered D      
  options(digits=3)
  D.mean = colMeans(D) #(1xM) column means of D 	
  # or apply(D, 2, mean)
  D.hat = array(0, dim=c(nrow(D), ncol(D), r)) #predicted D-values	  
  D.svd = svd(D)
  
  #obtaining the PCA parameters (scores and loadings) in matrix form
  V = D.svd$v[,1:r]    #(Nxr) D-loadings matrix
  Z = D %*% V  #(Nxr) D-scores matrix
  
  #predicted D-values using (1:r) components	     	
  D.hat[, , r] = (Z[, 1:r] %*% t(V[, 1:r])) + rep(D.mean, each=nrow(D))      
  
  dimnames(Z) = list(rownames(D),  paste("Comp", 1:r, seq=""))
  dimnames(V) = list(colnames(D),  paste("Comp", 1:r, seq=""))
  dimnames(D.hat) = list(rownames(D), colnames(D), paste(1:r, "Comps")) 
  list(D.scores=Z, D.loadings=V, D.hat=D.hat) 
}  	

#' The Principal Component Analysis (PCA) biplot
#' 
#' Takes in a samples by variables data matrix and produces a PCA biplot.
#' @param D A samples by variables data matrix    
#' @param method the mod.PCA algorithm
#' @param ax.tickvec.D tick marker length per axis in the PCA biplot
#' @param ... Other arguments. Currently ignored
#' @return The PCA biplot of D with some parameters 
#' @examples  
#' if(require(pls))
#' data(oliveoil, package="pls")  
#' Dmat = as.matrix(oliveoil)  #(16x11) overall original data matrix      
#' dimnames(Dmat) = list(paste(c("G1","G2","G3","G4","G5","I1","I2","I3","I4","I5",
#' "S1","S2","S3","S4","S5","S6")),
#' paste(c("Acidity","Peroxide","K232","K270","DK","Yellow",
#' "Green","Brown","Glossy","Transp","Syrup")))
#' PCA.biplot(D=Dmat, method=mod.PCA, ax.tickvec.D=c(8,5,5,7,6,4,5,5,8,7,7))  
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export
PCA.biplot = function(D, method=NULL, ax.tickvec.D=NULL, ... ) 
{ 
  #D=data matrix
  options(digits=3)  
  D.mean = apply(D, 2, mean)  #column means of D or colMeans(D)
  D.std = apply(D, 2, sd)  #column standard deviations of D     
  D.scal = scale(D, center=TRUE, scale=TRUE)  #scaled D  
  
  #biplot	
  r = 2
  main = method(D.scal,r) 	 
  Zmat = main$D.scores
  Vmat = main$D.loadings
  
  #overall quality of approximation		 
  D.hat = ((Zmat %*% t(Vmat)) * rep(D.std,each=nrow(D))) + rep(D.mean,each=nrow(D))
  dimnames(D.hat) = dimnames(D)      
  overall.quality = sum(diag(t(D.hat)%*%D.hat))  /  sum(diag(t(D)%*%D))     	
  
  #axes predictivity
  axis.pred = (diag(t(D.hat)%*%D.hat)) / (diag(t(D)%*%D))	    
  names(axis.pred) = colnames(D)         	
  
  #samples predictivity
  sample.pred = (diag(D.hat%*%t(D.hat))) / (diag(D%*%t(D)))	    
  names(sample.pred) = rownames(D) 		   
  
  #biplot points	 
  Z = Zmat #(Nx2) biplot points         
  dimnames(Z) = list(rownames(D), NULL)
  sample.style = biplot.sample.control(1, label=TRUE)
  sample.style$col = "red"
  sample.style$pch = 15
  Gmat = cbind(rep(1,nrow(D)))
  dimnames(Gmat) = list(NULL, "samples")
  classes = 1 
  
  #axes direction    
  #D.hat=ZV'
  D.axes.direction = (1 / (diag(Vmat%*%t(Vmat)))) * Vmat       	   
  
  #calibration of axes	 	
  z.axes.D = lapply(1:ncol(D), calibrate.axis, D, D.mean, D.std, D.axes.direction, 1:ncol(D), 
                    ax.tickvec.D, rep(0,ncol(D)), rep(0,ncol(D)), NULL)      
  z.axes = vector("list", ncol(D))
  for (i in 1:ncol(D)) 
  {
    z.axes[[i]] = z.axes.D[[i]]
  }
  
  #biplot axes
  ax.style = biplot.ax.control(ncol(D), c(colnames(D)))     
  ax.style$col[1:ncol(D)] = "blue"
  ax.style$label.col[1:ncol(D)] = "blue"
  ax.style$tick.col[1:ncol(D)] = "blue"
  ax.style$tick.label.col[1:ncol(D)] = "blue" 
  draw.biplot(Z, G=Gmat, classes=classes, sample.style=sample.style, z.axes=z.axes, 
              ax.style=ax.style, VV=D.axes.direction) 	 	
  list(overall.quality=overall.quality, axis.pred=axis.pred, sample.pred=sample.pred, D.hat=D.hat)    
}	      

#' The Principal Component Analysis (PCA) biplot with the labels of the samples points excluded 
#' 
#' Takes in a samples by variables data matrix and produces a PCA biplot, where the labels of the samples points excluded.
#' @param D A samples by variables data matrix    
#' @param method the mod.PCA algorithm
#' @param ax.tickvec.D tick marker length per axis in the PCA biplot
#' @param ... Other arguments. Currently ignored
#' @return The PCA biplot of D with some parameters 
#' @examples 
#' if(require(pls))
#' data(oliveoil, package="pls")  
#' Dmat = as.matrix(oliveoil)  #(16x11) overall original data matrix      
#' dimnames(Dmat) = list(paste(c("G1","G2","G3","G4","G5","I1","I2","I3","I4","I5",
#' "S1","S2","S3","S4","S5","S6")),
#' paste(c("Acidity","Peroxide","K232","K270","DK","Yellow","Green","Brown",
#' "Glossy","Transp","Syrup")))
#' PCA.biplot_no.SN(D=Dmat, method=mod.PCA, ax.tickvec.D=c(8,5,5,7,6,4,5,5,8,7,7))  
#'
#' #glass data
#' if(require(chemometrics))
#' data(glass, package="chemometrics")
#' Dmat = matrix(glass,ncol=13)  
#' dimnames(Dmat) = list(1:180, paste(c("Na2O", "MgO", "Al2O3", "SiO2", 
#' "P2O5", "SO3", "Cl", "K2O", "CaO", "MnO", "Fe2O3", "BaO", "PbO")))
#' PCA.biplot_no.SN(D=Dmat, method=mod.PCA, ax.tickvec.D=rep(5,ncol(Dmat)))   
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export
PCA.biplot_no.SN = function(D, method=NULL, ax.tickvec.D=NULL, ... ) 
{ 
  #D=data matrix
  options(digits=3)  
  D.mean = apply(D, 2, mean)  #column means of D or colMeans(D)
  D.std = apply(D, 2, sd)  #column standard deviations of D     
  D.scal = scale(D, center=TRUE, scale=TRUE)  #scaled D  
  
  #biplot	
  r = 2
  main = method(D.scal,r) 	 
  Zmat = main$D.scores
  Vmat = main$D.loadings
  
  #overall quality of approximation		 
  D.hat = ((Zmat %*% t(Vmat)) * rep(D.std,each=nrow(D))) + rep(D.mean,each=nrow(D))
  dimnames(D.hat) = dimnames(D)      
  overall.quality = sum(diag(t(D.hat)%*%D.hat))  /  sum(diag(t(D)%*%D))     	
  
  #axes predictivity
  axis.pred = (diag(t(D.hat)%*%D.hat)) / (diag(t(D)%*%D))	    
  names(axis.pred) = colnames(D)         	
  
  #samples predictivity
  sample.pred = (diag(D.hat%*%t(D.hat))) / (diag(D%*%t(D)))	    
  names(sample.pred) = rownames(D) 		   
  
  #biplot points	 
  Z = Zmat #(Nx2) biplot points         
  dimnames(Z) = list(rownames(D), NULL)
  sample.style = biplot.sample.control(1, label=FALSE)
  sample.style$col = "red"
  sample.style$pch = 15
  Gmat = cbind(rep(1,nrow(D)))
  dimnames(Gmat) = list(NULL, "samples")
  classes = 1 
  
  #axes direction    
  #D.hat=ZV'
  D.axes.direction = (1 / (diag(Vmat%*%t(Vmat)))) * Vmat       	   
  
  #calibration of axes	 	
  z.axes.D = lapply(1:ncol(D), calibrate.axis, D, D.mean, D.std, D.axes.direction, 1:ncol(D), 
                    ax.tickvec.D, rep(0,ncol(D)), rep(0,ncol(D)), NULL)      
  z.axes = vector("list", ncol(D))
  for (i in 1:ncol(D)) 
  {
    z.axes[[i]] = z.axes.D[[i]]
  }
  
  #biplot axes
  ax.style = biplot.ax.control(ncol(D), c(colnames(D)))     
  ax.style$col[1:ncol(D)] = "blue"
  ax.style$label.col[1:ncol(D)] = "blue"
  ax.style$tick.col[1:ncol(D)] = "blue"
  ax.style$tick.label.col[1:ncol(D)] = "blue"
  draw.biplot(Z, G=Gmat, classes=classes, sample.style=sample.style, z.axes=z.axes, 
              ax.style=ax.style, VV=D.axes.direction)  	
  list(overall.quality=overall.quality, axis.pred=axis.pred, sample.pred=sample.pred, D.hat=D.hat)   
}


# --------------------------------------------------------------------------------------------------------------------


#A.3  Chapter 3: PLSR
# The discussed NIPALS, Kernel, and SIMPLS algorithms. 
# Example of Partial Least Squares Regression (PLSR).
# Variable Importance in the Projection (VIP) to determine which X-variables 
# are important/relevant.
# Magnitude of PLSR coefficients.
# Example of Multivariate Multiple Linear Regression (MMLR) and Principal 
# Component Regression (PCR). 
# One data set: 
#              - olive oil. 
# Users can choose their preferred PLS algorithm.  
# RGui(32-bit).

# PLSR algorithms 
#' The Nonlinear Iterative PArtial Least Squares (NIPALS) algorithm 
#' 
#' Takes in a set of predictor variables and a set of response variables and gives the Partial Least Squares (PLS) parameters.
#' @param X A (NxP) predictor matrix 
#' @param Y A (NxM) response matrix
#' @param A The number of PLS components 
#' @param ... Other arguments. Currently ignored
#' @return The PLS parameters using the NIPALS algorithm
#' @examples
#' if(require(pls)) 
#' data(oliveoil, package="pls")
#' X = as.matrix(oliveoil$chemical, ncol=5)  #predictors 
#' dimnames(X) = list(paste(c("G1","G2","G3","G4","G5","I1","I2","I3","I4","I5",
#' "S1","S2","S3","S4","S5","S6")),
#' paste(c("Acidity","Peroxide","K232","K270","DK")))   
#' Y = as.matrix(oliveoil$sensory, ncol=6)  #responses 
#' dimnames(Y) = list(paste(c("G1","G2","G3","G4","G5","I1","I2","I3","I4","I5",
#' "S1","S2","S3","S4","S5","S6")),
#' paste(c("Yellow","Green","Brown","Glossy","Transp","Syrup")))     
#' mod.NIPALS(X, Y, A=5)  
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export
mod.NIPALS = function(X, Y, A, ...)
{
  #A=initial number of components    
  X.cen = scale(X, center=TRUE, scale=FALSE)  #centered X
  Y.cen = scale(Y, center=TRUE, scale=FALSE)  #centered Y
  Y.mean = colMeans(Y)  #(1xM) column means of Y
  Y.hat = array(0, dim=c(nrow(Y), ncol(Y), A))  #predicted Y-values     
  U = matrix(0, nrow=nrow(Y), ncol=A)  #(NxA) Y.scores matrix	
  W = P = matrix(0, nrow=ncol(X), ncol=A)  #(PxA) X.weights and X.loadings
  #matrices	
  C = Q = matrix(0, nrow=ncol(Y), ncol=A)  #(MxA) Y.weights and Y.loadings
  # matrices
  rmsep = array(NA, dim=c(1,A)) 	#(1XA) Root Mean Squared Error of Prediction (RMSEP) value								 
  
  #steps 2 to 9
  X.a = X.cen
  Y.a = Y.cen
  for (a in 1:A) 
  {	
    u.a = Y.a[, which.max(colSums(Y.a * Y.a))]	   
    t.a.old = 0
    repeat 
    {
      w.a = t(X.a) %*% u.a  / as.numeric(sqrt(t(t(X.a)%*%u.a)%*%t(X.a) %*% u.a)) 
      t.a.new = X.a %*% w.a  / as.numeric(sqrt(t(X.a %*% w.a)%*%X.a %*% w.a))  
      c.a = t(Y.a) %*% t.a.new  / as.numeric(sqrt(t(t(Y.a)%*%t.a.new)%*%t(Y.a)%*%t.a.new)) 
      if (sum(abs(t.a.new - t.a.old)/abs(t.a.new)) < 10^-8)
        break
      else 
      { 
        u.a = Y.a %*% c.a  / as.numeric(sqrt(t(Y.a %*% c.a)%*%Y.a %*% c.a)) 
        t.a.old = t.a.new						 	         
      } 					  
    }
    p.a = t(X.a) %*% t.a.new  #(Px1) X.loading
    q.a = t(Y.a) %*% t.a.new  #(Mx1) Y.loading
    #update X.a and Y.a 
    X.a = X.a - t.a.new%*%t(p.a)
    Y.a = Y.a - t.a.new%*%t(c.a)
    
    W[, a] = w.a   #(PxA) X.weights matrix       
    C[, a] = c.a   #(MxA) Y.weights matrix
    P[, a] = p.a   #(PxA) X.loadings matrix
    Q[, a] = q.a   #(MxA) Y.loadings matrix
    U[, a] = u.a   #(NxA) Y.scores matrix	
  }
  R = W %*% solve(t(P)%*%W)  #(PxA) transformed X.weights matrix 
  T = X.cen %*% R   #(NxA) X.scores matrix
  for (a in 1:A)
  {
    Y.hat[, , a] = (T[, 1:a] %*% solve(t(T[, 1:a]) %*% T[, 1:a]) %*% t(T[, 1:a]) %*% Y.cen) + 
      rep(Y.mean, each=nrow(Y))  #predicted Y-values using (1:A) components
    rmsep[,a] = sqrt(sum((Y - Y.hat[,,a])^2) / nrow(Y))  #Root Mean Squared Error of Prediction (RMSEP)
    dimnames(rmsep) = list(paste("RMSEP.value"), paste(1:A, "Comps"))
  }	
  
  dimnames(Y.hat) = list(rownames(Y), colnames(Y), paste(1:A, "Comps")) 
  dimnames(T) = list(rownames(X), paste("Comp", 1:A, seq=""))
  dimnames(U) = list(rownames(Y), paste("Comp", 1:A, seq=""))
  dimnames(W) = dimnames(P) = dimnames(R) = list(colnames(X), paste("Comp", 1:A, seq=""))
  dimnames(Q) = list (colnames(Y), paste("Comp", 1:A, seq="")) 	   
  
  list(X.scores=T, Y.scores=U, X.weights=W, RMSEP=rmsep, X.weights.trans=R, 
       X.loadings=P, Y.loadings=Q, Y.hat=Y.hat) 
} 

#' The Kernel algorithm for few(er) samples but large variables by Rannar et al. (1994)
#' 
#' Takes in a set of predictor variables and a set of response variables and gives the Partial Least Squares (PLS) parameters.
#' @param X A (NxP) predictor matrix 
#' @param Y A (NxM) response matrix
#' @param A The number of PLS components 
#' @param ... Other arguments. Currently ignored
#' @return The PLS parameters using the Kernel algorithm by RC$nnar et al. (1994)
#' @examples 
#' if(require(pls)) 
#' data(oliveoil, package="pls")
#' X = as.matrix(oliveoil$chemical, ncol=5) 
#' dimnames(X) = list(paste(c("G1","G2","G3","G4","G5","I1","I2","I3","I4","I5",
#' "S1","S2","S3","S4","S5","S6")),
#' paste(c("Acidity","Peroxide","K232","K270","DK")))   
#' Y = as.matrix(oliveoil$sensory, ncol=6)  
#' dimnames(Y) = list(paste(c("G1","G2","G3","G4","G5","I1","I2","I3","I4","I5",
#' "S1","S2","S3","S4","S5","S6")),
#' paste(c("Yellow","Green","Brown","Glossy","Transp","Syrup")))     
#' mod.KernelPLS_R(X, Y, A=5)  
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export
mod.KernelPLS_R = function(X, Y, A, ...)
{
  #A=initial number of components    
  X.cen = scale(X, center=TRUE, scale=FALSE)  #centered X
  Y.cen = scale(Y, center=TRUE, scale=FALSE)  #centered Y
  Y.mean = colMeans(Y)  #(1xM) column means of Y   
  Y.hat = array(0, dim=c(nrow(Y), ncol(Y), A)) #predicted Y-values	   
  U = matrix(0, nrow=nrow(Y), ncol=A)	#(NxA) Y.scores matrix	
  W = P = matrix(0, nrow=ncol(X), ncol=A) #(PxA) X.weights and X.loadings 
  #matrices	
  Q = matrix(0, nrow=ncol(Y), ncol=A)   #(MxA) Y.loadings matrix	
  rmsep = array(NA, dim=c(1,A)) 	#(1XA) Root Mean Squared Error of Prediction (RMSEP) value	
  X.c = X.cen 
  Y.c = Y.cen
  XXt = X.c %*% t(X.c)
  YYt = Y.c %*% t(Y.c)	  
  #steps 2 to 6
  for (a in 1:A) 
  {	    
    t.a = cbind(as.numeric(eigen(XXt%*%YYt)$vector[,1])) #(Nx1) X.score
    u.a = YYt %*% t.a / as.numeric(sqrt(t(YYt%*%t.a)%*%YYt%*%t.a))	#(Nx1) Y.score       
    #update XXt and YYt 
    G.a = diag(1,nrow=nrow(X)) - t.a%*%t(t.a)	
    XXt = G.a %*% XXt %*% G.a
    YYt = G.a %*% YYt %*% G.a   
    
    w.a = t(X.c) %*% u.a  #(Px1) X.weight
    p.a = t(X.c) %*% t.a  #(Px1) X.loading
    q.a = t(Y.c) %*% t.a  #(Mx1) Y.loading
    U[, a] = u.a  #(NxA) Y.score matrix
    W[, a] = w.a  #(PxA) X.weight matrix   
    P[, a] = p.a  #(PxA) X.loading matrix    
    Q[, a] = q.a  #(MxA) Y.loading matrix   	     		   
  }  
  R = W %*% solve(t(P)%*%W)  #(PxA) transformed X.weights matrix 
  T = X.cen %*% R   #(NxA) X.scores matrix
  for (a in 1:A)
  {
    Y.hat[, , a] = (T[, 1:a] %*% solve(t(T[, 1:a]) %*% T[, 1:a]) %*% t(T[, 1:a]) %*% Y.cen) + 
      rep(Y.mean, each=nrow(Y))  #predicted Y-values using (1:A) components
    rmsep[,a] = sqrt(sum((Y - Y.hat[,,a])^2) / nrow(Y))  #Root Mean Squared Error of Prediction (RMSEP)
    dimnames(rmsep) = list(paste("RMSEP.value"), paste(1:A, "Comps"))					   
  }	
  dimnames(Y.hat) = list(rownames(Y), colnames(Y), paste(1:A, "Comps")) 
  dimnames(T) = list(rownames(X), paste("Comp", 1:A, seq=""))
  dimnames(U) = list(rownames(Y), paste("Comp", 1:A, seq=""))
  dimnames(W) = dimnames(P) = dimnames(R) = list(colnames(X), paste("Comp", 1:A, seq=""))
  dimnames(Q) = list(colnames(Y), paste("Comp", 1:A, seq="")) 	     
  
  list(X.scores=T, Y.scores=U, X.weights=W, RMSEP=rmsep, X.weights.trans=R, 
       X.loadings=P, Y.loadings=Q, Y.hat=Y.hat) 
}   

#' The Kernel algorithm for few(er) variables but large samples by Lindgren et al. (1993) 
#' 
#' Takes in a set of predictor variables and a set of response variables and gives the Partial Least Squares (PLS) parameters.
#' @param X A (NxP) predictor matrix 
#' @param Y A (NxM) response matrix
#' @param A The number of PLS components 
#' @param ... Other arguments. Currently ignored
#' @return The PLS parameters using the kernel algorithm by Lindgren et al. (1993)
#' @examples 
#' if(require(pls))
#' data(oliveoil, package="pls")
#' X = as.matrix(oliveoil$chemical, ncol=5) 
#' dimnames(X) = list(paste(c("G1","G2","G3","G4","G5","I1","I2","I3","I4","I5",
#' "S1","S2","S3","S4","S5","S6")),
#' paste(c("Acidity","Peroxide","K232","K270","DK")))   
#' Y = as.matrix(oliveoil$sensory, ncol=6)  
#' dimnames(Y) = list(paste(c("G1","G2","G3","G4","G5","I1","I2","I3","I4","I5",
#' "S1","S2","S3","S4","S5","S6")),
#' paste(c("Yellow","Green","Brown","Glossy","Transp","Syrup")))     
#' mod.KernelPLS_L(X, Y, A=5)  
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export
mod.KernelPLS_L = function(X, Y, A, ...)
{
  #A=initial number of components    
  X.cen = scale(X, center=TRUE, scale=FALSE)  #centered X
  Y.cen = scale(Y, center=TRUE, scale=FALSE)  #centered Y
  Y.mean = colMeans(Y)  #(1xM) column means of Y   
  Y.hat = array(0, dim=c(nrow(Y), ncol(Y), A)) #predicted Y-values    
  W = P = matrix(0, nrow=ncol(X), ncol=A) #(PxA) X.weights and X.loadings 
  #matrices	
  Q = matrix(0, nrow=ncol(Y), ncol=A)   #(MxA) Y.loadings matrix	
  rmsep = array(NA, dim=c(1,A)) 	#(1XA) Root Mean Squared Error of Prediction (RMSEP) value	
  X.c = X.cen 
  Y.c = Y.cen
  XtX = t(X.c) %*% (X.c)
  XtY = t(X.c) %*% (Y.c)	  
  #steps 2 to 6
  for (a in 1:A) 
  {	    
    w.a = cbind(as.numeric(eigen(XtY%*%t(XtY))$vector[,1])) #(Nx1) X.score 
    p.a = XtX %*% w.a  / as.numeric(sqrt(t(XtX%*%w.a)%*%XtX%*%w.a)) #(Px1) X.loading
    q.a = t(XtY) %*% w.a / as.numeric(sqrt(t(t(XtY)%*%w.a)%*%t(XtY)%*%w.a))  #(Mx1) Y.loading	   
    #update XtX and XtY 
    O.a = diag(1,nrow=ncol(X)) - w.a%*%t(p.a)	
    XtX = t(O.a) %*% XtX %*% O.a
    XtY = t(O.a) %*% XtY 
    
    W[, a] = w.a  #(PxA) X.weights matrix   
    P[, a] = p.a  #(PxA) X.loadings matrix    
    Q[, a] = q.a  #(MxA) Y.loadings matrix   	     		   
  }  
  R = W %*% solve(t(P)%*%W)  #(PxA) transformed X.weights matrix 
  T = X.cen %*% R   #(NxA) X.scores matrix
  for (a in 1:A)
  {
    Y.hat[, , a] = (T[, 1:a] %*% solve(t(T[, 1:a]) %*% T[, 1:a]) %*% t(T[, 1:a]) %*% Y.cen) + 
      rep(Y.mean, each=nrow(Y))  #predicted Y-values using (1:A) components
    rmsep[,a] = sqrt(sum((Y - Y.hat[,,a])^2) / nrow(Y))  #Root Mean Squared Error of Prediction (RMSEP)
    dimnames(rmsep) = list(paste("RMSEP.value"), paste(1:A, "Comps"))					  
  }	     
  dimnames(Y.hat) = list(rownames(Y), colnames(Y), paste(1:A, "Comps")) 
  dimnames(T) = list(rownames(X), paste("Comp", 1:A, seq=""))
  dimnames(W) = dimnames(P) = dimnames(R) = list(colnames(X), paste("Comp", 1:A, seq=""))
  dimnames(Q) = list(colnames(Y), paste("Comp", 1:A, seq="")) 	     
  
  list(X.scores=T, X.weights=W, RMSEP=rmsep, X.weights.trans=R, X.loadings=P,
       Y.loadings=Q, Y.hat=Y.hat) 
} 

#' The Statistical Inspired Modification to Partial Least Squares (SIMPLS) algorithm 
#' 
#' Takes in a set of predictor variables and a set of response variables and gives the Partial Least Squares (PLS) parameters.
#' @param X A (NxP) predictor matrix 
#' @param Y A (NxM) response matrix
#' @param A The number of PLS components 
#' @param ... Other arguments. Currently ignored
#' @return The PLS parameters using the SIMPLS algorithm
#' @examples  
#' if(require(pls)) 
#' data(oliveoil, package="pls")
#' X = as.matrix(oliveoil$chemical, ncol=5)   
#' dimnames(X) = list(paste(c("G1","G2","G3","G4","G5","I1","I2","I3","I4","I5",
#' "S1","S2","S3","S4","S5","S6")),
#' paste(c("Acidity","Peroxide","K232","K270","DK")))   
#' Y = as.matrix(oliveoil$sensory, ncol=6) 
#' dimnames(Y) = list(paste(c("G1","G2","G3","G4","G5","I1","I2","I3","I4","I5",
#' "S1","S2","S3","S4","S5","S6")),
#' paste(c("Yellow","Green","Brown","Glossy","Transp","Syrup")))   
#' #final number of PLS components  
#' RMSEP = mod.SIMPLS(X, Y, A=5)$RMSEP #RMSEP values
#' plot(t(RMSEP), type = "b", xlab="Number of components", ylab="RMSEP  values")     
#' A.final = 2 #from the RMSEP plot 
#' #PLS matrices R, P, T, Q, and Y.hat from SIMPLS algorithm
#' options(digits=3)
#' mod.SIMPLS(X, Y, A=A.final)
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export
mod.SIMPLS = function(X, Y, A, ...)
{
  #A=initial number of components    
  Xc = scale(X, center=TRUE, scale=FALSE)  #centered X
  Yc = scale(Y, center=TRUE, scale=FALSE)  #centered Y
  Y.mean = colMeans(Y)  #(1xM) column means of Y   
  Y.hat = array(0, dim=c(nrow(Y), ncol(Y), A)) #predicted Y-values	
  T = matrix(0, nrow=nrow(X), ncol=A)   #(NxA) X.scores matrix
  R = P = matrix(0, nrow=ncol(X), ncol=A)  #(PxA) transformed X.weights and 
  #X.loadings matrices 
  Q = matrix(0, nrow=ncol(Y), ncol=A)  #(MxA) Y.loadings matrix 	
  rmsep = array(NA, dim=c(1,A)) 	#(1XA) Root Mean Squared Error of Prediction (RMSEP) value	
  S = t(Xc) %*% Yc	  
  #steps 2 to 4
  for (a in 1:A) 
  {	    
    c.a = svd(S)$v[,1]  #(Mx1) Y.weight
    r.a = svd(S)$u[,1]  #(Px1) X.weight         
    t.a = Xc %*% r.a   #(Nx1) X.score
    t.a.norm = as.numeric(sqrt(t(t.a)%*%t.a))
    t.a = t.a  / t.a.norm   #normalized
    p.a = t(Xc) %*% t.a  #(Px1) X.loading
    q.a = t(Yc) %*% t.a  #(Mx1) Y.loading 
    
    #update S
    S = S - p.a%*%solve(t(p.a)%*%p.a)%*%t(p.a)%*%S	   
    
    #step 5: obtaining the PLS parameters (scores, loadings, weights and 
    #PLSR coefficients) in matrix form
    R[, a] = r.a   #(PxA) X.weights matrix   
    T[, a] = t.a   #(NxA) X.scores matrix	         	   
    P[, a] = p.a   #(PxA) X.loadings matrix
    Q[, a] = q.a   #(MxA) Y.loadings matrix
    Y.hat[, , a] = (T[, 1:a] %*% solve(t(T[, 1:a]) %*% T[, 1:a]) %*% t(T[, 1:a]) %*% Yc) + 
      rep(Y.mean, each=nrow(Y))   #predicted Y-values using (1:A) components
    rmsep[,a] = sqrt(sum((Y - Y.hat[,,a])^2) / nrow(Y))  #Root Mean Squared Error of Prediction (RMSEP)
    dimnames(rmsep) = list(paste("RMSEP.value"), paste(1:A, "Comps"))
  }
  
  dimnames(Y.hat) = list(rownames(Y), colnames(Y), paste(1:A, "Comps")) 
  dimnames(R) = dimnames(P) = list(colnames(X), paste("Comp", 1:A, seq=""))	 
  dimnames(T) = list(rownames(X), paste("Comp", 1:A, seq=""))
  dimnames(Q) = list(colnames(Y), paste("Comp", 1:A, seq="")) 
  list(X.scores=T, X.weights.trans=R, X.loadings=P, Y.loadings=Q, 
       Y.hat=Y.hat, RMSEP=rmsep) 
} 

#' The Variable Importance in the Projection (VIP) values
#' 
#' Takes in a set of predictor variables and a set of response variables and gives the VIP values for the predictor variables.
#' @param X A (NxP) predictor matrix 
#' @param Y A (NxM) response matrix 
#' @param A The number of Partial Least Squares (PLS) components 
#' @param algorithm Any of the PLS algorithms ("mod.SIMPLS","mod.NIPALS", "mod.KernelPLS_R", "mod.KernelPLS_L") 
#' @param cutoff desired cut off value to use for selecting the important X-variables
#' @param ... Other arguments. Currently ignored
#' @return The VIP value for each of the X-variables
#' @examples 
#' if(require(chemometrics))
#' data(cereal, package="chemometrics")
#' X = as.matrix(cbind(cereal$X))  
#' Y = as.matrix(cbind(cereal$Y))  
#' main2 = mod.VIP(X=X, Y=Y, algorithm=mod.SIMPLS, A=2, cutoff=0.8) 
#' main2
#' X.new = X[,c(main2$X.impor)]  #important X-variables  
#' X.new
#'
#' #nutrimouse data
#' if(require(mixOmics))
#' data(nutrimouse, package="mixOmics")  
#' X1 = as.matrix(nutrimouse$lipid, ncol=21) 
#' Y1 = as.matrix(nutrimouse$gene, ncol=120) 
#' main = mod.SIMPLS(X=X1, Y=Y1, A=17) #using the SIMPLS algorithm
#' #RMSEP
#' RMSEP = main$RMSEP 
#' plot(t(RMSEP), type = "b", xlab="Number of components", ylab="RMSEP  values")   
#' A.final = 9 #from the RMSEP plot 
#' #Final PLSR
#' mod.SIMPLS(X=X1, Y=Y1, A=A.final)    
#' #VIP 
#' main2 = mod.VIP(X=X1, Y=Y1, algorithm=mod.SIMPLS, A=A.final, cutoff=0.8) 
#' main2
#' X.new = X1[,c(main2$X.impor)]  #important X-variables  
#' X.new
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export 
mod.VIP = function(X, Y, algorithm=NULL, A, cutoff=NULL, ... ) 
{ 
  options(digits=3)        	   	   
  X.scal = scale(X, center=TRUE, scale=TRUE)  #scaled X      	
  Y.scal = scale(Y, center=TRUE, scale=TRUE)  #scaled Y	
  
  #PLS and PLSR parameters
  main = algorithm(X.scal,Y.scal,A)	
  Rmat = main$X.weights.trans  #(PxA) transformed X.weights matrix
  Tmat = main$X.scores	#(NxA) X.scores matrix   
  Qmat = main$Y.loadings #(MxA) Y.loadings matrix
  Bmat = Rmat %*% t(Qmat)  #(PxM) estimated PLSR coefficients matrix	     
  dimnames(Bmat) = list(colnames(X), colnames(Y))     	  
  
  #VIP values  
  P = ncol(X) #or nrow(Bmat)
  VIP = array(0, dim=P)  #(Px1) VIP values vector 
  for (i in 1:P)
  {
    vip = sqrt( (P / (t(Bmat[i,]) %*% Bmat[i,] %*% sum(t(Tmat[,1:A]) %*% Tmat[,1:A]))) * 
                  (t(Bmat[i,]) %*% Bmat[i,] %*% sum(t(Tmat[,1:A]) %*% Tmat[,1:A] %*% ((
                    Rmat[i,1:A])^2))) )
    VIP[i] = vip  #(Px1) VIP value vector
  }
  names(VIP)=rownames(Bmat)     
  list(VIP.values=VIP, X.impor=which(VIP >= cutoff))
} 

#' Magnitude of the Partial Least Squares Regression (PLSR) coefficients matrix
#' 
#' Takes in a set of predictor variables and a set of response variables and produces the mean plot of the absolute values of the PLSR coefficients matrix.
#' @param X A (NxP) predictor matrix 
#' @param Y A (NxM) response matrix
#' @param A The number of Partial Least Squares (PLS) components
#' @param algorithm Any of the PLS algorithms ("mod.NIPALS", "mod.KernelPLS_R", "mod.KernelPLS_L", "mod.SIMPLS") 
#' @param ... Other arguments. Currently ignored
#' @return The mean plot of the absolute values of the PLSR coefficients matrix
#' @examples 
#' if(require(pls))
#' data(oliveoil, package="pls")
#' X = as.matrix(oliveoil$chemical, ncol=5) 
#' dimnames(X) = list(paste(c("G1","G2","G3","G4","G5","I1","I2","I3","I4","I5",
#' "S1","S2","S3","S4","S5","S6")),
#' paste(c("Acidity","Peroxide","K232","K270","DK")))   
#' Y = as.matrix(oliveoil$sensory, ncol=6)  
#' dimnames(Y) = list(paste(c("G1","G2","G3","G4","G5","I1","I2","I3","I4","I5",
#' "S1","S2","S3","S4","S5","S6")),
#' paste(c("Yellow","Green","Brown","Glossy","Transp","Syrup")))     
#' Mag.Bmat.plot(X, Y, algorithm=mod.SIMPLS, A=2) 
#'
#' #nutrimouse data
#' if(require(mixOmics))
#' data(nutrimouse, package="mixOmics")  
#' X1 = as.matrix(nutrimouse$lipid, ncol=21) 
#' Y1 = as.matrix(nutrimouse$gene, ncol=120)  
#' #VIP 
#' A.final = 9
#' main2 = mod.VIP(X=X1, Y=Y1, algorithm=mod.SIMPLS, A=A.final, cutoff=0.8) 
#' X.new = X1[,c(main2$X.impor)]  #important X-variables  
#' Mag.Bmat.plot(X=X.new, Y1, algorithm=mod.SIMPLS, A=A.final)
#' #alternatively
#' X.scal = scale(X.new, center=TRUE, scale=TRUE)
#' Y.scal = scale(Y1, center=TRUE, scale=TRUE)
#' main3 = mod.SIMPLS(X.scal, Y.scal, A.final)
#' Bmat = main3$X.weights.trans %*% t(main3$Y.loadings)  #PLSR coefficients matrix
#' dimnames(Bmat) = list(colnames(X.new), colnames(Y1))
#' Abs.Bmat = abs(Bmat) #absolute values of the coefficients
#' rowMeans(Abs.Bmat)
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export 
Mag.Bmat.plot = function(X, Y, algorithm=NULL, A, ... ) 
{ 
  options(digits=3)        	   	   
  X.scal = scale(X, center=TRUE, scale=TRUE)  #scaled X      	
  Y.scal = scale(Y, center=TRUE, scale=TRUE)  #scaled Y	
  
  #PLS and PLSR parameters
  main = algorithm(X.scal,Y.scal,A)	
  Rmat = main$X.weights.trans  #(PxA) transformed X.weights matrix
  Qmat = main$Y.loadings #(MxA) Y.loadings matrix
  Bmat = Rmat %*% t(Qmat)  #(PxM) estimated PLSR coefficients matrix	     
  dimnames(Bmat) = list(colnames(X), colnames(Y)) 	
  
  #Absolute values of the PLSR coefficients matrix	 
  Abs.Bmat = abs(Bmat) #absolute values	 
  #graphically
  plot(rowMeans(Abs.Bmat), type="b", ylim=c(range(Abs.Bmat)*1.5), xaxt="n", xlab="", 
       ylab="Means values", main="", asp=1)      
  lines(rowMeans(Abs.Bmat), col="green", type="b")
  text(rowMeans(Abs.Bmat), lab=rownames(Abs.Bmat), pos=3, cex = 0.9, col="green") 
  legend("topright", legend="absolute B.PLS", col="green", pch=1, lty=1)	
}

#' Multivariate Multiple Linear Regression (MMLR) 
#' 
#' Takes in a set of predictor variables and a set of response variables and gives the MMLR parameters.
#' @param X A (NxP) predictor matrix 
#' @param Y A (NxM) response matrix
#' @param ... Other arguments. Currently ignored
#' @return The MMLR parameters
#' @examples
#' if(require(pls))
#' data(oliveoil, package="pls")
#' X = as.matrix(oliveoil$chemical, ncol=5) 
#' dimnames(X) = list(paste(c("G1","G2","G3","G4","G5","I1","I2","I3","I4","I5",
#' "S1","S2","S3","S4","S5","S6")),
#' paste(c("Acidity","Peroxide","K232","K270","DK")))   
#' Y = as.matrix(oliveoil$sensory, ncol=6)  
#' dimnames(Y) = list(paste(c("G1","G2","G3","G4","G5","I1","I2","I3","I4","I5",
#' "S1","S2","S3","S4","S5","S6")),
#' paste(c("Yellow","Green","Brown","Glossy","Transp","Syrup")))     
#' mod.MMLR(X, Y)
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export
mod.MMLR = function(X, Y, ...)
{      
  X.cen = scale(X, center=TRUE, scale=FALSE)  #centered X
  Y.cen = scale(Y, center=TRUE, scale=FALSE)  #centered Y
  Y.mean = colMeans(Y)  #(1xM) column means of Y   
  
  Y.hat = (X.cen %*% solve(t(X.cen)%*%X.cen) %*% t(X.cen) %*% Y.cen) + 
    rep(Y.mean, each=nrow(Y))   #predicted Y-values 
  B.mmlr = solve(t(X.cen)%*%X.cen) %*% t(X.cen) %*% Y.cen  #(PxM) MMLR coefficients matrix			 
  
  dimnames(Y.hat) = list(rownames(Y), colnames(Y))  
  dimnames(B.mmlr) = list(colnames(X), colnames(Y))  
  list(Y.hat=Y.hat, B.mmlr=B.mmlr) 
}	

#' Principal Component Regression (PCR)
#' 
#' Takes in a set of predictor variables and a set of response variables and gives the PCR parameters.
#' @param X A (NxP) predictor matrix 
#' @param Y A (NxM) response matrix
#' @param r The number of PCA components
#' @param ... Other arguments. Currently ignored
#' @return The PCR parameters
#' @examples 
#' if(require(pls))
#' data(oliveoil, package="pls")
#' X = as.matrix(oliveoil$chemical, ncol=5) 
#' dimnames(X) = list(paste(c("G1","G2","G3","G4","G5","I1","I2","I3","I4","I5",
#' "S1","S2","S3","S4","S5","S6")),
#' paste(c("Acidity","Peroxide","K232","K270","DK")))   
#' Y = as.matrix(oliveoil$sensory, ncol=6)  
#' dimnames(Y) = list(paste(c("G1","G2","G3","G4","G5","I1","I2","I3","I4","I5",
#' "S1","S2","S3","S4","S5","S6")),
#' paste(c("Yellow","Green","Brown","Glossy","Transp","Syrup")))     
#' mod.PCR(X, Y, r=2)
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export
mod.PCR = function(X, Y, r, ...)
{
  #r = initial number of components    
  X.cen = scale(X, center=TRUE, scale=FALSE)  #centered X
  Y.cen = scale(Y, center=TRUE, scale=FALSE)  #centered Y
  Y.mean = colMeans(Y)  #(1xM) column means of Y     
  Y.hat = array(0, dim=c(nrow(Y), ncol(Y), r)) #predicted Y-values
  B.pcr = array(0, dim=c(ncol(X), ncol(Y), r)) #(PxM) predicted PCR 
  # coefficients matrix 	  
  X.svd = svd(X.cen)
  V = X.svd$v   #(Pxr) X.loadings matrix   
  Z = X.cen %*% V  #(Nxr) X.scores matrix
  
  for (i in 1:r) 
  {	    
    Y.hat[, , i] = (Z[, 1:i] %*% solve(t(Z[, 1:i]) %*% Z[, 1:i]) %*% t(Z[, 1:i]) %*% Y.cen) + 
      rep(Y.mean, each=nrow(Y))   #predicted Y-values using (1:r) components
    B.pcr[, , i] = V[, 1:i] %*% solve(t(Z[, 1:i]) %*% Z[, 1:i]) %*% t(Z[, 1:i]) %*% Y.cen  
    #(PxM) PCR coefficients matrix	   
  }	 
  
  dimnames(Y.hat) = list(rownames(Y), colnames(Y), paste(1:r, "Comps")) 
  dimnames(Z) = list(rownames(X), paste("Comp", 1:ncol(Z), seq=""))	   
  dimnames(V) = list(colnames(X), paste("Comp", 1:ncol(V), seq=""))	  
  dimnames(B.pcr) = list(colnames(X), colnames(Y), paste(1:r, "Comps")) 
  list(X.scores=Z, X.loadings=V, Y.hat=Y.hat, B.pcr=B.pcr) 
}	 	


# --------------------------------------------------------------------------------------------------------------------


#A.4  Chapter  4: Covariance biplots
# Example to illustrate the covariance biplot. 
# Example to illustrate the covariance monoplot.
# One data sets:
#               - olive oil.
# R (32-bit).

#' The covariance biplot
#' 
#' Takes in a set of predictor variables and a set of response variables and produces a covariance biplot.
#' @param X A (NxP) predictor matrix 
#' @param Y A (NxM) response matrix
#' @param ... Other arguments. Currently ignored
#' @return The covariance biplot of X and Y
#' @examples 
#' if(require(pls))
#' data(oliveoil, package="pls")
#' X = as.matrix(oliveoil$chemical, ncol=5) 
#' dimnames(X) = list(paste(c("G1","G2","G3","G4","G5","I1","I2","I3","I4",
#' "I5","S1","S2","S3","S4","S5","S6")),
#' paste(c("Acidity","Peroxide","K232","K270","DK")))   
#' Y = as.matrix(oliveoil$sensory, ncol=6)  
#' dimnames(Y) = list(paste(c("G1","G2","G3","G4","G5","I1","I2","I3","I4",
#' "I5","S1","S2","S3","S4","S5","S6")),
#' paste(c("Yellow","Green","Brown","Glossy","Transp","Syrup")))     
#' cov.biplot(X, Y)
#'
#' #cocktail data
#' if(require(SensoMineR))
#' data(cocktail, package="SensoMineR")
#' X3 = as.matrix(compo.cocktail, ncol=4)  
#' Y3 = as.matrix(senso.cocktail, ncol=13) 
#' cov.biplot(X3,Y3)
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export
#X'Y = UDV'
#G=UD^0.5 and H=VD^0.5 ==> sub-optimal representation for both X and Y
#covariance of two sets of variables at a time = covariance BETWEEN two 
#sets of variables
cov.biplot = function(X, Y, ...)
{
  options(digits=3)  
  X.scal = scale(X, center=TRUE, scale=TRUE)  #scaled X  
  Y.scal = scale(Y, center=TRUE, scale=TRUE)  #scaled Y  
  S = cov(X.scal,Y.scal)  #(PxM) covariance matrix
  S.svd = svd(S)
  G = S.svd$u %*% (diag(S.svd$d, nrow=ncol(t(S.svd$d))))^0.5
  dimnames(G) = list(colnames(X), paste("Comp", 1:ncol(G), seq=""))      
  H = S.svd$v %*% (diag(S.svd$d, nrow=ncol(t(S.svd$d))))^0.5
  dimnames(H) = list(colnames(Y), paste("Comp", 1:ncol(H), seq=""))         
  plot(rbind(G,H), xlim=range(rbind(G,H))*1.5, ylim=range(rbind(G,H))*1.5, xlab="", ylab="", 
       xaxs="i", yaxs="i", xaxt="n", yaxt="n", type="n", pch=15, asp=1, main="")
  #for X-variables
  arrows(G[, 1]*0, G[, 2]*0, G[, 1], G[, 2], length=0.09, col="blue")
  text(G[, 1], G[, 2], labels=rownames(G), col="blue", pos="2", cex=0.6)
  arrows(G[, 1]*0, G[, 2]*0, -G[, 1], -G[, 2], length=0, col="blue")
  #for Y-variables
  arrows(H[, 1]*0, H[, 2]*0, H[, 1], H[, 2], length=0.09, col="green")
  text(H[, 1], H[, 2], labels=rownames(H), col="green", pos="2", cex=0.6) 
  arrows(H[, 1]*0, H[, 2]*0, -H[, 1], -H[, 2], length=0, col="green")        
  list(G__UDhalf=G[,1:2], H__VDhalf=H[,1:2]) 
}

#' The covariance monoplot
#' 
#' Takes in only one set of variables (e.g., predictors) and produces a covariance monoplot.
#' @param X A (NxP) predictor matrix 
#' @param ... Other arguments. Currently ignored
#' @return The covariance monoplot of X
#' @examples 
#' if(require(pls))
#' data(oliveoil, package="pls")
#' Y = as.matrix(oliveoil$sensory, ncol=6)  
#' dimnames(Y) = list(paste(c("G1","G2","G3","G4","G5","I1","I2","I3","I4",
#' "I5","S1","S2","S3","S4","S5","S6")),
#' paste(c("Yellow","Green","Brown","Glossy","Transp","Syrup")))     
#' cov.monoplot(Y)
#'
#' #cocktail data
#' if(require(SensoMineR))
#' data(cocktail, package="SensoMineR")
#' Y3 = as.matrix(senso.cocktail, ncol=13) 
#' cov.monoplot(Y3)
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export
#for X=UDV', X'X = VDDV'
#G=VD and H=VD ==> optimal representation for X
#covariance of one set of variables at a time = covariance WITHIN one set 
#of variables = Variances
cov.monoplot = function(X, ...)
{
  options(digits=3)  
  X.scal = scale(X, center=TRUE, scale=TRUE)  #scaled X 
  S = cov(X.scal,X.scal)  #(PxP) Variance matrix
  S.svd = svd(S)
  G = S.svd$u %*% (diag(S.svd$d, nrow=nrow(S)))
  dimnames(G) = list(colnames(X), paste("Comp", 1:ncol(G), seq=""))      
  H = S.svd$v %*% (diag(S.svd$d, nrow=nrow(S)))
  dimnames(H) = list(colnames(X), paste("Comp", 1:ncol(H), seq=""))  
  plot(G, xlim=range(G)*1.5, ylim=range(G)*1.5, xlab="",ylab="", xaxs="i", yaxs="i", 
       xaxt="n", yaxt="n", type="n", pch=15, asp=1, main="")
  #for X-variables
  arrows(G[, 1]*0, G[, 2]*0, G[, 1], G[, 2], length=0.09, col="green")
  text(G[, 1], G[, 2], labels=rownames(G), col="green", pos="2", cex=0.6)
  arrows(G[, 1]*0, G[, 2]*0, -G[, 1], -G[, 2], length=0, col="green")          
  list(G__VD=G[,1:2], H__VD=H[,1:2]) 
}    


# --------------------------------------------------------------------------------------------------------------------


#A.5  Chapter  5: PLS biplots
# Example to illustrate the SIMPLS and Kernel PLS biplots. 
# One data set: 
#              - olive oil. 
# Area biplot triangle for reading off the approximated values for the 
# coefficient bi in the PLS biplot, where, i=(1,2,...P).
# PLS biplot with No sample names
# PLS biplot with No sample and tick markers names
# Users can choose their:
#                     (i) preferred PLS algorithm, 
#                    (ii) length of tick markers (i.e. "ax.tickvec") for 
#                         variables in the PLS biplots, 
#                   (iii) preferred PLSR coefficient values bi to approximate
#                         using area biplot triangles, and   
#                    (iv) preferred Y-variable(s) to use as base(s). 
# R (32-bit).

# PLS biplots
#' The Partial Least Squares (PLS) biplot 
#' 
#' Takes in a set of predictor variables and a set of response variables and produces a PLS biplot.
#' @param X A (NxP) predictor matrix 
#' @param Y A (NxM) response matrix
#' @param algorithm Any of the PLS algorithms ("mod.NIPALS", "mod.KernelPLS_R", "mod.KernelPLS_L", "mod.SIMPLS") 
#' @param ax.tickvec.X tick marker length for each X-variable axis in the PLS biplot
#' @param ax.tickvec.Y tick marker length for each Y-variable axis in the PLS biplot
#' @param ... Other arguments. Currently ignored
#' @return The PLS biplot of D=[X Y] with some parameters
#' @examples 
#' if(require(pls))
#' data(oliveoil, package="pls")
#' X = as.matrix(oliveoil$chemical, ncol=5) 
#' dimnames(X) = list(paste(c("G1","G2","G3","G4","G5","I1","I2","I3","I4","I5",
#' "S1","S2","S3","S4","S5","S6")),
#' paste(c("Acidity","Peroxide","K232","K270","DK")))   
#' Y = as.matrix(oliveoil$sensory, ncol=6)  
#' dimnames(Y) = list(paste(c("G1","G2","G3","G4","G5","I1","I2","I3","I4","I5",
#' "S1","S2","S3","S4","S5","S6")),
#' paste(c("Yellow","Green","Brown","Glossy","Transp","Syrup")))     
#' #SIMPLS biplot
#' PLS.biplot(X, Y, algorithm=mod.SIMPLS, ax.tickvec.X=c(8,5,5,5,5), ax.tickvec.Y=c(5,8,5,6,9,8))  
#' #Kernel PLS biplot
#' PLS.biplot(X, Y, algorithm=mod.KernelPLS_R, ax.tickvec.X=c(3,3,4,5,2), ax.tickvec.Y=c(3,3,5,6,7,6))  
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export
PLS.biplot =  function(X, Y, algorithm=NULL, ax.tickvec.X=NULL, 
                       ax.tickvec.Y=NULL, ... ) 
{
  options(digits=3)
  D = cbind(X,Y)
  X.mean = apply(X, 2, mean)  #column means of X or colMeans(X)
  Y.mean = apply(Y, 2, mean)  #column means of Y or colMeans(Y)
  X.std = apply(X, 2, sd)  #column standard deviations of X 
  Y.std = apply(Y, 2, sd)  #column standard deviations of Y
  X.scal = scale(X, center=TRUE, scale=TRUE)  #scaled X      	
  Y.scal = scale(Y, center=TRUE, scale=TRUE)  #scaled Y    
  
  #biplot	
  A = 2
  main = algorithm(X.scal,Y.scal,A)
  Tmat = main$X.scores
  Pmat = main$X.loadings
  Qmat = main$Y.loadings
  Rmat = main$X.weights.trans	
  Bmat = Rmat %*% t(Qmat)  #(PxM) estimated PLS coefficients matrix	 
  dimnames(Bmat) = list(colnames(X), colnames(Y))	    
  
  #overall quality of approximation		 
  X.hat = ((Tmat %*% t(Pmat)) * rep(X.std,each=nrow(X))) + rep(X.mean,each=nrow(X))
  dimnames(X.hat) = dimnames(X)
  Y.hat = ((Tmat %*% t(Qmat)) * rep(Y.std,each=nrow(Y))) + rep(Y.mean,each=nrow(Y))
  dimnames(Y.hat) = dimnames(Y)
  D.hat = cbind(X.hat, Y.hat)	   
  overall.quality = sum(diag(t(D.hat)%*%D.hat))  /  sum(diag(t(D)%*%D))
  
  #axes predictivity
  axis.pred = (diag(t(D.hat)%*%D.hat)) / (diag(t(D)%*%D))	    
  names(axis.pred) = colnames(D)		 
  
  #biplot points
  Z = rbind(Tmat, Rmat)  #((N+P)x2) biplot points
  dimnames(Z) = list(c(rownames(X),paste("b",1:ncol(X),sep="")), NULL)
  sample.style = biplot.sample.control(2, label=TRUE)
  sample.style$col = c("red", "purple")
  sample.style$pch = c(15,16)
  Gmat = cbind(c(rep(1,nrow(X)),rep(0,ncol(X))), c(rep(0,nrow(X)),rep(1,ncol(X))))
  dimnames(Gmat) = list(NULL, c("samples","coefficients"))
  classes = 1:2	 
  
  #axes direction    
  #X.hat=TP'
  X.axes.direction = (1 / (diag(Pmat%*%t(Pmat)))) * Pmat       	   
  #Y.hat=TQ'
  Y.axes.direction = (1 / (diag(Qmat%*%t(Qmat)))) * Qmat
  #B.pls=RQ'
  B.axes.direction = (1 / (diag(Qmat%*%t(Qmat)))) * Qmat	
  
  #calibration of axes	 	
  z.axes.X = lapply(1:ncol(X), calibrate.axis, X, X.mean, X.std, X.axes.direction, 
                    1:ncol(X), ax.tickvec.X, rep(0,ncol(X)), rep(0,ncol(X)), NULL)  
  z.axes.Y = lapply(1:ncol(Y), calibrate.axis, Y, Y.mean, Y.std, Y.axes.direction, 
                    1:ncol(Y), ax.tickvec.Y, rep(0,ncol(Y)), rep(0,ncol(Y)), NULL)
  z.axes.B = lapply(1:ncol(Bmat), calibrate.axis, Y.scal, rep(0,ncol(Y)), rep(1,ncol(Y)),
                    B.axes.direction, 1:ncol(Bmat), rep(4,ncol(Bmat)), rep(0,ncol(Bmat)), 
                    rep(0,ncol(Bmat)), NULL)	
  z.axes = vector("list", ncol(X)+ncol(Y)+ncol(Bmat))
  for (i in 1:ncol(X)) 
  {
    z.axes[[i]] = z.axes.X[[i]]
  }
  for (i in 1:ncol(Y)) 
  {
    z.axes[[i+ncol(X)]] = z.axes.Y[[i]]
  }
  for (i in 1:ncol(Bmat)) 
  {
    z.axes[[i+ncol(X)+ncol(Y)]] = z.axes.B[[i]]  
  }			  
  
  #biplot axes 
  ax.style = biplot.ax.control(ncol(X)+ncol(Y)+ncol(Bmat), c(colnames(X), colnames(Y), colnames(Bmat)))
  #for X	 
  ax.style$col[1:ncol(X)] = "blue"
  ax.style$label.col[1:ncol(X)] = "blue"
  ax.style$tick.col[1:ncol(X)] = "blue"
  ax.style$tick.label.col[1:ncol(X)] = "blue"  
  #for Y	 
  ax.style$col[ncol(X)+1:ncol(Y)] = "black"
  ax.style$label.col[ncol(X)+1:ncol(Y)] = "black"
  ax.style$tick.col[ncol(X)+1:ncol(Y)] = "black"
  ax.style$tick.label.col[ncol(X)+1:ncol(Y)] = "black"  
  #for B.pls
  ax.style$col[ncol(X)+ncol(Y)+1:ncol(Bmat)] = "black"
  ax.style$label.col[ncol(X)+ncol(Y)+1:ncol(Bmat)] = "black" 
  ax.style$tick.col[ncol(X)+ncol(Y)+1:ncol(Bmat)] = "purple"
  ax.style$tick.label.col[ncol(X)+ncol(Y)+1:ncol(Bmat)] = "purple"   
  draw.biplot(Z, G=Gmat, classes=classes, sample.style=sample.style, z.axes=z.axes, 
              ax.style=ax.style, VV=rbind(X.axes.direction,Y.axes.direction))
  list(overall.quality=overall.quality, axis.pred=axis.pred, D.hat=D.hat, Bmat=Bmat)  	 
}

#' The Partial Least Squares (PLS) biplot with triangles for estimating the Partial Least Squares Regression (PLSR) coefficients 
#' 
#' Takes in a set of predictor variables and a set of response variables and produces a PLS biplot, but with rotated coefficient points.
#' @param X A (NxP) predictor matrix 
#' @param Y A (NxM) response matrix
#' @param algorithm Any of the PLS algorithms ("mod.NIPALS", "mod.KernelPLS_R", "mod.KernelPLS_L", "mod.SIMPLS") 
#' @param ax.tickvec.X tick marker length for each X-variable axis in the PLS biplot
#' @param ax.tickvec.Y tick marker length for each Y-variable axis in the PLS biplot
#' @param base.tri The desired Y-variable axis to use as the base for the triangle
#' @param bi.value The desired rotated coefficient points (bi) to approximate
#' @param ... Other arguments. Currently ignored
#' @return The PLS biplot of D=[X Y] with rotated coefficient points
#' @examples 
#' if(require(pls))
#' data(oliveoil, package="pls")
#' X = as.matrix(oliveoil$chemical, ncol=5) 
#' dimnames(X) = list(paste(c("G1","G2","G3","G4","G5","I1","I2","I3","I4","I5",
#' "S1","S2","S3","S4","S5","S6")),
#' paste(c("Acidity","Peroxide","K232","K270","DK")))   
#' Y = as.matrix(oliveoil$sensory, ncol=6)  
#' dimnames(Y) = list(paste(c("G1","G2","G3","G4","G5","I1","I2","I3","I4","I5",
#' "S1","S2","S3","S4","S5","S6")),
#' paste(c("Yellow","Green","Brown","Glossy","Transp","Syrup")))     
#' #with 1 triangle
#' PLS.biplot.area(X, Y, algorithm=mod.SIMPLS, ax.tickvec.X=c(8,5,5,5,5), 
#' ax.tickvec.Y=c(5,10,5,6,7,10), base.tri=3, bi.value=4) 
#' #with 4 triangles	   
#' PLS.biplot.area(X, Y, algorithm=mod.SIMPLS, ax.tickvec.X=c(8,5,5,5,5), 
#' ax.tickvec.Y=c(5,10,5,6,7,10), base.tri=2, bi.value=c(1,2,3,4,5)) 
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export
PLS.biplot.area =  function(X, Y, algorithm=NULL, ax.tickvec.X=NULL, 
                            ax.tickvec.Y=NULL, base.tri, bi.value, ... ) 
{ 
  #base.tri = desired Y-variable axis to use as the base for the triangle 
  #bi.value = desired rotated coefficient points (bi) to approximate	 
  options(digits=3)
  X.mean = apply(X, 2, mean)  #column means of X or colMeans(X)
  Y.mean = apply(Y, 2, mean)  #column means of Y or colMeans(Y)
  X.std = apply(X, 2, sd)  #column standard deviations of X 
  Y.std = apply(Y, 2, sd)  #column standard deviations of Y	 
  X.scal = scale(X, center=TRUE, scale=TRUE)  #scaled X      	
  Y.scal = scale(Y, center=TRUE, scale=TRUE)  #scaled Y	
  
  #area biplot idea
  #90 degrees rotation matrix for the area triangles
  theta = (pi)/2
  rotation.through.90.degrees = matrix(c(cos(theta), -sin(theta), sin(theta), 
                                         cos(theta)), ncol=2)
  #biplot
  A = 2
  main = algorithm(X.scal,Y.scal,A)
  Tmat = main$X.scores
  Pmat = main$X.loadings
  Qmat = main$Y.loadings
  Rmat = main$X.weights.trans			  
  Bmat = Rmat %*% t(Qmat)  #(PxM) estimated PLS coefficients matrix	
  dimnames(Bmat) = list(colnames(X), colnames(Y))	
  
  #biplot points
  Z = rbind(Tmat, Rmat%*%rotation.through.90.degrees)  #((N+P)x2) biplot points       
  dimnames(Z) = list(c(rownames(X),paste("b",1:ncol(X),sep="")), NULL)
  sample.style = biplot.sample.control(2, label=c(FALSE,TRUE))
  sample.style$col = c("red", "purple")
  sample.style$pch = c(3,16)
  Gmat = cbind(c(rep(1,nrow(X)),rep(0,ncol(X))), c(rep(0,nrow(X)), rep(1,ncol(X))))
  dimnames(Gmat) = list(NULL, c("samples","coefficients"))
  classes = 1:2
  
  #axes direction    
  #X.hat=TP'
  X.axes.direction = (1 / (diag(Pmat%*%t(Pmat)))) * Pmat       	   
  #Y.hat=TQ'
  Y.axes.direction = (1 / (diag(Qmat%*%t(Qmat)))) * Qmat
  
  #calibration of axes	 	
  z.axes.X = lapply(1:ncol(X), calibrate.axis, X, X.mean, X.std, X.axes.direction, 1:ncol(X), 
                    ax.tickvec.X, rep(0,ncol(X)), rep(0,ncol(X)), NULL)      
  z.axes.Y = lapply(1:ncol(Y), calibrate.axis, Y, Y.mean, Y.std, Y.axes.direction, 1:ncol(Y), 
                    ax.tickvec.Y, rep(0,ncol(Y)), rep(0,ncol(Y)), NULL)
  z.axes = vector("list", ncol(X)+ncol(Y))
  for (i in 1:ncol(X)) 
  {
    z.axes[[i]] = z.axes.X[[i]]
  }
  for (i in 1:ncol(Y)) 
  {
    z.axes[[i+ncol(X)]] = z.axes.Y[[i]]
  }					  
  
  #biplot axes  
  ax.style = biplot.ax.control(ncol(X)+ncol(Y), c(colnames(X),colnames(Y)))
  #for X	 
  ax.style$col[1:ncol(X)] = "blue"
  ax.style$label.col[1:ncol(X)] = "blue"
  ax.style$tick.col[1:ncol(X)] = "blue"
  ax.style$tick.label.col[1:ncol(X)] = "blue"  
  #for Y	 
  ax.style$col[ncol(X)+1:ncol(Y)] = "black"
  ax.style$label.col[ncol(X)+1:ncol(Y)] = "black"
  ax.style$tick.col[ncol(X)+1:ncol(Y)] = "black"
  ax.style$tick.label.col[ncol(X)+1:ncol(Y)] = "black"  
  draw.biplot.with.triangles(Z, G=Gmat, classes=classes, sample.style=sample.style, 
                             z.axes=z.axes, ax.style=ax.style, QQ=Y.axes.direction,
                             BB=Rmat%*%rotation.through.90.degrees, base.tri=base.tri, 
                             bi.value=bi.value, VV=rbind(X.axes.direction,Y.axes.direction))     
  list(Bmat=Bmat)  	 
} 

#' The Partial Least Squares (PLS) biplot with no sample points names 
#' 
#' Takes in a set of predictor variables and a set of response variables and produces a PLS biplot, but with no sample points names.
#' @param X A (NxP) predictor matrix 
#' @param Y A (NxM) response matrix
#' @param algorithm Any of the PLS algorithms ("mod.NIPALS", "mod.KernelPLS_R", "mod.KernelPLS_L", "mod.SIMPLS") 
#' @param ax.tickvec.X tick marker length for each X-variable axis in the PLS biplot
#' @param ax.tickvec.Y tick marker length for each Y-variable axis in the PLS biplot
#' @param ... Other arguments. Currently ignored
#' @return The PLS biplot of D=[X Y] with no sample points names
#' @examples  
#' if(require(pls)) 
#' data(oliveoil, package="pls")
#' X = as.matrix(oliveoil$chemical, ncol=5) 
#' dimnames(X) = list(paste(c("G1","G2","G3","G4","G5","I1","I2","I3","I4","I5",
#' "S1","S2","S3","S4","S5","S6")),
#' paste(c("Acidity","Peroxide","K232","K270","DK")))   
#' Y = as.matrix(oliveoil$sensory, ncol=6)  
#' dimnames(Y) = list(paste(c("G1","G2","G3","G4","G5","I1","I2","I3","I4","I5",
#' "S1","S2","S3","S4","S5","S6")),
#' paste(c("Yellow","Green","Brown","Glossy","Transp","Syrup")))     
#' PLS.biplot_no.SN(X, Y, algorithm=mod.SIMPLS, ax.tickvec.X=c(8,5,5,5,5),
#'  ax.tickvec.Y=c(5,8,5,6,9,8))  
#'
#' #cocktail data
#' if(require(SensoMineR))
#' data(cocktail, package="SensoMineR")
#' X3 = as.matrix(compo.cocktail, ncol=4)  
#' Y3 = as.matrix(senso.cocktail, ncol=13) 
#' PLS.biplot_no.SN(X=X3, Y3, algorithm=mod.SIMPLS, ax.tickvec.X=rep(2,ncol(X3)),
#'  ax.tickvec.Y=rep(3,ncol(Y3)))  
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export
PLS.biplot_no.SN = function(X, Y, algorithm=NULL, ax.tickvec.X=NULL, ax.tickvec.Y=NULL, ... ) 
{ 
  options(digits=3)
  D = cbind(X,Y)
  X.mean = apply(X, 2, mean)  #column means of X or colMeans(X)
  Y.mean = apply(Y, 2, mean)  #column means of Y or colMeans(Y)
  X.std = apply(X, 2, sd)  #column standard deviations of X 
  Y.std = apply(Y, 2, sd)  #column standard deviations of Y
  X.scal = scale(X, center=TRUE, scale=TRUE)  #scaled X        
  Y.scal = scale(Y, center=TRUE, scale=TRUE)  #scaled Y    
  
  #biplot	
  A = 2
  main = algorithm(X.scal,Y.scal,A)
  Tmat = main$X.scores
  Pmat = main$X.loadings
  Qmat = main$Y.loadings
  Rmat = main$X.weights.trans	
  Bmat = Rmat %*% t(Qmat)  #(PxM) estimated PLS coefficients matrix	 
  dimnames(Bmat) = list(colnames(X), colnames(Y))	    
  
  #overall quality of approximation		 
  X.hat = ((Tmat %*% t(Pmat)) * rep(X.std,each=nrow(X))) + rep(X.mean,each=nrow(X))
  dimnames(X.hat) = dimnames(X)
  Y.hat = ((Tmat %*% t(Qmat)) * rep(Y.std,each=nrow(Y))) + rep(Y.mean,each=nrow(Y))
  dimnames(Y.hat) = dimnames(Y)
  D.hat = cbind(X.hat, Y.hat)	   
  overall.quality = sum(diag(t(D.hat)%*%D.hat))  /  sum(diag(t(D)%*%D))
  
  #axes predictivity
  axis.pred = (diag(t(D.hat)%*%D.hat)) / (diag(t(D)%*%D))	    
  names(axis.pred) = colnames(D)		 	 
  
  #biplot points
  Z = rbind(Tmat, Rmat) #((N+P)x2) biplot points
  dimnames(Z) = list(c(rownames(X),paste("b",1:ncol(X),sep="")), NULL)
  sample.style = biplot.sample.control(2, label=c(FALSE,TRUE)) 
  #No sample names in the biplot
  sample.style$col = c("red", "purple")
  sample.style$pch = c(15,16)
  Gmat = cbind(c(rep(1,nrow(X)),rep(0,ncol(X))), c(rep(0,nrow(X)),rep(1,ncol(X))))
  dimnames(Gmat) = list(NULL, c("samples","coefficients"))
  classes = 1:2	 
  
  #axes direction    
  #X.hat=TP'
  X.axes.direction = (1 / (diag(Pmat%*%t(Pmat)))) * Pmat       	   
  #Y.hat=TQ'
  Y.axes.direction = (1 / (diag(Qmat%*%t(Qmat)))) * Qmat
  #B.pls=RQ'
  B.axes.direction = (1 / (diag(Qmat%*%t(Qmat)))) * Qmat	
  
  #calibration of axes	 	
  z.axes.X = lapply(1:ncol(X), calibrate.axis, X, X.mean, X.std, X.axes.direction, 1:ncol(X), 
                    ax.tickvec.X, rep(0,ncol(X)), rep(0,ncol(X)), NULL)      
  z.axes.Y = lapply(1:ncol(Y), calibrate.axis, Y, Y.mean, Y.std, Y.axes.direction, 1:ncol(Y), 
                    ax.tickvec.Y, rep(0,ncol(Y)), rep(0,ncol(Y)), NULL)
  z.axes.B = lapply(1:ncol(Bmat), calibrate.axis, Y.scal, rep(0,ncol(Y)), rep(1,ncol(Y)), 
                    B.axes.direction, 1:ncol(Bmat), rep(4,ncol(Bmat)), rep(0,ncol(Bmat)), 
                    rep(0,ncol(Bmat)), NULL)	
  z.axes = vector("list", ncol(X)+ncol(Y)+ncol(Bmat))
  for (i in 1:ncol(X)) 
  {
    z.axes[[i]] = z.axes.X[[i]]
  }
  for (i in 1:ncol(Y)) 
  {
    z.axes[[i+ncol(X)]] = z.axes.Y[[i]]
  }
  for (i in 1:ncol(Bmat)) 
  {
    z.axes[[i+ncol(X)+ncol(Y)]] = z.axes.B[[i]]  
  }			  
  
  #biplot axes 
  ax.style = biplot.ax.control(ncol(X)+ncol(Y)+ncol(Bmat), c(colnames(X),colnames(Y),colnames(Bmat)))
  #for X	 
  ax.style$col[1:ncol(X)] = "blue"
  ax.style$label.col[1:ncol(X)] = "blue"
  ax.style$tick.col[1:ncol(X)] = "blue"
  ax.style$tick.label.col[1:ncol(X)] = "blue"  
  #for Y	 
  ax.style$col[ncol(X)+1:ncol(Y)] = "black"
  ax.style$label.col[ncol(X)+1:ncol(Y)] = "black"
  ax.style$tick.col[ncol(X)+1:ncol(Y)] = "black"
  ax.style$tick.label.col[ncol(X)+1:ncol(Y)] = "black"  
  #for B.pls
  ax.style$col[ncol(X)+ncol(Y)+1:ncol(Bmat)] = "black"
  ax.style$label.col[ncol(X)+ncol(Y)+1:ncol(Bmat)] = "black" 
  ax.style$tick.col[ncol(X)+ncol(Y)+1:ncol(Bmat)] = "purple"
  ax.style$tick.label.col[ncol(X)+ncol(Y)+1:ncol(Bmat)] = "purple"   
  draw.biplot(Z, G=Gmat, classes=classes, sample.style=sample.style, z.axes=z.axes, ax.style=ax.style, VV=rbind(X.axes.direction,Y.axes.direction))  
  list(overall.quality=overall.quality, axis.pred=axis.pred, D.hat=D.hat, Bmat=Bmat)  	 
}

#' The Partial Least Squares (PLS) biplot with the labels of the samples, coefficient points and tick markers excluded
#' 
#' Takes in a set of predictor variables and a set of response variables and produces a PLS biplot, but with no sample points names.
#' @param X A (NxP) predictor matrix 
#' @param Y A (NxM) response matrix
#' @param algorithm Any of the PLS algorithms ("mod.NIPALS", "mod.KernelPLS_R", "mod.KernelPLS_L", "mod.SIMPLS") 
#' @param ax.tickvec.X tick marker length for each X-variable axis in the PLS biplot
#' @param ax.tickvec.Y tick marker length for each Y-variable axis in the PLS biplot
#' @param ... Other arguments. Currently ignored
#' @return The PLS biplot of D=[X Y] with no sample points names but with fainted tick markers and labels
#' @examples  
#' if(require(pls))
#' data(oliveoil, package="pls")
#' X = as.matrix(oliveoil$chemical, ncol=5) 
#' dimnames(X) = list(paste(c("G1","G2","G3","G4","G5","I1","I2","I3","I4","I5",
#' "S1","S2","S3","S4","S5","S6")),
#' paste(c("Acidity","Peroxide","K232","K270","DK")))   
#' Y = as.matrix(oliveoil$sensory, ncol=6)  
#' dimnames(Y) = list(paste(c("G1","G2","G3","G4","G5","I1","I2","I3","I4","I5",
#' "S1","S2","S3","S4","S5","S6")),
#' paste(c("Yellow","Green","Brown","Glossy","Transp","Syrup")))     
#' PLS.biplot_no_labels(X, Y, algorithm=mod.SIMPLS, ax.tickvec.X=c(8,5,5,5,5), 
#' ax.tickvec.Y=c(5,8,5,6,9,8))  
#'
#' #cocktail data
#' if(require(SensoMineR))
#' data(cocktail, package="SensoMineR")
#' X3 = as.matrix(compo.cocktail, ncol=4)  
#' Y3 = as.matrix(senso.cocktail, ncol=13) 
#' PLS.biplot_no_labels(X=X3, Y3, algorithm=mod.SIMPLS, ax.tickvec.X=rep(2,ncol(X3)), 
#' ax.tickvec.Y=rep(3,ncol(Y3)))  
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export
PLS.biplot_no_labels = function(X, Y, algorithm=NULL, ax.tickvec.X=NULL, ax.tickvec.Y=NULL, ... ) 
{ 
  options(digits=3)
  D = cbind(X,Y)
  X.mean = apply(X, 2, mean)  #column means of X or colMeans(X)
  Y.mean = apply(Y, 2, mean)  #column means of Y or colMeans(Y)
  X.std = apply(X, 2, sd)  #column standard deviations of X 
  Y.std = apply(Y, 2, sd)  #column standard deviations of Y
  X.scal = scale(X, center=TRUE, scale=TRUE)  #scaled X        
  Y.scal = scale(Y, center=TRUE, scale=TRUE)  #scaled Y    
  
  #biplot	
  A = 2
  main = algorithm(X.scal,Y.scal,A)
  Tmat = main$X.scores
  Pmat = main$X.loadings
  Qmat = main$Y.loadings
  Rmat = main$X.weights.trans	
  Bmat = Rmat %*% t(Qmat)  #(PxM) estimated PLS coefficients matrix	 
  dimnames(Bmat) = list(colnames(X), colnames(Y))	    
  
  #overall quality of approximation		 
  X.hat = ((Tmat %*% t(Pmat)) * rep(X.std,each=nrow(X))) + rep(X.mean,each=nrow(X))
  dimnames(X.hat) = dimnames(X)
  Y.hat = ((Tmat %*% t(Qmat)) * rep(Y.std,each=nrow(Y))) + rep(Y.mean,each=nrow(Y))
  dimnames(Y.hat) = dimnames(Y)
  D.hat = cbind(X.hat, Y.hat)	   
  overall.quality = sum(diag(t(D.hat)%*%D.hat))  /  sum(diag(t(D)%*%D))
  
  #axes predictivity
  axis.pred = (diag(t(D.hat)%*%D.hat)) / (diag(t(D)%*%D))	    
  names(axis.pred) = colnames(D)		 	 
  
  #biplot points
  Z = rbind(Tmat, Rmat) #((N+P)x2) biplot points
  dimnames(Z) = list(c(rownames(X),paste("b",1:ncol(X),sep="")), NULL)
  sample.style = biplot.sample.control(2, label=FALSE) 
  #No samples and coefficients names in the biplot
  sample.style$col = c("red", "purple")
  sample.style$pch = c(15,16)
  Gmat = cbind(c(rep(1,nrow(X)),rep(0,ncol(X))), c(rep(0,nrow(X)),rep(1,ncol(X))))
  dimnames(Gmat) = list(NULL, c("samples","coefficients"))
  classes = 1:2	 
  
  #axes direction    
  #X.hat=TP'
  X.axes.direction = (1 / (diag(Pmat%*%t(Pmat)))) * Pmat       	   
  #Y.hat=TQ'
  Y.axes.direction = (1 / (diag(Qmat%*%t(Qmat)))) * Qmat
  #B.pls=RQ'
  B.axes.direction = (1 / (diag(Qmat%*%t(Qmat)))) * Qmat	
  
  #calibration of axes	 	
  z.axes.X = lapply(1:ncol(X), calibrate.axis, X, X.mean, X.std, X.axes.direction, 1:ncol(X), 
                    ax.tickvec.X, rep(0,ncol(X)), rep(0,ncol(X)), NULL)      
  z.axes.Y = lapply(1:ncol(Y), calibrate.axis, Y, Y.mean, Y.std, Y.axes.direction, 1:ncol(Y), 
                    ax.tickvec.Y, rep(0,ncol(Y)), rep(0,ncol(Y)), NULL)
  z.axes.B = lapply(1:ncol(Bmat), calibrate.axis, Y.scal, rep(0,ncol(Y)), rep(1,ncol(Y)), 
                    B.axes.direction, 1:ncol(Bmat), rep(0,ncol(Bmat)), rep(0,ncol(Bmat)), 
                    rep(0,ncol(Bmat)), NULL)	
  z.axes = vector("list", ncol(X)+ncol(Y)+ncol(Bmat))
  for (i in 1:ncol(X)) 
  {
    z.axes[[i]] = z.axes.X[[i]]
  }
  for (i in 1:ncol(Y)) 
  {
    z.axes[[i+ncol(X)]] = z.axes.Y[[i]]
  }
  for (i in 1:ncol(Bmat)) 
  {
    z.axes[[i+ncol(X)+ncol(Y)]] = z.axes.B[[i]]  
  }			  
  
  #biplot axes 
  ax.style = biplot.ax.control(ncol(X)+ncol(Y)+ncol(Bmat), c(colnames(X),colnames(Y),colnames(Bmat)))
  ax.style$tick.col[1:ncol(X)] = "azure"  
  ax.style$tick.label.col[1:ncol(X)] = "azure" 
  ax.style$tick.col[ncol(X)+1:ncol(Y)] = "azure"  
  ax.style$tick.label.col[ncol(X)+1:ncol(Y)] = "azure" 
  ax.style$tick.col[ncol(X)+ncol(Y)+1:ncol(Bmat)] = "azure" 
  ax.style$tick.label.col[ncol(X)+ncol(Y)+1:ncol(Bmat)] = "azure"
  #for X	 
  ax.style$col[1:ncol(X)] = "blue"
  ax.style$label.col[1:ncol(X)] = "blue"   
  #for Y	 
  ax.style$col[ncol(X)+1:ncol(Y)] = "black"
  ax.style$label.col[ncol(X)+1:ncol(Y)] = "black" 
  #for B.pls
  ax.style$col[ncol(X)+ncol(Y)+1:ncol(Bmat)] = "black"
  ax.style$label.col[ncol(X)+ncol(Y)+1:ncol(Bmat)] = "black"     
  draw.biplot(Z, G=Gmat, classes=classes, sample.style=sample.style, z.axes=z.axes, ax.style=ax.style, VV=rbind(X.axes.direction,Y.axes.direction))  
  list(overall.quality=overall.quality, axis.pred=axis.pred, D.hat=D.hat, Bmat=Bmat)  	 
}  


# --------------------------------------------------------------------------------------------------------------------


#A.6  Chapter  6: PLS biplots for Generalized Linear Models 
# Partial Least Squares-Generalized Linear Models (PLS-GLMs) algorithms.
# PLS biplot for the (univariate) Poisson and Binomial PLS-GLMs. 
# RGui(32-bit).

#' Partial Least Squares-Generalized Linear Model (PLS-GLM) algorithm 
#' 
#' Takes in a set of predictor variables and a set of response variables and gives the PLS-GLM parameters.
#' @param X A (NxP) predictor matrix 
#' @param y A (Nx1) Poisson-distributed response vector 
#' @param A The number of PLS components 
#' @param PLS_algorithm The mod.NIPALS algorithm 
#' @param eps Cut off value for convergence step 
#' @param ... Other arguments. Currently ignored
#' @return The PLS-GLM parameters of D=[X y]
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export
PLS.GLM = function(X, y, A, PLS_algorithm=NULL, eps=1e-04, ...)
{
  #For M=1
  #y is Poisson distributed
  #A=initial number of components    
  #1. Initialization
  eta = log(y+2)  #link function 
  N = nrow(X)
  P = ncol(X)
  one = matrix(rep(1,each=N), nrow=N)  #(Nx1) vector of ones
  Vmat = diag(c(exp(eta)), ncol=N, nrow=N)  #(NxN) diagonal weight matrix for the weighted least squares
  z.weighted.mean = ( (one%*%t(one)%*%Vmat%*%eta)/N ) 
  z0 = eta - z.weighted.mean
  X.weighted.mean = ( (one%*%t(one)%*%Vmat%*%X)/N )
  X0 = X - X.weighted.mean
  W.star =  P.star = matrix(0, nrow=P, ncol=A)  #(PxA) (deflated) X.weights and X.loadings
  #matrices	
  R.star = matrix(0, nrow=P, ncol=A)  #(PxA) (transformed) X.weights matrix 											  
  q.star_vec = matrix(0, nrow=1, ncol=A)  #(1xA) y.loadings vector
  T.star = matrix(0, nrow=N, ncol=A)  #(NxA) X.scores matrix 
  main0 = PLS_algorithm(X,y,A)
  b.pls_glm = main0$X.weights.trans %*% t(main0$Y.loadings)
  #2. Weighted PLS 
  repeat
  {
    X.a = X0
    z.a = z0
    for (a in 1:A)
    {
      w.a = t(X.a) %*% Vmat %*% z.a  / as.numeric(sqrt(t(t(X.a)%*%Vmat%*%z.a)%*%t(X.a)%*%Vmat%*%z.a))  #(Px1) X.weight vector
      t.a = X.a %*% w.a  / as.numeric(sqrt(t(X.a%*%w.a)%*%X.a%*%w.a))     #(Nx1) X.score vector
      q.a = t(z.a) %*% Vmat %*% t.a  /  as.numeric(t(t.a)%*%Vmat%*%t.a)    #(Mx1) y.loading vector
      p.a = t(X.a) %*% Vmat %*% t.a  /  as.numeric(t(t.a)%*%Vmat%*%t.a)   #(Px1) X.loading vector
      r.a = w.a%*%solve(t(p.a)%*%w.a)  #(Px1) (transformed) X.weight vector
      t.a = X0 %*% r.a / as.numeric(sqrt(t(X0%*%r.a)%*%X0%*%r.a))
      
      #update X.a and z.a 
      X.a = X.a - t.a%*%t(p.a)
      z.a = z.a - t.a%*%t(q.a)  
      
      W.star[, a] = w.a   #(PxA) X.weights matrix       
      P.star[, a] = p.a   #(PxA) X.loadings matrix
      q.star_vec[, a] = q.a   #(1xA) y.loadings vector
      R.star[, a] = r.a    #(PxA) transformed X.weights matrix 
      T.star[, a] = t.a    #(NxA) X.scores matrix
      
      dimnames(P.star) = dimnames(W.star) = dimnames(R.star) = list(colnames(X), paste("Comp", 1:A, seq=""))
      dimnames(q.star_vec) = list(colnames(y), paste("Comp", 1:A, seq=""))
      dimnames(T.star) = list(rownames(X), paste("Comp", 1:A, seq=""))
      #t(T.star)%*%Vmat%*%T.star  #non-diagonal elements of 0
      #t(T.star)%*%T.star  #diagonal elements of 1
    }     
    #3. Update eta 
    eta = T.star%*%t(q.star_vec) +  z.weighted.mean 	
    #4. Update Vmat and z0
    Vmat = diag(c(exp(eta)), ncol=N, nrow=N) 
    z0 = eta + diag(c(1/exp(eta)), ncol=N, nrow=N) %*% (y-exp(eta))  #the linearized form     
    #5-6. Check that changes are sufficiently small, else re-run from step 2  
    b.pls_glm.old = b.pls_glm
    b.pls_glm = R.star %*% t(q.star_vec)  #(Px1) PLS-GLM coefficient vector 
    if ( (sum((b.pls_glm) - (b.pls_glm.old)/(b.pls_glm))) < eps )           		   
      break
    else
    {
      W.star = W.star   #(PxA) X.weights matrix       
      P.star = P.star   #(PxA) X.loadings matrix
      q.star_vec = q.star_vec   #(1xA) y.loadings vector
      R.star = R.star    #(PxA) transformed X.weights matrix 
      T.star = T.star    #(NxA) X.scores matrix
      b.pls_glm = R.star %*% t(q.star_vec)  #(Px1) PLS-GLM coefficient vector 
    }
  }		 
  list(X.weights=W.star, X.loadings=P.star, y.loadings=q.star_vec, X.weights.trans=R.star, 
       X.scores=T.star, b.pls_glm_vec=b.pls_glm)
} 

#' Partial Least Squares-Generalized Linear Model (PLS-GLM) fitted using the SIMPLS algorithm 
#' 
#' Takes in a set of predictor variables and a set of response variables and gives the PLS-GLM parameters.
#' @param X A (NxP) predictor matrix 
#' @param y A (Nx1) Poisson-distributed response vector 
#' @param A The number of PLS components 
#' @param PLS_algorithm The mod.SIMPLS algorithm  
#' @param eps Cut off value for convergence step
#' @param ... Other arguments. Currently ignored
#' @return The PLS-GLM parameters of D=[X y]
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export 
PLS.GLM_SIMPLS = function(X, y, A, PLS_algorithm=NULL, eps=1e-04, ...)
{
  #For M=1
  #A=initial number of components 
  #1. Initialization
  eta = log(y+2)  #link function 
  N = nrow(X)
  P = ncol(X)
  one = matrix(rep(1,each=N), nrow=N)  #(Nx1) vector of ones
  Vmat = diag(c(exp(eta)), ncol=N, nrow=N)  #(NxN) diagonal weight matrix for the weighted least squares
  z.weighted.mean = ( (one%*%t(one)%*%Vmat%*%eta)/N ) 
  z0 = eta - z.weighted.mean
  X.weighted.mean = ( (one%*%t(one)%*%Vmat%*%X)/N )
  X0 = X - X.weighted.mean
  P.star = matrix(0, nrow=P, ncol=A)  #(PxA) (deflated) X.loadings matrices  
  R.star = matrix(0, nrow=P, ncol=A)  #(PxA) (transformed) X.weights matrix 											  
  q.star_vec = matrix(0, nrow=1, ncol=A)  #(1xA) y.loadings vector
  T.star = matrix(0, nrow=N, ncol=A)  #(NxA) X.scores matrix 
  main0 = PLS_algorithm(X,y,A)
  b.pls_glm = main0$X.weights.trans %*% t(main0$Y.loadings)
  #2. Weighted PLS 	  
  repeat
  {
    S = t(X0) %*% Vmat %*% z0 
    for (a in 1:A)
    {
      c.a = svd(S)$v[,1]  #(Mx1) y.weight
      r.a = svd(S)$u[,1]  #(Px1) X.weight         
      t.a = X0 %*% r.a   #(Nx1) X.score
      t.a.norm = as.numeric(sqrt(t(t.a)%*%t.a))
      t.a = t.a  / t.a.norm   #normalized
      q.a = t(z0) %*% Vmat %*% t.a  /  as.numeric(t(t.a)%*%Vmat%*%t.a)    #(Mx1) y-loading vector
      p.a = t(X0) %*% Vmat %*% t.a  /  as.numeric(t(t.a)%*%Vmat%*%t.a)   #(Px1) X-loading vector
      
      #update S
      S = S - p.a%*%solve(t(p.a)%*%p.a)%*%t(p.a)%*%S	     			   
      
      P.star[, a] = p.a   #(PxA) X.loadings matrix
      q.star_vec[, a] = q.a   #(1xA) y.loadings vector
      R.star[, a] = r.a    #(PxA) transformed X.weights matrix 
      T.star[, a] = t.a    #(NxA) X.scores matrix
      
      dimnames(P.star) = dimnames(R.star) = list(colnames(X), paste("Comp", 1:A, seq=""))
      dimnames(q.star_vec) = list(colnames(y), paste("Comp", 1:A, seq=""))
      dimnames(T.star) = list(rownames(X), paste("Comp", 1:A, seq=""))
      #t(T.star)%*%Vmat%*%T.star  #non-diagonal elements of 0
      #t(T.star)%*%T.star  #diagonal elements of 1
    }     
    #3. Update eta 
    eta = T.star%*%t(q.star_vec) +  z.weighted.mean 	
    #4. Update Vmat and z0
    Vmat = diag(c(exp(eta)), ncol=N, nrow=N) 
    z0 = eta + diag(c(1/exp(eta)), ncol=N, nrow=N) %*% (y-exp(eta))  #the linearized form     
    #5-6. Check that changes are sufficiently small, else re-run from step 2  
    b.pls_glm.old = b.pls_glm
    b.pls_glm = R.star %*% t(q.star_vec)  #(Px1) PLS-GLM coefficient vector 
    if ( (sum((b.pls_glm) - (b.pls_glm.old)/(b.pls_glm))) < eps)           		   
      break
    else
    {                             
      P.star = P.star   #(PxA) X.loadings matrix
      q.star_vec = q.star_vec   #(1xA) y.loadings vector
      R.star = R.star    #(PxA) transformed X.weights matrix 
      T.star = T.star    #(NxA) X.scores matrix
      b.pls_glm = R.star %*% t(q.star_vec)  #(Px1) PLS-GLM coefficient vector 
    }
  }		 
  list(X.loadings=P.star, y.loadings=q.star_vec, X.weights.trans=R.star, 
       X.scores=T.star, b.pls_glm_vec=b.pls_glm)
} 

#' Partial Least Squares-Generalized Linear Model (PLS-GLM) fitted for Binomial y 
#' 
#' Takes in a set of predictor variables and a set of response variables and gives the PLS-GLM parameters.
#' @param X A (NxP) predictor matrix 
#' @param y A (Nx1) Binomial-distributed response vector 
#' @param A The number of PLS components 
#' @param PLS_algorithm The mod.NIPALS algorithm  
#' @param eps Cut off value for convergence step
#' @param ... Other arguments. Currently ignored
#' @return The PLS-GLM parameters of D=[X y]
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export
PLS.binomial.GLM = function(X, y, A, PLS_algorithm=NULL, eps=1e-04, ...)
{
  #For M=1 
  #y is Binomial distributed
  #A=initial number of components   
  #1. Initialization  
  main0 = PLS_algorithm(X,y,A)
  b.pls_glm = main0$X.weights.trans %*% t(main0$Y.loadings)
  eta = X %*% b.pls_glm  #link function 
  N = nrow(X)
  P = ncol(X)
  one = matrix(rep(1,each=N), nrow=N)  #(Nx1) vector of ones
  pii = exp(X%*%b.pls_glm) / (1+exp(X%*%b.pls_glm))
  Vmat = diag(c(N*pii*(1-pii)), ncol=N, nrow=N)  #(NxN) diagonal weight matrix for the weighted least squares
  z.weighted.mean = ( (one%*%t(one)%*%Vmat%*%eta)/N ) 
  z0 = eta - z.weighted.mean
  X.weighted.mean = ( (one%*%t(one)%*%Vmat%*%X)/N )
  X0 = X - X.weighted.mean
  W.star =  P.star = matrix(0, nrow=P, ncol=A)  #(PxA) (deflated) X.weights and X.loadings
  #matrices	
  R.star = matrix(0, nrow=P, ncol=A)  #(PxA) (transformed) X.weights matrix 											  
  q.star_vec = matrix(0, nrow=1, ncol=A)  #(1xA) y.loadings vector
  T.star = matrix(0, nrow=N, ncol=A)  #(NxA) X.scores matrix 
  
  #2. Weighted PLS     
  repeat
  {
    X.a = X0
    z.a = z0
    for (a in 1:A)
    {
      w.a = t(X.a) %*% Vmat %*% z.a  / as.numeric(sqrt(t(t(X.a)%*%Vmat%*%z.a)%*%t(X.a)%*%Vmat%*%z.a))    #(Px1) X-weight vector
      t.a = X.a %*% w.a  / as.numeric(sqrt(t(X.a%*%w.a)%*%X.a%*%w.a))     #(Nx1) X-score vector
      q.a = t(z.a) %*% Vmat %*% t.a  /  as.numeric(t(t.a)%*%Vmat%*%t.a)    #(Mx1) Y-loading vector
      p.a = t(X.a) %*% Vmat %*% t.a  /  as.numeric(t(t.a)%*%Vmat%*%t.a)   #(Px1) X-loading vector
      r.a = w.a%*%solve(t(p.a)%*%w.a)  #(Px1) (transformed) X-weight vector
      t.a = X0 %*% r.a / as.numeric(sqrt(t(X0 %*% r.a)%*%X0 %*% r.a))
      
      #update X.a and z.a 
      X.a = X.a - t.a%*%t(p.a)
      z.a = z.a - t.a%*%t(q.a)  
      
      W.star[, a] = w.a   #(PxA) X.weights matrix       
      P.star[, a] = p.a   #(PxA) X.loadings matrix
      q.star_vec[, a] = q.a   #(1xA) y.loadings vector
      R.star[, a] = r.a    #(PxA) transformed X.weights matrix 
      T.star[, a] = t.a    #(NxA) X.scores matrix
      
      dimnames(P.star) = dimnames(W.star) = dimnames(R.star) = list(colnames(X), paste("Comp", 1:A, seq=""))
      dimnames(q.star_vec) = list(colnames(y), paste("Comp", 1:A, seq=""))
      dimnames(T.star) = list(rownames(X), paste("Comp", 1:A, seq=""))
      #t(T.star)%*%Vmat%*%T.star  #non-diagonal elements of 0
      #t(T.star)%*%T.star  #diagonal elements of 1
    }     
    #3. Update eta 
    eta = T.star%*%t(q.star_vec) +  z.weighted.mean 
    
    #4. Update Vmat and z0
    pii = exp(X0%*%R.star%*%t(q.star_vec)) / (1+exp(X0%*%R.star%*%t(q.star_vec)))
    Vmat = diag(c(N*pii*(1-pii)), ncol=N, nrow=N) 
    z0 = eta + ((y-(N*pii))/(N*pii*(1-pii)))  #the linearized form 
    
    #5-6. Check that changes are sufficiently small, else re-run from step 2  
    b.pls_glm.old = b.pls_glm
    b.pls_glm = R.star %*% t(q.star_vec)  #(Px1) PLS-GLM coefficient vector 
    if ( (sum((b.pls_glm) - (b.pls_glm.old)/(b.pls_glm))) < eps )           		   
      break
    else
    {
      W.star = W.star   #(PxA) X.weights matrix       
      P.star = P.star   #(PxA) X.loadings matrix
      q.star_vec = q.star_vec   #(1xA) y.loadings vector
      R.star = R.star    #(PxA) transformed X.weights matrix 
      T.star = T.star    #(NxA) X.scores matrix
      b.pls_glm = R.star %*% t(q.star_vec)  #(Px1) PLS-GLM coefficient vector 
    }
  }		 
  list(X.weights=W.star, X.loadings=P.star, y.loadings=q.star_vec, X.weights.trans=R.star, 
       X.scores=T.star, b.pls_glm_vec=b.pls_glm)
}  

# PLS-GLM biplots
#' The Partial Least Squares (PLS) biplot for Generalized Linear Model (GLM) 
#' 
#' Takes in a set of predictor variables and a set of response variables and produces a PLS biplot for the (univariate) GLMs.
#' @param X A (NxP) predictor matrix 
#' @param y A (Nx1) response vector 
#' @param algorithm The PLS.GLM algorithm 
#' @param ax.tickvec.X tick marker length for each X-variable axis in the biplot
#' @param ax.tickvec.y tick marker length for the y-variable axis in the biplot
#' @param ax.tickvec.b (purple) tick marker length for the y-variable axis in the biplot
#' @param ... Other arguments. Currently ignored
#' @return The PLS biplot of a GLM of D=[X y] with some parameters
#' @examples
#' if(require(robustbase))
#' possum.mat #data matrix
#' y = as.matrix(possum.mat[,1], ncol=1) 
#' dimnames(y) = list(paste("S", 1:nrow(possum.mat), seq=""), "Diversity")
#' X = as.matrix(possum.mat[,2:14], ncol=13) 
#' dimnames(X) = list(paste("S", 1:nrow(possum.mat), seq=""), colnames(possum.mat[,2:14]))  
#' PLS.GLM.biplot(X, y, algorithm=PLS.GLM, ax.tickvec.X=rep(5,ncol(X)),
#'  ax.tickvec.y=10, ax.tickvec.b=7)  
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export
PLS.GLM.biplot = function(X, y, algorithm=NULL, ax.tickvec.X=NULL, ax.tickvec.y=NULL, ax.tickvec.b=NULL, ... ) 
{    
  options(digits=3)
  D = cbind(X,y)
  N = nrow(D)
  X.mean = apply(X, 2, mean)  #column means of X or colMeans(X)
  y.mean = mean(y)  #mean of y 
  X.std = apply(X, 2, sd)  #column standard deviations of X 
  y.std = sd(y)  #standard deviation of y   
  X.scal = scale(X, center=TRUE, scale=TRUE)  #scaled X
  y.scal = scale(y, center=TRUE, scale=TRUE)  #scaled y  
  
  #biplot 
  A = 2
  main = algorithm(X.scal,y.scal,A,PLS_algorithm=mod.NIPALS,eps=1e-04)
  Tmat = main$X.scores
  Pmat = main$X.loadings
  qvec = main$y.loadings
  Rmat = main$X.weights.trans	
  bvec = Rmat %*% t(qvec)  #(Px1) estimated PLS-GLM coefficient vector 	 
  dimnames(bvec) = list(colnames(X), colnames(y))	
  
  #approximation		 
  X.hat = ( (Tmat %*% t(Pmat)) * rep(X.std,each=N) ) + rep(X.mean,each=N)
  dimnames(X.hat) = dimnames(X)
  y.hat = ( (Tmat %*% t(qvec)) * rep(y.std,each=N) ) + rep(y.mean,each=N)
  dimnames(y.hat) = list(rownames(y), paste("Expected_",colnames(y),sep=""))
  D.hat = cbind(X.hat, y.hat)	  
  
  #biplot points
  Z = rbind(Tmat, Rmat)  #((N+P)x2) biplot points
  dimnames(Z) = list(c(rownames(X),paste("b",1:ncol(X),sep="")), NULL)
  sample.style = biplot.sample.control(2, label=TRUE)
  sample.style$col = c("red", "purple")
  sample.style$pch = c(15,16)
  Gmat = cbind(c(rep(1,nrow(X)),rep(0,ncol(X))), c(rep(0,nrow(X)),rep(1,ncol(X))))
  dimnames(Gmat) = list(NULL, c("samples","coefficients"))
  classes = 1:2	 
  
  #axes direction    
  #X.hat=TP'
  X.axes.direction = (1 / (diag(Pmat%*%t(Pmat)))) * Pmat       	   
  #y.hat=Tq
  qvec = matrix(qvec, nrow=1)  #since for M=1, j = 1:1
  y.axes.direction = (1 / (diag(qvec%*%t(qvec)))) * qvec
  #b.pls-glm=RQ'
  b.axes.direction = (1 / (diag(qvec%*%t(qvec)))) * qvec	
  
  #calibration of axes	 	
  z.axes.X = lapply(1:ncol(X), calibrate.axis, X, X.mean, X.std, X.axes.direction, 
                    1:ncol(X), ax.tickvec.X, rep(0,ncol(X)), rep(0,ncol(X)), NULL)  
  z.axes.y = lapply(1, calibrate.axis, y, y.mean, y.std, y.axes.direction, 
                    1, ax.tickvec.y, rep(0,1), rep(0,1), NULL)
  z.axes.b = lapply(1, calibrate.axis, y.scal, rep(0,1), rep(1,1), b.axes.direction, 
                    1, ax.tickvec.b, rep(0,1), rep(0,1), NULL)		
  z.axes = vector("list", ncol(X)+ncol(y)+ncol(bvec))
  for (i in 1:ncol(X)) 
  {
    z.axes[[i]] = z.axes.X[[i]]
  }
  for (i in 1:ncol(y)) 
  {
    z.axes[[i+ncol(X)]] = z.axes.y[[i]]
  }
  for (i in 1:ncol(bvec)) 
  {
    z.axes[[i+ncol(X)+ncol(y)]] = z.axes.b[[i]]  
  }			  
  
  #biplot axes 
  ax.style = biplot.ax.control(ncol(X)+ncol(y)+ncol(bvec), c(colnames(X), colnames(y), colnames(bvec)))
  #for X	 
  ax.style$col[1:ncol(X)] = "blue"
  ax.style$label.col[1:ncol(X)] = "blue"
  ax.style$tick.col[1:ncol(X)] = "blue"
  ax.style$tick.label.col[1:ncol(X)] = "blue"  
  #for y	 
  ax.style$col[ncol(X)+1:ncol(y)] = "black"
  ax.style$label.col[ncol(X)+1:ncol(y)] = "black"
  ax.style$tick.col[ncol(X)+1:ncol(y)] = "black"
  ax.style$tick.label.col[ncol(X)+1:ncol(y)] = "black"  
  #for b.pls-glm
  ax.style$col[ncol(X)+ncol(y)+1:ncol(bvec)] = "black"
  ax.style$label.col[ncol(X)+ncol(y)+1:ncol(bvec)] = "black" 
  ax.style$tick.col[ncol(X)+ncol(y)+1:ncol(bvec)] = "purple"
  ax.style$tick.label.col[ncol(X)+ncol(y)+1:ncol(bvec)] = "purple"   
  draw.biplot(Z, G=Gmat, classes=classes, sample.style=sample.style, z.axes=z.axes, ax.style=ax.style, VV=rbind(X.axes.direction,y.axes.direction)) 
  
  list(D.hat=D.hat, bvec=bvec)  	 
}

#' The Partial Least Squares (PLS) biplot for Generalized Linear Model (GLM) fitted using the SIMPLS algorithm   
#' 
#' Takes in a set of predictor variables and a set of response variables and produces a PLS biplot for the (univariate) GLMs.
#' @param X A (NxP) predictor matrix 
#' @param y A (Nx1) response vector 
#' @param algorithm The PLS.GLM_SIMPLS algorithm
#' @param ax.tickvec.X tick marker length for each X-variable axis in the biplot
#' @param ax.tickvec.y tick marker length for the y-variable axis in the biplot
#' @param ax.tickvec.b (purple) tick marker length for the y-variable axis in the biplot
#' @param ... Other arguments. Currently ignored
#' @return The PLS biplot of a GLM (fitted using the SIMPLS algorithm) of D=[X y] with some parameters
#' @examples
#' if(require(robustbase))
#' possum.mat   
#' y = as.matrix(possum.mat[,1], ncol=1) 
#' dimnames(y) = list(paste("S", 1:nrow(possum.mat), seq=""), "Diversity")
#' X = as.matrix(possum.mat[,2:14], ncol=13) 
#' dimnames(X) = list(paste("S", 1:nrow(possum.mat), seq=""), colnames(possum.mat[,2:14]))
#' PLS.GLM.biplot_SIMPLS(X, y, algorithm=PLS.GLM, 
#' ax.tickvec.X=rep(5,ncol(X)), ax.tickvec.y=10, ax.tickvec.b=7)  
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export
PLS.GLM.biplot_SIMPLS = function(X, y, algorithm=NULL, ax.tickvec.X=NULL, ax.tickvec.y=NULL, ax.tickvec.b=NULL, ... ) 
{    
  options(digits=3)
  D = cbind(X,y)
  N = nrow(D)
  X.mean = apply(X, 2, mean)  #column means of X or colMeans(X)
  y.mean = mean(y)  #mean of y 
  X.std = apply(X, 2, sd)  #column standard deviations of X 
  y.std = sd(y)  #standard deviation of y   
  X.scal = scale(X, center=TRUE, scale=TRUE)  #scaled X
  y.scal = scale(y, center=TRUE, scale=TRUE)  #scaled y  
  
  #biplot 
  A = 2
  main = algorithm(X.scal,y.scal,A,PLS_algorithm=mod.SIMPLS,eps=1e-04)
  Tmat = main$X.scores
  Pmat = main$X.loadings
  qvec = main$y.loadings
  Rmat = main$X.weights.trans	
  bvec = Rmat %*% t(qvec)  #(Px1) estimated PLS-GLM coefficient vector 	 
  dimnames(bvec) = list(colnames(X), colnames(y))	
  
  #approximation		 
  X.hat = ( (Tmat %*% t(Pmat)) * rep(X.std,each=N) ) + rep(X.mean,each=N)
  dimnames(X.hat) = dimnames(X)
  y.hat = ( (Tmat %*% t(qvec)) * rep(y.std,each=N) ) + rep(y.mean,each=N)
  dimnames(y.hat) = list(rownames(y), paste("Expected_",colnames(y),sep=""))
  D.hat = cbind(X.hat, y.hat)	  
  
  #biplot points
  Z = rbind(Tmat, Rmat)  #((N+P)x2) biplot points
  dimnames(Z) = list(c(rownames(X),paste("b",1:ncol(X),sep="")), NULL)
  sample.style = biplot.sample.control(2, label=TRUE)
  sample.style$col = c("red", "purple")
  sample.style$pch = c(15,16)
  Gmat = cbind(c(rep(1,nrow(X)),rep(0,ncol(X))), c(rep(0,nrow(X)),rep(1,ncol(X))))
  dimnames(Gmat) = list(NULL, c("samples","coefficients"))
  classes = 1:2	 
  
  #axes direction    
  #X.hat=TP'
  X.axes.direction = (1 / (diag(Pmat%*%t(Pmat)))) * Pmat       	   
  #y.hat=Tq
  qvec = matrix(qvec, nrow=1)  #since for M=1, j = 1:1
  y.axes.direction = (1 / (diag(qvec%*%t(qvec)))) * qvec
  #b.pls-glm=RQ'
  b.axes.direction = (1 / (diag(qvec%*%t(qvec)))) * qvec	
  
  #calibration of axes	 	
  z.axes.X = lapply(1:ncol(X), calibrate.axis, X, X.mean, X.std, X.axes.direction, 
                    1:ncol(X), ax.tickvec.X, rep(0,ncol(X)), rep(0,ncol(X)), NULL)  
  z.axes.y = lapply(1, calibrate.axis, y, y.mean, y.std, y.axes.direction, 
                    1, ax.tickvec.y, rep(0,1), rep(0,1), NULL)
  z.axes.b = lapply(1, calibrate.axis, y.scal, rep(0,1), rep(1,1), b.axes.direction, 
                    1, ax.tickvec.b, rep(0,1), rep(0,1), NULL)		
  z.axes = vector("list", ncol(X)+ncol(y)+ncol(bvec))
  for (i in 1:ncol(X)) 
  {
    z.axes[[i]] = z.axes.X[[i]]
  }
  for (i in 1:ncol(y)) 
  {
    z.axes[[i+ncol(X)]] = z.axes.y[[i]]
  }
  for (i in 1:ncol(bvec)) 
  {
    z.axes[[i+ncol(X)+ncol(y)]] = z.axes.b[[i]]  
  }			  
  
  #biplot axes 
  ax.style = biplot.ax.control(ncol(X)+ncol(y)+ncol(bvec), c(colnames(X), colnames(y), colnames(bvec)))
  #for X	 
  ax.style$col[1:ncol(X)] = "blue"
  ax.style$label.col[1:ncol(X)] = "blue"
  ax.style$tick.col[1:ncol(X)] = "blue"
  ax.style$tick.label.col[1:ncol(X)] = "blue"  
  #for y	 
  ax.style$col[ncol(X)+1:ncol(y)] = "black"
  ax.style$label.col[ncol(X)+1:ncol(y)] = "black"
  ax.style$tick.col[ncol(X)+1:ncol(y)] = "black"
  ax.style$tick.label.col[ncol(X)+1:ncol(y)] = "black"  
  #for b.pls-glm
  ax.style$col[ncol(X)+ncol(y)+1:ncol(bvec)] = "black"
  ax.style$label.col[ncol(X)+ncol(y)+1:ncol(bvec)] = "black" 
  ax.style$tick.col[ncol(X)+ncol(y)+1:ncol(bvec)] = "purple"
  ax.style$tick.label.col[ncol(X)+ncol(y)+1:ncol(bvec)] = "purple"   
  draw.biplot(Z, G=Gmat, classes=classes, sample.style=sample.style, z.axes=z.axes, ax.style=ax.style, VV=rbind(X.axes.direction,y.axes.direction)) 
  
  list(D.hat=D.hat, bvec=bvec)  	 
}

#' The Partial Least Squares (PLS) biplot for Generalized Linear Model (GLM),  with the labels of the samples, coefficient points and tick markers excluded 
#' 
#' Takes in a set of predictor variables and a set of response variables and produces a PLS biplot for the (univariate) GLMs, with the labels of the samples, coefficient points and tick markers excluded.
#' @param X A (NxP) predictor matrix 
#' @param y A (Nx1) response vector 
#' @param algorithm The PLS.GLM algorithm 
#' @param ax.tickvec.X tick marker length for each X-variable axis in the biplot
#' @param ax.tickvec.y tick marker length for the y-variable axis in the biplot
#' @param ax.tickvec.b (purple) tick marker length for the y-variable axis in the biplot
#' @param ... Other arguments. Currently ignored
#' @return The PLS biplot of a GLM of D=[X y] with some parameters
#' @examples
#' if(require(robustbase))
#' possum.mat 
#' y = as.matrix(possum.mat[,1], ncol=1) 
#' dimnames(y) = list(paste("S", 1:nrow(possum.mat), seq=""), "Diversity")
#' X = as.matrix(possum.mat[,2:14], ncol=13) 
#' dimnames(X) = list(paste("S", 1:nrow(possum.mat), seq=""), colnames(possum.mat[,2:14]))  
#' PLS.GLM.biplot_no_labels(X, y, algorithm=PLS.GLM, ax.tickvec.X=rep(5,ncol(X)), 
#' ax.tickvec.y=10, ax.tickvec.b=7)  
#'
#' #Pima.tr data
#' if(require(MASS))
#' data(Pima.tr, package="MASS")
#' X = as.matrix(cbind(Pima.tr[,1:7]))  
#' dimnames(X) = list(1:nrow(X), colnames(X))
#' y = as.matrix(as.numeric(Pima.tr$type)-1, ncol=1)
#' #0=No and 1=Yes
#' dimnames(y) = list(1:nrow(y), paste("type"))
#' PLS.GLM.biplot_no_labels(X, y, algorithm=PLS.binomial.GLM, 
#' ax.tickvec.X=c(3,3,8,7,8,5,2), ax.tickvec.y=3, ax.tickvec.b=3)  
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export
PLS.GLM.biplot_no_labels = function(X, y, algorithm=NULL, ax.tickvec.X=NULL, ax.tickvec.y=NULL, ax.tickvec.b=NULL, ... ) 
{    
  options(digits=3)
  D = cbind(X,y)
  N = nrow(D)
  X.mean = apply(X, 2, mean)  #column means of X or colMeans(X)
  y.mean = mean(y)  #mean of y 
  X.std = apply(X, 2, sd)  #column standard deviations of X 
  y.std = sd(y)  #standard deviation of y   
  X.scal = scale(X, center=TRUE, scale=TRUE)  #scaled X
  y.scal = scale(y, center=TRUE, scale=TRUE)  #scaled y  
  
  #biplot 
  A = 2
  main = algorithm(X.scal,y.scal,A,PLS_algorithm=mod.NIPALS,eps=1e-04)
  Tmat = main$X.scores
  Pmat = main$X.loadings
  qvec = main$y.loadings
  Rmat = main$X.weights.trans	
  bvec = Rmat %*% t(qvec)  #(Px1) estimated PLS-GLM coefficient vector 	 
  dimnames(bvec) = list(colnames(X), colnames(y))	
  
  #approximation		 
  X.hat = ( (Tmat %*% t(Pmat)) * rep(X.std,each=N) ) + rep(X.mean,each=N)
  dimnames(X.hat) = dimnames(X)
  y.hat = ( (Tmat %*% t(qvec)) * rep(y.std,each=N) ) + rep(y.mean,each=N)
  dimnames(y.hat) = list(rownames(y), paste("Expected_",colnames(y),sep=""))
  D.hat = cbind(X.hat, y.hat)	  
  
  #biplot points
  Z = rbind(Tmat, Rmat)  #((N+P)x2) biplot points
  dimnames(Z) = list(c(rownames(X),paste("b",1:ncol(X),sep="")), NULL)
  sample.style = biplot.sample.control(2, label=FALSE)
  sample.style$col = c("red", "purple")
  sample.style$pch = c(15,16)
  Gmat = cbind(c(rep(1,nrow(X)),rep(0,ncol(X))), c(rep(0,nrow(X)),rep(1,ncol(X))))
  dimnames(Gmat) = list(NULL, c("samples","coefficients"))
  classes = 1:2	 
  
  #axes direction    
  #X.hat=TP'
  X.axes.direction = (1 / (diag(Pmat%*%t(Pmat)))) * Pmat       	   
  #y.hat=Tq
  qvec = matrix(qvec, nrow=1)  #since for M=1, j = 1:1
  y.axes.direction = (1 / (diag(qvec%*%t(qvec)))) * qvec
  #b.pls-glm=RQ'
  b.axes.direction = (1 / (diag(qvec%*%t(qvec)))) * qvec	
  
  #calibration of axes	 	
  z.axes.X = lapply(1:ncol(X), calibrate.axis, X, X.mean, X.std, X.axes.direction, 
                    1:ncol(X), ax.tickvec.X, rep(0,ncol(X)), rep(0,ncol(X)), NULL)  
  z.axes.y = lapply(1, calibrate.axis, y, y.mean, y.std, y.axes.direction, 
                    1, ax.tickvec.y, rep(0,1), rep(0,1), NULL)
  z.axes.b = lapply(1, calibrate.axis, y.scal, rep(0,1), rep(1,1), b.axes.direction, 
                    1, ax.tickvec.b, rep(0,1), rep(0,1), NULL)		
  z.axes = vector("list", ncol(X)+ncol(y)+ncol(bvec))
  for (i in 1:ncol(X)) 
  {
    z.axes[[i]] = z.axes.X[[i]]
  }
  for (i in 1:ncol(y)) 
  {
    z.axes[[i+ncol(X)]] = z.axes.y[[i]]
  }
  for (i in 1:ncol(bvec)) 
  {
    z.axes[[i+ncol(X)+ncol(y)]] = z.axes.b[[i]]  
  }			  
  
  #biplot axes 
  ax.style = biplot.ax.control(ncol(X)+ncol(y)+ncol(bvec), c(colnames(X), colnames(y), colnames(bvec)))
  ax.style$tick.col[1:ncol(X)] = "azure"   
  ax.style$tick.label.col[1:ncol(X)] = "azure"  
  ax.style$tick.col[ncol(X)+1:ncol(y)] = "azure"  
  ax.style$tick.label.col[ncol(X)+1:ncol(y)] = "azure" 
  ax.style$tick.col[ncol(X)+ncol(y)+1:ncol(bvec)] = "azure"  
  ax.style$tick.label.col[ncol(X)+ncol(y)+1:ncol(bvec)] = "azure"
  #for X	 
  ax.style$col[1:ncol(X)] = "blue"
  ax.style$label.col[1:ncol(X)] = "blue" 
  #for y	 
  ax.style$col[ncol(X)+1:ncol(y)] = "black"
  ax.style$label.col[ncol(X)+1:ncol(y)] = "black"  
  #for b.pls-glm
  ax.style$col[ncol(X)+ncol(y)+1:ncol(bvec)] = "black"
  ax.style$label.col[ncol(X)+ncol(y)+1:ncol(bvec)] = "black"    
  draw.biplot(Z, G=Gmat, classes=classes, sample.style=sample.style, z.axes=z.axes, ax.style=ax.style, VV=rbind(X.axes.direction,y.axes.direction)) 
  
  list(D.hat=D.hat, bvec=bvec)  	 
}      

#' The Partial Least Squares (PLS) biplot for Generalized Linear Model (GLM) with the labels of the sample points excluded 
#' 
#' Takes in a set of predictor variables and a set of response variables and produces a PLS biplot for the (univariate) GLMs, with the labels of the sample points excluded.
#' @param X A (NxP) predictor matrix 
#' @param y A (Nx1) response vector 
#' @param algorithm The PLS.GLM algorithm 
#' @param ax.tickvec.X tick marker length for each X-variable axis in the biplot
#' @param ax.tickvec.y tick marker length for the y-variable axis in the biplot
#' @param ax.tickvec.b (purple) tick marker length for the y-variable axis in the biplot
#' @param ... Other arguments. Currently ignored
#' @return The PLS biplot of a GLM of D=[X y] with some parameters
#' @examples
#' if(require(robustbase))
#' possum.mat 
#' y = as.matrix(possum.mat[,1], ncol=1) 
#' dimnames(y) = list(paste("S", 1:nrow(possum.mat), seq=""), "Diversity")
#' X = as.matrix(possum.mat[,2:14], ncol=13) 
#' dimnames(X) = list(paste("S", 1:nrow(possum.mat), seq=""), colnames(possum.mat[,2:14]))  
#' #Poisson-fitted
#' PLS.GLM.biplot_no.SN(X, y, algorithm=PLS.GLM, ax.tickvec.X=rep(5,ncol(X)),
#'  ax.tickvec.y=10, ax.tickvec.b=7)  
#'
#' #Pima.tr data
#' if(require(MASS))
#' data(Pima.tr, package="MASS")
#' X = as.matrix(cbind(Pima.tr[,1:7]))  
#' dimnames(X) = list(1:nrow(X), colnames(X))
#' y = as.matrix(as.numeric(Pima.tr$type)-1, ncol=1)
#' #0=No and 1=Yes
#' dimnames(y) = list(1:nrow(y), paste("type"))
#' PLS.GLM.biplot_no.SN(X, y, algorithm=PLS.binomial.GLM, 
#' ax.tickvec.X=c(3,3,8,7,8,5,2), ax.tickvec.y=3, ax.tickvec.b=3)    
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export
PLS.GLM.biplot_no.SN = function(X, y, algorithm=NULL, ax.tickvec.X=NULL, ax.tickvec.y=NULL, ax.tickvec.b=NULL, ... ) 
{    
  options(digits=3)
  D = cbind(X,y)
  N = nrow(D)
  X.mean = apply(X, 2, mean)  #column means of X or colMeans(X)
  y.mean = mean(y)  #mean of y 
  X.std = apply(X, 2, sd)  #column standard deviations of X 
  y.std = sd(y)  #standard deviation of y   
  X.scal = scale(X, center=TRUE, scale=TRUE)  #scaled X
  y.scal = scale(y, center=TRUE, scale=TRUE)  #scaled y
  
  #biplot 
  A = 2
  main = algorithm(X.scal,y.scal,A,PLS_algorithm=mod.NIPALS,eps=1e-04)
  Tmat = main$X.scores
  Pmat = main$X.loadings
  qvec = main$y.loadings
  Rmat = main$X.weights.trans
  bvec = Rmat %*% t(qvec)  #(Px1) estimated PLS-GLM coefficient vector  
  dimnames(bvec) = list(colnames(X), colnames(y))
  
  #approximation 
  X.hat = ( (Tmat %*% t(Pmat)) * rep(X.std,each=N) ) + rep(X.mean,each=N)
  dimnames(X.hat) = dimnames(X)
  y.hat = ( (Tmat %*% t(qvec)) * rep(y.std,each=N) ) + rep(y.mean,each=N)
  dimnames(y.hat) = list(rownames(y), paste("Expected_",colnames(y),sep=""))
  D.hat = cbind(X.hat, y.hat)  
  
  #biplot points
  Z = rbind(Tmat, Rmat)  #((N+P)x2) biplot points
  dimnames(Z) = list(c(rownames(X),paste("b",1:ncol(X),sep="")), NULL)
  sample.style = biplot.sample.control(2, label=c(FALSE,TRUE))
  sample.style$col = c("red", "purple")
  sample.style$pch = c(15,16)
  Gmat = cbind(c(rep(1,nrow(X)),rep(0,ncol(X))), c(rep(0,nrow(X)),rep(1,ncol(X))))
  dimnames(Gmat) = list(NULL, c("samples","coefficients"))
  classes = 1:2 
  
  #axes direction    
  #X.hat=TP'
  X.axes.direction = (1 / (diag(Pmat%*%t(Pmat)))) * Pmat          
  #y.hat=Tq
  qvec = matrix(qvec, nrow=1)  #since for M=1, j = 1:1
  y.axes.direction = (1 / (diag(qvec%*%t(qvec)))) * qvec
  #b.pls-glm=RQ'
  b.axes.direction = (1 / (diag(qvec%*%t(qvec)))) * qvec
  
  #calibration of axes 
  z.axes.X = lapply(1:ncol(X), calibrate.axis, X, X.mean, X.std, X.axes.direction, 
                    1:ncol(X), ax.tickvec.X, rep(0,ncol(X)), rep(0,ncol(X)), NULL)  
  z.axes.y = lapply(1, calibrate.axis, y, y.mean, y.std, y.axes.direction, 
                    1, ax.tickvec.y, rep(0,1), rep(0,1), NULL)
  z.axes.b = lapply(1, calibrate.axis, y.scal, rep(0,1), rep(1,1), b.axes.direction, 
                    1, ax.tickvec.b, rep(0,1), rep(0,1), NULL)
  z.axes = vector("list", ncol(X)+ncol(y)+ncol(bvec))
  for (i in 1:ncol(X)) 
  {
    z.axes[[i]] = z.axes.X[[i]]
  }
  for (i in 1:ncol(y)) 
  {
    z.axes[[i+ncol(X)]] = z.axes.y[[i]]
  }
  for (i in 1:ncol(bvec)) 
  {
    z.axes[[i+ncol(X)+ncol(y)]] = z.axes.b[[i]]  
  }  
  
  #biplot axes 
  ax.style = biplot.ax.control(ncol(X)+ncol(y)+ncol(bvec), c(colnames(X), colnames(y), colnames(bvec)))
  #for X 
  ax.style$col[1:ncol(X)] = "blue"
  ax.style$label.col[1:ncol(X)] = "blue"
  ax.style$tick.col[1:ncol(X)] = "blue"
  ax.style$tick.label.col[1:ncol(X)] = "blue"  
  #for y 
  ax.style$col[ncol(X)+1:ncol(y)] = "black"
  ax.style$label.col[ncol(X)+1:ncol(y)] = "black"
  ax.style$tick.col[ncol(X)+1:ncol(y)] = "black"
  ax.style$tick.label.col[ncol(X)+1:ncol(y)] = "black"  
  #for b.pls-glm
  ax.style$col[ncol(X)+ncol(y)+1:ncol(bvec)] = "black"
  ax.style$label.col[ncol(X)+ncol(y)+1:ncol(bvec)] = "black" 
  ax.style$tick.col[ncol(X)+ncol(y)+1:ncol(bvec)] = "purple"
  ax.style$tick.label.col[ncol(X)+ncol(y)+1:ncol(bvec)] = "purple"   
  draw.biplot(Z, G=Gmat, classes=classes, sample.style=sample.style, z.axes=z.axes, ax.style=ax.style, VV=rbind(X.axes.direction,y.axes.direction)) 
  
  list(D.hat=D.hat, bvec=bvec)   
}

#' The Partial Least Squares (PLS) biplot for Generalized Linear Model (GLM) fitted using the SIMPLS algorithm, with the labels of the sample points excluded 
#' 
#' Takes in a set of predictor variables and a set of response variables and produces a PLS biplot for the (univariate) GLMs, with the labels of the sample points excluded.
#' @param X A (NxP) predictor matrix 
#' @param y A (Nx1) response vector 
#' @param algorithm The PLS.GLM_SIMPLS algorithm
#' @param ax.tickvec.X tick marker length for each X-variable axis in the biplot
#' @param ax.tickvec.y tick marker length for the y-variable axis in the biplot
#' @param ax.tickvec.b (purple) tick marker length for the y-variable axis in the biplot
#' @param ... Other arguments. Currently ignored
#' @return The PLS biplot of a GLM (fitted using the SIMPLS algorithm) of D=[X y] with some parameters
#' @examples
#' if(require(robustbase))
#' possum.mat 
#' y = as.matrix(possum.mat[,1], ncol=1) 
#' dimnames(y) = list(paste("S", 1:nrow(possum.mat), seq=""), "Diversity")
#' X = as.matrix(possum.mat[,2:14], ncol=13) 
#' dimnames(X) = list(paste("S", 1:nrow(possum.mat), seq=""), colnames(possum.mat[,2:14]))  
#' #Poisson-fitted
#' PLS.GLM.biplot_SIMPLS_no.SN(X, y, algorithm=PLS.GLM, 
#' ax.tickvec.X=rep(5,ncol(X)), ax.tickvec.y=10, ax.tickvec.b=7)  
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export
PLS.GLM.biplot_SIMPLS_no.SN = function(X, y, algorithm=NULL, ax.tickvec.X=NULL, ax.tickvec.y=NULL, ax.tickvec.b=NULL, ... ) 
{    
  options(digits=3)
  D = cbind(X,y)
  N = nrow(D)
  X.mean = apply(X, 2, mean)  #column means of X or colMeans(X)
  y.mean = mean(y)  #mean of y 
  X.std = apply(X, 2, sd)  #column standard deviations of X 
  y.std = sd(y)  #standard deviation of y   
  X.scal = scale(X, center=TRUE, scale=TRUE)  #scaled X
  y.scal = scale(y, center=TRUE, scale=TRUE)  #scaled y  
  
  #biplot 
  A = 2
  main = algorithm(X.scal,y.scal,A,PLS_algorithm=mod.SIMPLS,eps=1e-04)
  Tmat = main$X.scores
  Pmat = main$X.loadings
  qvec = main$y.loadings
  Rmat = main$X.weights.trans	
  bvec = Rmat %*% t(qvec)  #(Px1) estimated PLS-GLM coefficient vector 	 
  dimnames(bvec) = list(colnames(X), colnames(y))	
  
  #approximation		 
  X.hat = ( (Tmat %*% t(Pmat)) * rep(X.std,each=N) ) + rep(X.mean,each=N)
  dimnames(X.hat) = dimnames(X)
  y.hat = ( (Tmat %*% t(qvec)) * rep(y.std,each=N) ) + rep(y.mean,each=N)
  dimnames(y.hat) = list(rownames(y), paste("Expected_",colnames(y),sep=""))
  D.hat = cbind(X.hat, y.hat)	  
  
  #biplot points
  Z = rbind(Tmat, Rmat)  #((N+P)x2) biplot points
  dimnames(Z) = list(c(rownames(X),paste("b",1:ncol(X),sep="")), NULL)
  sample.style = biplot.sample.control(2, label=c(FALSE,TRUE))
  sample.style$col = c("red", "purple")
  sample.style$pch = c(15,16)
  Gmat = cbind(c(rep(1,nrow(X)),rep(0,ncol(X))), c(rep(0,nrow(X)),rep(1,ncol(X))))
  dimnames(Gmat) = list(NULL, c("samples","coefficients"))
  classes = 1:2	 
  
  #axes direction    
  #X.hat=TP'
  X.axes.direction = (1 / (diag(Pmat%*%t(Pmat)))) * Pmat       	   
  #y.hat=Tq
  qvec = matrix(qvec, nrow=1)  #since for M=1, j = 1:1
  y.axes.direction = (1 / (diag(qvec%*%t(qvec)))) * qvec
  #b.pls-glm=RQ'
  b.axes.direction = (1 / (diag(qvec%*%t(qvec)))) * qvec	
  
  #calibration of axes	 	
  z.axes.X = lapply(1:ncol(X), calibrate.axis, X, X.mean, X.std, X.axes.direction, 
                    1:ncol(X), ax.tickvec.X, rep(0,ncol(X)), rep(0,ncol(X)), NULL)  
  z.axes.y = lapply(1, calibrate.axis, y, y.mean, y.std, y.axes.direction, 
                    1, ax.tickvec.y, rep(0,1), rep(0,1), NULL)
  z.axes.b = lapply(1, calibrate.axis, y.scal, rep(0,1), rep(1,1), b.axes.direction, 
                    1, ax.tickvec.b, rep(0,1), rep(0,1), NULL)		
  z.axes = vector("list", ncol(X)+ncol(y)+ncol(bvec))
  for (i in 1:ncol(X)) 
  {
    z.axes[[i]] = z.axes.X[[i]]
  }
  for (i in 1:ncol(y)) 
  {
    z.axes[[i+ncol(X)]] = z.axes.y[[i]]
  }
  for (i in 1:ncol(bvec)) 
  {
    z.axes[[i+ncol(X)+ncol(y)]] = z.axes.b[[i]]  
  }			  
  
  #biplot axes 
  ax.style = biplot.ax.control(ncol(X)+ncol(y)+ncol(bvec), c(colnames(X), colnames(y), colnames(bvec)))
  #for X	 
  ax.style$col[1:ncol(X)] = "blue"
  ax.style$label.col[1:ncol(X)] = "blue"
  ax.style$tick.col[1:ncol(X)] = "blue"
  ax.style$tick.label.col[1:ncol(X)] = "blue"  
  #for y	 
  ax.style$col[ncol(X)+1:ncol(y)] = "black"
  ax.style$label.col[ncol(X)+1:ncol(y)] = "black"
  ax.style$tick.col[ncol(X)+1:ncol(y)] = "black"
  ax.style$tick.label.col[ncol(X)+1:ncol(y)] = "black"  
  #for b.pls-glm
  ax.style$col[ncol(X)+ncol(y)+1:ncol(bvec)] = "black"
  ax.style$label.col[ncol(X)+ncol(y)+1:ncol(bvec)] = "black" 
  ax.style$tick.col[ncol(X)+ncol(y)+1:ncol(bvec)] = "purple"
  ax.style$tick.label.col[ncol(X)+ncol(y)+1:ncol(bvec)] = "purple"   
  draw.biplot(Z, G=Gmat, classes=classes, sample.style=sample.style, z.axes=z.axes, ax.style=ax.style, VV=rbind(X.axes.direction,y.axes.direction)) 
  
  list(D.hat=D.hat, bvec=bvec)  	 
} 

#' A zoomed-in display of the coefficient points in the Partial Least Squares (PLS) biplot for Generalized Linear Model (GLM)   
#' 
#' Takes in a set of predictor variables and a set of response variables and produces a zoomed-in display of the coefficient points in the PLS biplot for the (univariate) GLMs.
#' @param X A (NxP) predictor matrix 
#' @param y A (Nx1) response vector 
#' @param algorithm The PLS.GLM_SIMPLS algorithm
#' @param ax.tickvec.b (purple) tick marker length for the y-variable axis in the biplot
#' @param ... Other arguments. Currently ignored
#' @return A zoomed-in display of the coefficient points in the PLS biplot of a GLM of D=[X y] with some parameters
#' @examples  
#' if(require(robustbase))
#' possum.mat 
#' y = as.matrix(possum.mat[,1], ncol=1) 
#' dimnames(y) = list(paste("S", 1:nrow(possum.mat), seq=""), "Diversity")
#' X = as.matrix(possum.mat[,2:14], ncol=13) 
#' dimnames(X) = list(paste("S", 1:nrow(possum.mat), seq=""), colnames(possum.mat[,2:14]))
#' PLS.GLM.biplot_bvec(X, y, algorithm=PLS.GLM, ax.tickvec.b=10)
#' 
#' #Pima.tr data
#' if(require(MASS))
#' data(Pima.tr, package="MASS")
#' X = as.matrix(cbind(Pima.tr[,1:7]))  
#' dimnames(X) = list(1:nrow(X), colnames(X))
#' y = as.matrix(as.numeric(Pima.tr$type)-1, ncol=1)
#' #0=No and 1=Yes
#' dimnames(y) = list(1:nrow(y), paste("type"))
#' PLS.GLM.biplot_bvec(X, y, algorithm=PLS.binomial.GLM,ax.tickvec.b=10)   
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export
PLS.GLM.biplot_bvec = function(X, y, algorithm=NULL, ax.tickvec.b=NULL, ... ) 
{    
  options(digits=3)
  X.scal = scale(X, center=TRUE, scale=TRUE)  #scaled X
  y.scal = scale(y, center=TRUE, scale=TRUE)  #scaled y  
  
  #biplot 
  A = 2
  main = algorithm(X.scal,y.scal,A,PLS_algorithm=mod.NIPALS,eps=0.001)
  qvec = main$y.loadings
  Rmat = main$X.weights.trans	
  bvec = Rmat %*% t(qvec)  #(Px1) estimated PLS-GLM coefficient vector 	 
  dimnames(bvec) = list(colnames(X), colnames(y))	    
  
  #biplot points
  Z = Rmat  #
  dimnames(Z) = list(paste("b",1:ncol(X),sep=""), NULL)
  sample.style = biplot.sample.control(1, label=TRUE)
  sample.style$col = "purple"
  sample.style$pch = 16
  Gmat = cbind(rep(1,ncol(X)))
  dimnames(Gmat) = list(NULL, c("coefficients"))
  classes = 1	 
  
  #axes direction   
  qvec = matrix(qvec, nrow=1)  #since for M=1, j = 1:1
  #b.pls-glm=RQ'
  b.axes.direction = (1 / (diag(qvec%*%t(qvec)))) * qvec	
  
  #calibration of axes	 	
  z.axes.b = lapply(1, calibrate.axis, y.scal, rep(0,1), rep(1,1), b.axes.direction, 
                    1, ax.tickvec.b, rep(0,1), rep(0,1), NULL)		
  z.axes = vector("list", ncol(bvec))      
  for (i in 1:ncol(bvec)) 
  {
    z.axes[[i]] = z.axes.b[[i]]  
  }			  
  
  #biplot axes 
  ax.style = biplot.ax.control(ncol(bvec), c(colnames(bvec)))     
  #for b.pls-glm
  ax.style$col[1:ncol(bvec)] = "black"
  ax.style$label.col[1:ncol(bvec)] = "black" 
  ax.style$tick.col[1:ncol(bvec)] = "purple"
  ax.style$tick.label.col[1:ncol(bvec)] = "purple"   
  draw.biplot(Z, G=Gmat, classes=classes, sample.style=sample.style, z.axes=z.axes, ax.style=ax.style, VV=b.axes.direction) 
  
  list(bvec=bvec)  	 
}

# --------------------------------------------------------------------------------------------------------------------



#A.7  Chapter 7:  Biplots for Sparse Partial Least Squares.
# Sparse Partial Least Squares (SPLS) and SPLS-GLMs algorithms.
# PLS biplot for the SPLS and SPLS-GLMs. 
# RGui(32-bit).

#' Sparse Partial Least Squares (SPLS) algorithm
#' 
#' Takes in a set of predictor variables and a set of response variables and gives the SPLS parameters.
#' @param X A (NxP) predictor matrix 
#' @param Y A (NxM) response matrix
#' @param A The number of PLS components 
#' @param lambdaY A value for the penalty parameters for the soft-thresholding penalization function for Y-weights 
#' @param lambdaX A value for the penalty parameters for the soft-thresholding penalization function for X-weights
#' @param eps Cut off value for convergence step
#' @param ... Other arguments. Currently ignored
#' @return The SPLS parameters of D=[X Y]
#' @examples 
#' if(require(chemometrics))
#' data(ash, package="chemometrics")
#' X1 = as.matrix(ash[,10:17], ncol=8)  
#' Y1 = as.matrix(ash$SOT)   
#' colnames(Y1) = paste("SOT")
#' mod.SPLS(X=scale(X1), Y=scale(Y1), A=2, lambdaY=0, lambdaX=10.10, eps=1e-5)
#' #lambdaX and lambdaY value are determined using function opt.penalty.values 
#' #for more details, see opt.penalty.values help file
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export
mod.SPLS = function(X, Y, A, lambdaY, lambdaX, eps, ...)
{
  #lambdaY & lambdaX are the value for the penalty parameters for the soft-thresholding penalization functions
  X.cen = scale(X, center=TRUE, scale=FALSE)  #centered X
  Y.cen = scale(Y, center=TRUE, scale=FALSE)  #centered Y
  Y.mean = colMeans(Y)  #(1xM) column means of Y
  P = ncol(X)
  M = ncol(Y) 
  N = nrow(X)
  R.ash = P.ash = matrix(0, nrow=P, ncol=A)  #(PxA) X.weights and X.loadings                                                                                                                                                                            #matrices             
  C.ash = Q.ash = matrix(0, nrow=M, ncol=A)  #(MxA) Y.weights and Y.loadings
  T.ash = matrix(0, nrow=N, ncol=A)   #(NxA) X.scores matrix
  
  #soft thresholding penalization
  g.vec = function(u, lambda=NULL, ...)
  {
    sign(u) * ifelse ((abs(u)-lambda>0),(abs(u)-lambda), 0)
  } 
  #steps 2 to 9
  S = t(X.cen) %*% Y.cen             
  for (a in 1:A)
  { 
    r.a = svd(S)$u[,1]  #(Px1) X.weight
    c.a = svd(S)$v[,1]  #(Mx1) Y.weight     
    r.a.old = r.a
    c.a.old = c.a      
    conv.r.a = conv.c.a = 99
    while (conv.r.a > eps & conv.c.a > eps)
    {
      r.a.new = g.vec(u=S%*%c.a.old, lambda=lambdaX) 
      if (!all(r.a.new==0)) r.a.new = r.a.new / as.numeric(sqrt(t(r.a.new)%*%r.a.new))
      c.a.new = g.vec(u=t(S)%*%r.a.old, lambda=lambdaY)  
      if (!all(c.a.new==0)) c.a.new = c.a.new / as.numeric(sqrt(t(c.a.new)%*%c.a.new))
      if (all(r.a.new==0)) conv.r.a = 0 else { old.non.zero = (c(r.a.old)!=0)
                                               r.a.old1 = r.a.old[old.non.zero]
                                               r.a.new1 = r.a.new[old.non.zero]
                                               if (length(r.a.old1)>0)
                                                 conv.r.a = sum(abs(r.a.new1 - r.a.old1)/abs(r.a.old1))
                                               else conv.r.a = 0    
      }
      if (all(c.a.new==0)) conv.c.a = 0 else { old.non.zero = (c(c.a.old)!=0)
                                               c.a.old1 = c.a.old[old.non.zero]
                                               c.a.new1 = c.a.new[old.non.zero]
                                               if (length(c.a.old1)>0)
                                                 conv.c.a = sum(abs(c.a.new1 - c.a.old1)/abs(c.a.old1))
                                               else conv.c.a = 0
      }
      r.a = r.a.new                                                                              
      r.a.old = r.a.new      
      c.a = c.a.new                                                                             
      c.a.old = c.a.new   
    }   
    t.a = X.cen %*% r.a  #(Nx1) X.score
    if (!all(t.a==0))  t.a = t.a / as.numeric(sqrt(t(X.cen%*%r.a)%*%X.cen%*%r.a))                
    p.a = t(X.cen) %*% t.a  #(Px1) X.loading
    q.a = t(Y.cen) %*% t.a  #(Mx1) Y.loading
    
    #update S
    if (!all(p.a==0)) S = S - p.a%*%solve(t(p.a)%*%p.a)%*%t(p.a)%*%S                  
    
    R.ash[, a] = r.a   #(PxA) X.weights matrix      
    C.ash[, a] = c.a   #(MxA) Y.weights matrix
    P.ash[, a] = p.a   #(PxA) X.loadings matrix
    Q.ash[, a] = q.a   #(MxA) Y.loadings matrix  
    T.ash[, a] = t.a   #(NxA) X.scores matrix  
  }                          
  dimnames(R.ash) = dimnames(P.ash) = list(colnames(X), paste("Comp", 1:A, seq=""))
  dimnames(Q.ash) = dimnames(C.ash) = list (colnames(Y), paste("Comp", 1:A, seq=""))               
  dimnames(T.ash) = list(rownames(X), paste("Comp", 1:A, seq=""))
  
  B.hat = R.ash %*% t(Q.ash)
  if (all(T.ash[,1]==0) | all(T.ash[,2]==0)) {
    Y.hat = (T.ash%*%t(Q.ash)) + rep(Y.mean,each=N) #predicted Y-values
  }
  else {
    if (!all(T.ash[,1]==0) | !all(T.ash[,2]==0)) {
      Y.hat = (T.ash%*%solve(t(T.ash)%*%T.ash)%*%t(Q.ash)) + rep(Y.mean,each=N) #predicted Y-values 
    }   
  }
  rmsep = sqrt(sum((Y-Y.hat)^2)/N)  #Root Mean Squared Error of Prediction (RMSEP)  
  
  #selecting non-zero values
  X.select = which(!R.ash[,1]==0, arr.ind=TRUE)
  Y.select = which(!C.ash[,1]==0, arr.ind=TRUE)
  
  list(X.scores=T.ash, X.weights.trans=R.ash, Y.weights=C.ash, X.loadings=P.ash, Y.loadings=Q.ash, B.mat=B.hat, Y.hat=Y.hat, RMSEP=rmsep, 
       X.select=X.select, Y.select=Y.select) 
} 

#' Sparse Partial Least Squares-Generalized Linear Model (SPLS-GLM) algorithm
#' 
#' Takes in a set of predictor variables and a set of response variables and gives the SPLS-GLM parameters.
#' @param X A (NxP) predictor matrix 
#' @param y A (Nx1) Poisson-distributed response vector 
#' @param A The number of PLS components 
#' @param lambdaY A value for the penalty parameters for the soft-thresholding penalization function for Y-weights 
#' @param lambdaX A value for the penalty parameters for the soft-thresholding penalization function for X-weights 
#' @param eps Cut off value for convergence step
#' @param ... Other arguments. Currently ignored
#' @return The SPLS-GLM parameters of D=[X y]
#' @examples  
#' if(require(robustbase))
#' possum.mat 
#' y = as.matrix(possum.mat[,1], ncol=1) 
#' dimnames(y) = list(paste("S", 1:nrow(possum.mat), seq=""), "Diversity")
#' X = as.matrix(possum.mat[,2:14], ncol=13) 
#' dimnames(X) = list(paste("S", 1:nrow(possum.mat), seq=""), colnames(possum.mat[,2:14])) 
#' SPLS.GLM(scale(X), scale(y), A=2, lambdaY=0, lambdaX=3.3, eps=1e-3)
#' #lambdaX and lambdaY value are determined using function opt.penalty.values 
#' #for more details, see opt.penalty.values help file
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export 
SPLS.GLM = function(X, y, A, lambdaY, lambdaX, eps=0.001, ...)
{
  #For M=1
  #y is Poisson distributed
  #A=initial number of components  
  #lambdaY & lambdaX are the value for the penalty parameters for the soft-thresholding penalization functions
  #1. Initialization
  N = nrow(X)
  P = ncol(X)  
  M = 1 #ncol(y)
  eta = matrix(c(log(y+2)), ncol=1)  #link function 
  one = matrix(rep(1,each=N), nrow=N)  #(Nx1) vector of ones
  Vmat = diag(c(exp(eta)), ncol=N, nrow=N)  #(NxN) diagonal weight matrix for the weighted least squares
  z.weighted.mean = ( (one%*%t(one)%*%Vmat%*%eta)/N ) 
  z0 = eta - z.weighted.mean
  X.weighted.mean = ( (one%*%t(one)%*%Vmat%*%X)/N )
  X0 = X - X.weighted.mean
  W.plus =  P.plus = matrix(0, nrow=P, ncol=A)  #(PxA) (deflated) X.weights and X.loadings
  #matrices											  
  q.plus_vec = c.plus_vec = matrix(0, nrow=M, ncol=A)  #(1xA) y.loadings and y.weights vectors
  T.plus = matrix(0, nrow=N, ncol=A)  #(NxA) X.scores matrix 
  b.spls_glm  = mod.SPLS(X,y,A,lambdaY,lambdaX,eps)$B.mat
  #soft thresholding penalization
  g.vec = function(u, lambda=NULL, ...)
  {
    sign(u) * ifelse ((abs(u)-lambda>0),(abs(u)-lambda), 0)
  }   
  #generalized inverse adapted from "MASS" package
  ginv = function (X, tol=sqrt(.Machine$double.eps), ...) 
  {
    if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X))) {
      stop("'X' must be a numeric or complex matrix")
    }
    if (!is.matrix(X)){ 
      X = as.matrix(X)
    }
    X.svd = svd(X)
    
    if (is.complex(X)){
      X.svd$u = Conj(X.svd$u)
    }
    Positive = X.svd$d > max(tol * X.svd$d[1L], 0)
    
    if (all(Positive)){ 
      X.svd$v %*% (1/X.svd$d * t(X.svd$u))
    }
    else if (!any(Positive)) {
      array(0, dim(X)[2L:1L])
    }
    else X.svd$v[, Positive, drop = FALSE] %*% ((1/X.svd$d[Positive]) * t(X.svd$u[, Positive, drop = FALSE]))
  }
  #2. Weighted PLS 
  repeat
  {
    #steps 2 to 9
    X.a = X0
    z.a = z0   
    for (a in 1:A)
    {       
      S = t(X.a) %*% Vmat %*% z.a
      w.a = svd(S)$u[,1]  #(Px1) X.weight
      c.a = svd(S)$v[,1]  #(Mx1) Y.weight     
      w.a.old = w.a
      c.a.old = c.a      
      conv.w.a = conv.c.a = 99
      while (conv.w.a > eps & conv.c.a > eps)
      {
        w.a.new = g.vec(u=S%*%c.a.old, lambda=lambdaX)
        if (!all(w.a.new==0)) w.a.new = w.a.new / as.numeric(sqrt(t(w.a.new)%*%w.a.new))
        c.a.new = g.vec(u=t(S)%*%w.a.old, lambda=lambdaY)             
        if (!all(c.a.new==0)) c.a.new = c.a.new / as.numeric(sqrt(t(c.a.new)%*%c.a.new))
        if (all(w.a.new==0)) conv.w.a = 0 else { old.non.zero = (c(w.a.old)!=0)
                                                 w.a.old1 = w.a.old[old.non.zero]
                                                 w.a.new1 = w.a.new[old.non.zero]
                                                 if (length(w.a.old1)>0)
                                                   conv.w.a = sum(abs(w.a.new1 - w.a.old1)/abs(w.a.old1))
                                                 else conv.w.a = 0                                         
        }
        if (all(c.a.new==0)) conv.c.a = 0 else { old.non.zero = (c(c.a.old)!=0)
                                                 c.a.old1 = c.a.old[old.non.zero]
                                                 c.a.new1 = c.a.new[old.non.zero]
                                                 if (length(c.a.old1)>0)
                                                   conv.c.a = sum(abs(c.a.new1 - c.a.old1)/abs(c.a.old1))
                                                 else conv.c.a = 0
        }
        w.a = w.a.new                                                                              
        w.a.old = w.a.new      
        c.a = c.a.new                                                                             
        c.a.old = c.a.new
      }				 
      t.a = X.a %*% w.a
      if (!all(t.a==0))  t.a = t.a / as.numeric(sqrt(t(X.a%*%w.a)%*%X.a%*%w.a))                
      p.a = t(X.a) %*% Vmat %*% t.a /  as.numeric(t(t.a)%*%Vmat%*%t.a)    #(Px1) X.loading
      q.a = t(z.a) %*% Vmat %*% t.a /  as.numeric(t(t.a)%*%Vmat%*%t.a)     #(Mx1) Y.loading
      
      #update X.a and Y.a
      X.a = X.a - t.a%*%t(p.a)
      z.a = z.a - t.a%*%t(q.a) 
      
      W.plus[, a] = w.a   #(PxA) X.weights matrix  
      c.plus_vec[, a] = c.a   #(1xA) y.weights vector
      q.plus_vec[, a] = q.a   #(1xA) y.loadings vector
      P.plus[, a] = p.a   #(PxA) X.loadings matrix                  
    }    
    dimnames(P.plus) = dimnames(W.plus) = list(colnames(X), paste("Comp", 1:A, seq=""))
    dimnames(q.plus_vec) = dimnames(c.plus_vec) = list(colnames(y), paste("Comp", 1:A, seq=""))
    
    R.plus = W.plus %*% ginv(t(P.plus)%*%W.plus)  #(PxA) transformed X.weights matrix
    if (any(is.na(ginv(t(P.plus)%*%W.plus)))) {
      break
    }
    else {
      T.plus = X0 %*% R.plus  #(NxA) X.scores matrix  
      b.spls_glm.old = b.spls_glm
      b.spls_glm = R.plus %*% t(q.plus_vec)  #(Px1) PLS-GLM coefficient vector 
      eta.hat = (X0 %*% b.spls_glm) + z.weighted.mean #expected y-values 
    }        
    rmsep = sqrt(sum((eta-eta.hat)^2)/N)  #Root Mean Squared Error of Prediction (RMSEP)     
    dimnames(T.plus) = list(rownames(X), paste("Comp", 1:A, seq=""))           	 
    
    #3. Update eta 
    eta = eta.hat	
    #4. Update Vmat and z0
    Vmat = diag(c(exp(eta)), ncol=N, nrow=N) 
    z0 = eta + diag(c(1/exp(eta)), ncol=N, nrow=N) %*% (y-exp(eta))  #the linearized form     
    #5-6. Check that changes are sufficiently small, else re-run from step 2  
    check.b.spls_glm = sum((b.spls_glm) - (b.spls_glm.old)/(b.spls_glm))
    if (check.b.spls_glm=="NaN") check.b.spls_glm = 0
    if ( check.b.spls_glm < eps )           		   
      break
    else
    {
      W.plus = W.plus   #(PxA) X.weights matrix       
      P.plus = P.plus   #(PxA) X.loadings matrix
      c.plus_vec = c.plus_vec   #(1xA) y.weights vector
      q.plus_vec = q.plus_vec   #(1xA) y.loadings vector
      R.plus = R.plus    #(PxA) transformed X.weights matrix 
      T.plus = T.plus    #(NxA) X.scores matrix
      b.spls_glm = b.spls_glm  #(Px1) PLS-GLM coefficient vector 
    }
  }    		
  #selecting non-zero values
  X.select = which(!R.plus[,1]==0, arr.ind=TRUE)     
  list(X.loadings=P.plus, y.weights=c.plus_vec, y.loadings=q.plus_vec, X.weights.trans=R.plus, 
       X.scores=T.plus, b.spls_glm_vec=b.spls_glm, RMSEP=rmsep, X.select=X.select, 
       Y.select="Since M=1, there's no need to perform variable selection on the Y-variable")
}   

#' Sparse Partial Least Squares-Generalized Linear Model (SPLS-GLM) algorithm for Binomial y 
#' 
#' Takes in a set of predictor variables and a set of response variables and gives the SPLS-GLM parameters.
#' @param X A (NxP) predictor matrix 
#' @param y A (Nx1) Binomial-distributed response vector
#' @param A The number of PLS components 
#' @param lambdaY A value for the penalty parameters for the soft-thresholding penalization function for Y-weights 
#' @param lambdaX A value for the penalty parameters for the soft-thresholding penalization function for X-weights 
#' @param eps Cut off value for convergence step
#' @param ... Other arguments. Currently ignored
#' @return The SPLS-GLM parameters of D=[X y]
#' @examples
#' if(require(MASS))
#' data(Pima.tr, package="MASS")
#' X = as.matrix(cbind(Pima.tr[,1:7]))  
#' dimnames(X) = list(1:nrow(X), colnames(X))
#' y = as.matrix(as.numeric(Pima.tr$type)-1, ncol=1)
#' #0=No and 1=Yes
#' dimnames(y) = list(1:nrow(y), paste("type"))
#' SPLS.binomial.GLM(scale(X), scale(y), A=2, lambdaY=0, lambdaX=0.96, eps=1e-3)
#' #lambdaX and lambdaY value are determined using function opt.penalty.values 
#' #for more details, see opt.penalty.values help file
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export 
SPLS.binomial.GLM = function(X, y, A, lambdaY, lambdaX, eps=0.001, ...)
{
  #For M=1 
  #y is Binomial distributed
  #A=initial number of components 
  #lambdaY & lambdaX are the value for the penalty parameters for the soft-thresholding penalization functions
  #1. Initialization  
  b.spls_glm = mod.SPLS(X,y,A,lambdaY,lambdaX,eps)$B.mat
  eta = matrix(c(X%*%b.spls_glm), ncol=1)  #link function 
  eta.old =eta
  N = nrow(X)
  P = ncol(X)
  M = 1 #ncol(y)
  one = matrix(rep(1,each=N), nrow=N)  #(Nx1) vector of ones
  pii = exp(X%*%b.spls_glm) / (1+exp(X%*%b.spls_glm))
  Vmat = diag(c(N*pii*(1-pii)), ncol=N, nrow=N)  #(NxN) diagonal weight matrix for the weighted least squares
  z.weighted.mean = ( (one%*%t(one)%*%Vmat%*%eta)/N ) 
  z0 = eta - z.weighted.mean
  X.weighted.mean = ( (one%*%t(one)%*%Vmat%*%X)/N )
  X0 = X - X.weighted.mean
  W.plus =  P.plus = matrix(0, nrow=P, ncol=A)  #(PxA) (deflated) X.weights and X.loadings
  #matrices 
  q.plus_vec = c.plus_vec = matrix(0, nrow=M, ncol=A)  #(1xA) y.loadings and y.weights vectors
  T.plus = matrix(0, nrow=N, ncol=A)  #(NxA) X.scores matrix  
  #soft thresholding penalization
  g.vec = function(u, lambda=NULL, ...)
  {
    sign(u) * ifelse ((abs(u)-lambda>0),(abs(u)-lambda), 0)
  }   
  #generalized inverse adapted from "MASS" package
  ginv = function (X, tol=sqrt(.Machine$double.eps), ...) 
  {
    if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X))) {
      stop("'X' must be a numeric or complex matrix")
    }
    if (!is.matrix(X)){ 
      X = as.matrix(X)
    }
    X.svd = svd(X)
    
    if (is.complex(X)){
      X.svd$u = Conj(X.svd$u)
    }
    Positive = X.svd$d > max(tol * X.svd$d[1L], 0)
    
    if (all(Positive)){ 
      X.svd$v %*% (1/X.svd$d * t(X.svd$u))
    }
    else if (!any(Positive)) {
      array(0, dim(X)[2L:1L])
    }
    else X.svd$v[, Positive, drop = FALSE] %*% ((1/X.svd$d[Positive]) * t(X.svd$u[, Positive, drop = FALSE]))
  }
  #2. Weighted PLS 
  repeat
  {
    #steps 2 to 9
    X.a = X0
    z.a = z0   
    for (a in 1:A)
    {       
      S = t(X.a) %*% Vmat %*% z.a
      w.a = svd(S)$u[,1]  #(Px1) X.weight
      c.a = svd(S)$v[,1]  #(Mx1) Y.weight     
      w.a.old = w.a
      c.a.old = c.a      
      conv.w.a = conv.c.a = 99
      while (conv.w.a > eps & conv.c.a > eps)
      {
        w.a.new = g.vec(u=S%*%c.a.old, lambda=lambdaX)
        if (!all(w.a.new==0)) w.a.new = w.a.new / as.numeric(sqrt(t(w.a.new)%*%w.a.new))
        c.a.new = g.vec(u=t(S)%*%w.a.old, lambda=lambdaY)             
        if (!all(c.a.new==0)) c.a.new = c.a.new / as.numeric(sqrt(t(c.a.new)%*%c.a.new))
        if (all(w.a.new==0)) conv.w.a = 0 else { old.non.zero = (c(w.a.old)!=0)
                                                 w.a.old1 = w.a.old[old.non.zero]
                                                 w.a.new1 = w.a.new[old.non.zero]
                                                 if (length(w.a.old1)>0)
                                                   conv.w.a = sum(abs(w.a.new1 - w.a.old1)/abs(w.a.old1))
                                                 else conv.w.a = 0                                         
        }
        if (all(c.a.new==0)) conv.c.a = 0 else { old.non.zero = (c(c.a.old)!=0)
                                                 c.a.old1 = c.a.old[old.non.zero]
                                                 c.a.new1 = c.a.new[old.non.zero]
                                                 if (length(c.a.old1)>0)
                                                   conv.c.a = sum(abs(c.a.new1 - c.a.old1)/abs(c.a.old1))
                                                 else conv.c.a = 0
        }
        w.a = w.a.new                                                                              
        w.a.old = w.a.new      
        c.a = c.a.new                                                                             
        c.a.old = c.a.new    
      }			
      t.a = X.a %*% w.a
      if (!all(t.a==0))  t.a = t.a / as.numeric(sqrt(t(X.a%*%w.a)%*%X.a%*%w.a))                
      p.a = t(X.a) %*% Vmat %*% t.a /  as.numeric(t(t.a)%*%Vmat%*%t.a)    #(Px1) X.loading
      q.a = t(z.a) %*% Vmat %*% t.a /  as.numeric(t(t.a)%*%Vmat%*%t.a)     #(Mx1) Y.loading
      
      #update X.a and Y.a
      X.a = X.a - t.a%*%t(p.a)
      z.a = z.a - t.a%*%t(q.a) 
      
      W.plus[, a] = w.a   #(PxA) X.weights matrix  
      c.plus_vec[, a] = c.a   #(1xA) y.weights vector
      q.plus_vec[, a] = q.a   #(1xA) y.loadings vector
      P.plus[, a] = p.a   #(PxA) X.loadings matrix   
    } 	    
    dimnames(P.plus) = dimnames(W.plus) = list(colnames(X), paste("Comp", 1:A, seq=""))
    dimnames(q.plus_vec) = dimnames(c.plus_vec) = list(colnames(y), paste("Comp", 1:A, seq=""))
    
    R.plus = W.plus %*% ginv(t(P.plus)%*%W.plus)  #(PxA) transformed X.weights matrix
    if (any(is.na(ginv(t(P.plus)%*%W.plus)))) {
      break
    }
    else {
      T.plus = X0 %*% R.plus  #(NxA) X.scores matrix  
      b.spls_glm.old = b.spls_glm
      b.spls_glm = R.plus %*% t(q.plus_vec)
      eta.hat = (X0 %*% b.spls_glm) + z.weighted.mean #expected y-values 
    }        
    rmsep = sqrt(sum((y-eta.hat)^2)/N)  #Root Mean Squared Error of Prediction (RMSEP)     
    dimnames(T.plus) = list(rownames(X), paste("Comp", 1:A, seq=""))    
    #3. Update eta 
    eta = eta.hat 
    #4. Update Vmat and z0
    pii = exp(eta) / (1+exp(eta))	   
    Vmat = diag(c(N*pii*(1-pii)), ncol=N, nrow=N)
    z0 = eta + ((y-(N*pii))/(N*pii*(1-pii)))  #the linearized form 
    #5-6. Check that changes are sufficiently small, else re-run from step 2  
    check.b.spls_glm = sum((b.spls_glm) - (b.spls_glm.old)/(b.spls_glm))
    if (check.b.spls_glm=="NaN") check.b.spls_glm = 0
    if ( check.b.spls_glm < eps )              
      break
    else
    {
      W.plus = W.plus   #(PxA) X.weights matrix       
      P.plus = P.plus   #(PxA) X.loadings matrix
      c.plus_vec = c.plus_vec   #(1xA) y.weights vector
      q.plus_vec = q.plus_vec   #(1xA) y.loadings vector
      R.plus = R.plus    #(PxA) transformed X.weights matrix 
      T.plus = T.plus    #(NxA) X.scores matrix
      b.spls_glm =  b.spls_glm  #(Px1) PLSR coefficient vector 
    }
  }   		
  #selecting non-zero values
  X.select = which(!R.plus[,1]==0, arr.ind=TRUE)
  list(X.loadings=P.plus, y.weights=c.plus_vec, y.loadings=q.plus_vec, X.weights.trans=R.plus, 
       X.scores=T.plus, b.spls_glm_vec=b.spls_glm, RMSEP=rmsep, X.select=X.select, 
       Y.select="Since M=1, there's no need to perform variable selection on the Y-variable")
}
   
#' Choosing a value for penalty parameters lambdaX and lambdaY for the Sparse Partial Least Squares (SPLS) and Sparse Partial Least Squares-Generalized Linear Model (SPLS-GLM) analyses    
#' 
#' Gives the value of the penalty parameters (lambdaX,lambdaY) having the minimum RMSEP value.
#' @param X A (NxP) predictor matrix 
#' @param Y A (NxM) response matrix; can also be a vector in the case of the SPLS-GLM
#' @param A The number of Partial Least Squares (PLS) components 
#' @param algorithm Any of the SPLS or SPLS-GLM algorithms ("mod.SPLS", "SPLS.GLM", "SPLS.binomial.GLM") 
#' @param eps Cut off value for convergence step 
#' @param from.value.X starting value for lambdaX
#' @param to.value.X ending value for lambdaX
#' @param from.value.Y starting value for lambdaY
#' @param to.value.Y ending value for lambdaY
#' @param lambdaY.len length of lambdaY value
#' @param lambdaX.len length of lambdaX value
#' @param ... Other arguments. Currently ignored
#' @return the value of the penalty parameters (lambdaX,lambdaY) having the minimum RMSEP value, as well as the RMSEP values obtained when the lambdaX and lambdaY values were paired together
#' @examples
#' if(require(chemometrics))
#' data(ash, package="chemometrics")
#' X1 = as.matrix(ash[,10:17], ncol=8)  
#' Y1 = as.matrix(ash$SOT)   
#' colnames(Y1) = paste("SOT")
#' #choosing a value for the penalty parameters lambdaY and lambdaX for this data 
#' opt.penalty.values(X=scale(X1), Y=scale(Y1), A=2, algorithm=mod.SPLS, eps=1e-5, 
#' from.value.X=0, to.value.X=500, from.value.Y=0, to.value.Y=0, lambdaY.len=1, lambdaX.len=100)
#' #thus, use lambdaX = 10.10 and lambdaY = 0 for the SPLS analysis of this data
#' 
#' #possum.mat data   
#' if(require(robustbase))
#' possum.mat 
#' y = as.matrix(possum.mat[,1], ncol=1) 
#' dimnames(y) = list(paste("S", 1:nrow(possum.mat), seq=""), "Diversity")
#' X = as.matrix(possum.mat[,2:14], ncol=13) 
#' dimnames(X) = list(paste("S", 1:nrow(possum.mat), seq=""), colnames(possum.mat[,2:14])) 
#' #choosing a value for the penalty parameters lambdaY and lambdaX for this data 
#' opt.penalty.values(X=scale(X), Y=scale(y), A=2, algorithm=SPLS.GLM, eps=1e-3, 
#' from.value.X=1, to.value.X=4, from.value.Y=0, to.value.Y=0, lambdaY.len=1, lambdaX.len=100)
#' #thus, use lambdaY = 0 and lambdaX = 3.3 for the (Poisson) SPLS-GLM analysis of this data
#' 
#' #Pima.tr data
#' if(require(MASS))
#' data(Pima.tr, package="MASS")
#' X = as.matrix(cbind(Pima.tr[,1:7]))  
#' dimnames(X) = list(1:nrow(X), colnames(X))
#' y = as.matrix(as.numeric(Pima.tr$type)-1, ncol=1)
#' #0=No and 1=Yes
#' dimnames(y) = list(1:nrow(y), paste("type"))
#' #choosing a value for the penalty parameters lambdaY and lambdaX for this data 
#' opt.penalty.values(X=scale(X), Y=scale(y), A=2, algorithm=SPLS.binomial.GLM, eps=1e-3, 
#' from.value.X=0, to.value.X=95, from.value.Y=0, to.value.Y=0, lambdaY.len=1, lambdaX.len=100)
#' #thus, use lambdaY = 0 and lambdaX = 0.96 for the (Binomial) SPLS-GLM analysis of this data
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export 
opt.penalty.values = function(X, Y, A, algorithm=NULL, eps, from.value.X, to.value.X, from.value.Y, to.value.Y, lambdaY.len, lambdaX.len, ...)  
{    
  lambdaY = seq(from=from.value.Y, to=to.value.Y, length.out=lambdaY.len)
  lambdaX = seq(from=from.value.X, to=to.value.X, length.out=lambdaX.len)
  L.mat = expand.grid(lambdaY,lambdaX)
  dimnames(L.mat) = list(1:nrow(L.mat),paste(c("lambdaY","lambdaX")))  
  rmsep = array(NA, dim=c(1,nrow(L.mat))) 
  for (i in 1:nrow(L.mat))
  {
    flush.console() 
    rmsep[,i] = algorithm(X, Y, A, lambdaY=L.mat[i,1], lambdaX=L.mat[i,2], eps)$RMSEP
    dimnames(rmsep) = list(paste("RMSEP.value"), paste("Combo", 1:nrow(L.mat)))
  }
  options(digits=7)
  To.use = cbind(L.mat,t(rmsep))	
  min.RMSEP.value = min(To.use[,3])	
  i.final = which.min(To.use[,3])
  lambdaY.new = To.use[i.final,1]
  lambdaX.new = To.use[i.final,2]
  list(RMSEP.values=To.use, min.RMSEP.value=min.RMSEP.value, lambdaY.to.use=lambdaY.new, lambdaX.to.use=lambdaX.new)
}

# PLS biplot for SPLS 
#' The Partial Least Squares (PLS) biplot for Sparse Partial Least Squares (SPLS) 
#' 
#' Takes in a set of predictor variables and a set of response variables and produces a PLS biplot for the SPLS.
#' @param X A (NxP) predictor matrix 
#' @param Y A (NxM) response matrix
#' @param algorithm The SPLS algorithm 
#' @param eps Cut off value for convergence step 
#' @param lambdaY A value for the penalty parameters for the soft-thresholding penalization function for Y-weights 
#' @param lambdaX A value for the penalty parameters for the soft-thresholding penalization function for X-weights 
#' @param ax.tickvec.X tick marker length for each X-variable axis in the biplot
#' @param ax.tickvec.Y tick marker length for the Y-variable axis in the biplot
#' @param ... Other arguments. Currently ignored
#' @return The PLS biplot of a SPLS of D=[X Y] with some parameters
#' @examples   
#' if(require(robustbase))
#' data(toxicity, package="robustbase")
#' Y1 = as.matrix(cbind(toxicity$toxicity))  
#' dimnames(Y1) = list(paste(1:nrow(Y1)), "toxicity")
#' X1 = as.matrix(cbind(toxicity[,2:10]))  
#' rownames(X1) = paste(1:nrow(X1))
#' #choosing a value for the penalty parameters lambdaY and lambdaX for this data 
#' main2 = opt.penalty.values(X=scale(X1), Y=scale(Y1), A=2, algorithm=mod.SPLS, eps=1e-5, 
#' from.value.X=0, to.value.X=10, from.value.Y=0, to.value.Y=0, lambdaY.len=1, lambdaX.len=100)
#' min.RMSEP.value = main2$min.RMSEP.value
#' lambdaY.to.use = main2$lambdaY.to.use
#' lambdaX.to.use = main2$lambdaX.to.use
#' list(lambdaY.to.use=lambdaY.to.use, lambdaX.to.use=lambdaX.to.use, min.RMSEP.value=min.RMSEP.value)
#' #SPLS analysis
#' main3 = mod.SPLS(X=scale(X1), Y=scale(Y1), A=2, 
#' lambdaY=lambdaY.to.use, lambdaX=lambdaX.to.use,
#' eps=1e-5)
#' X.to.use = main3$X.select
#' Y.to.use = main3$Y.select
#' X.new = as.matrix(X1[,X.to.use])
#' colnames(X.new)
#' Y.new = as.matrix(Y1[,Y.to.use])
#' colnames(Y.new) = colnames(Y1)
#' colnames(Y.new) 
#' SPLS.biplot(X.new, Y.new, algorithm=mod.SPLS, lambdaY=lambdaY.to.use, lambdaX=lambdaX.to.use, 
#' eps=1e-5, ax.tickvec.X=rep(3,ncol(X.new)), ax.tickvec.Y=rep(4,ncol(Y.new)))
#' 
#' #ash data 
#' if(require(chemometrics))
#' data(ash, package="chemometrics")
#' X1 = as.matrix(ash[,10:17], ncol=8)  
#' Y1 = as.matrix(ash$SOT)   
#' colnames(Y1) = paste("SOT")
#' #choosing a value for the penalty parameters lambdaY and lambdaX for this data 
#' main2 = opt.penalty.values(X=scale(X1), Y=scale(Y1), A=2, algorithm=mod.SPLS, eps=1e-5, 
#' from.value.X=0, to.value.X=500, from.value.Y=0, to.value.Y=0, lambdaY.len=1, lambdaX.len=100)
#' min.RMSEP.value = main2$min.RMSEP.value
#' lambdaY.to.use = main2$lambdaY.to.use
#' lambdaX.to.use = main2$lambdaX.to.use
#' list(lambdaY.to.use=lambdaY.to.use, lambdaX.to.use=lambdaX.to.use, min.RMSEP.value=min.RMSEP.value) 
#' #SPLS analysis 
#' main3 = mod.SPLS(X=scale(X1), Y=scale(Y1), A=2, 
#' lambdaY=lambdaY.to.use, lambdaX=lambdaX.to.use,
#' eps=1e-5)
#' X.to.use = main3$X.select
#' Y.to.use = main3$Y.select
#' X.new = as.matrix(X1[,X.to.use])
#' colnames(X.new)  #P=6
#' colnames(X1)  #P=8
#' Y.new = as.matrix(Y1[,Y.to.use])
#' colnames(Y.new) = colnames(Y1)
#' colnames(Y.new)
#' SPLS.biplot(X.new, Y.new, algorithm=mod.SPLS, lambdaY=lambdaY.to.use, lambdaX=lambdaX.to.use, 
#' eps=1e-5, ax.tickvec.X=rep(1,ncol(X.new)), ax.tickvec.Y=rep(5,ncol(Y.new)))
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export
SPLS.biplot = function(X, Y, algorithm=NULL, eps, lambdaY=NULL, lambdaX=NULL, ax.tickvec.X=NULL, ax.tickvec.Y=NULL, ... ) 
{    
  options(digits=3)
  D = cbind(X,Y)
  X.mean = apply(X, 2, mean)  #column means of X or colMeans(X)
  Y.mean = apply(Y, 2, mean)  #column means of Y or colMeans(Y)
  X.std = apply(X, 2, sd)  #column standard deviations of X 
  Y.std = apply(Y, 2, sd)  #column standard deviations of Y   
  X.scal = scale(X, center=TRUE, scale=TRUE)  #scaled X
  Y.scal = scale(Y, center=TRUE, scale=TRUE)  #scaled Y  	
  
  #biplot	
  A = 2
  main = algorithm(X.scal,Y.scal,A,lambdaY,lambdaX,eps)
  Tmat = main$X.scores
  Pmat = main$X.loadings
  Qmat = main$Y.loadings
  Rmat = main$X.weights.trans	
  Bmat = Rmat %*% t(Qmat)  #(PxM) estimated SPLS coefficients matrix	 
  dimnames(Bmat) = list(colnames(X), colnames(Y))	    
  
  #overall quality of approximation		 
  X.hat = ((Tmat %*% t(Pmat)) * rep(X.std,each=nrow(X))) + rep(X.mean,each=nrow(X))
  dimnames(X.hat) = dimnames(X)
  Y.hat = ((Tmat %*% t(Qmat)) * rep(Y.std,each=nrow(Y))) + rep(Y.mean,each=nrow(Y))
  dimnames(Y.hat) = dimnames(Y)
  D.hat = cbind(X.hat, Y.hat)	   
  overall.quality = sum(diag(t(D.hat)%*%D.hat))  /  sum(diag(t(D)%*%D))
  
  #axes predictivity
  axis.pred = (diag(t(D.hat)%*%D.hat)) / (diag(t(D)%*%D))	    
  names(axis.pred) = colnames(D)		 
  
  #biplot points
  Z = rbind(Tmat, Rmat)  #((N+P)x2) biplot points
  dimnames(Z) = list(c(rownames(X),paste("b",1:ncol(X),sep="")), NULL)
  sample.style = biplot.sample.control(2, label=TRUE)
  sample.style$col = c("red", "purple")
  sample.style$pch = c(15,16)
  Gmat = cbind(c(rep(1,nrow(X)),rep(0,ncol(X))), c(rep(0,nrow(X)),rep(1,ncol(X))))
  dimnames(Gmat) = list(NULL, c("samples","coefficients"))
  classes = 1:2	 
  
  #axes direction    
  #X.hat=TP'
  X.axes.direction = (1 / (diag(Pmat%*%t(Pmat)))) * Pmat       	   
  #Y.hat=TQ'
  Y.axes.direction = (1 / (diag(Qmat%*%t(Qmat)))) * Qmat
  #B.spls=RQ'
  B.axes.direction = (1 / (diag(Qmat%*%t(Qmat)))) * Qmat	
  
  #calibration of axes	 	
  z.axes.X = lapply(1:ncol(X), calibrate.axis, X, X.mean, X.std, X.axes.direction, 
                    1:ncol(X), ax.tickvec.X, rep(0,ncol(X)), rep(0,ncol(X)), NULL)  
  z.axes.Y = lapply(1:ncol(Y), calibrate.axis, Y, Y.mean, Y.std, Y.axes.direction, 
                    1:ncol(Y), ax.tickvec.Y, rep(0,ncol(Y)), rep(0,ncol(Y)), NULL)
  z.axes.B = lapply(1:ncol(Bmat), calibrate.axis, Y.scal, rep(0,ncol(Y)), rep(1,ncol(Y)),
                    B.axes.direction, 1:ncol(Bmat), rep(4,ncol(Bmat)), rep(0,ncol(Bmat)), 
                    rep(0,ncol(Bmat)), NULL)	
  z.axes = vector("list", ncol(X)+ncol(Y)+ncol(Bmat))
  for (i in 1:ncol(X)) 
  {
    z.axes[[i]] = z.axes.X[[i]]
  }
  for (i in 1:ncol(Y)) 
  {
    z.axes[[i+ncol(X)]] = z.axes.Y[[i]]
  }
  for (i in 1:ncol(Bmat)) 
  {
    z.axes[[i+ncol(X)+ncol(Y)]] = z.axes.B[[i]]  
  }			  
  
  #biplot axes 
  ax.style = biplot.ax.control(ncol(X)+ncol(Y)+ncol(Bmat), c(colnames(X), colnames(Y),
                                                             colnames(Bmat)))
  #for X	 
  ax.style$col[1:ncol(X)] = "blue"
  ax.style$label.col[1:ncol(X)] = "blue"
  ax.style$tick.col[1:ncol(X)] = "blue"
  ax.style$tick.label.col[1:ncol(X)] = "blue"  
  #for Y	 
  ax.style$col[ncol(X)+1:ncol(Y)] = "black"
  ax.style$label.col[ncol(X)+1:ncol(Y)] = "black"
  ax.style$tick.col[ncol(X)+1:ncol(Y)] = "black"
  ax.style$tick.label.col[ncol(X)+1:ncol(Y)] = "black"  
  #for B.spls
  ax.style$col[ncol(X)+ncol(Y)+1:ncol(Bmat)] = "black"
  ax.style$label.col[ncol(X)+ncol(Y)+1:ncol(Bmat)] = "black" 
  ax.style$tick.col[ncol(X)+ncol(Y)+1:ncol(Bmat)] = "purple"
  ax.style$tick.label.col[ncol(X)+ncol(Y)+1:ncol(Bmat)] = "purple"   
  draw.biplot(Z, G=Gmat, classes=classes, sample.style=sample.style, z.axes=z.axes, 
              ax.style=ax.style, VV=rbind(X.axes.direction,Y.axes.direction))
  list(overall.quality=overall.quality, axis.pred=axis.pred, D.hat=D.hat, Bmat=Bmat)  	 
}

#' The Partial Least Squares (PLS) biplot for Sparse Partial Least Squares (SPLS), with the labels of the samples, coefficient points and tick markers excluded
#' 
#' Takes in a set of predictor variables and a set of response variables and produces a PLS biplot for the SPLS with the labels of the samples, coefficient points and tick markers excluded.
#' @param X A (NxP) predictor matrix 
#' @param Y A (NxM) response matrix
#' @param algorithm The SPLS algorithm 
#' @param eps Cut off value for convergence step 
#' @param lambdaY A value for the penalty parameters for the soft-thresholding penalization function for Y-weights 
#' @param lambdaX A value for the penalty parameters for the soft-thresholding penalization function for X-weights 
#' @param ax.tickvec.X tick marker length for each X-variable axis in the biplot
#' @param ax.tickvec.Y tick marker length for the Y-variable axis in the biplot
#' @param ... Other arguments. Currently ignored
#' @return The PLS biplot of a SPLS of D=[X Y] with some parameters
#' @examples
#' if(require(robustbase))
#' data(toxicity, package="robustbase")
#' Y1 = as.matrix(cbind(toxicity$toxicity))  
#' dimnames(Y1) = list(paste(1:nrow(Y1)), "toxicity")
#' X1 = as.matrix(cbind(toxicity[,2:10]))  
#' rownames(X1) = paste(1:nrow(X1))
#' #choosing a value for the penalty parameters lambdaY and lambdaX for this data 
#' main2 = opt.penalty.values(X=scale(X1), Y=scale(Y1), A=2, algorithm=mod.SPLS, eps=1e-5,
#' from.value.X=0, to.value.X=10, from.value.Y=0, to.value.Y=0, lambdaY.len=1, lambdaX.len=100)
#' min.RMSEP.value = main2$min.RMSEP.value
#' lambdaY.to.use = main2$lambdaY.to.use
#' lambdaX.to.use = main2$lambdaX.to.use
#' list(lambdaY.to.use=lambdaY.to.use, lambdaX.to.use=lambdaX.to.use, min.RMSEP.value=min.RMSEP.value)
#' #SPLS analysis
#' main3 = mod.SPLS(X=scale(X1), Y=scale(Y1), A=2, 
#' lambdaY=lambdaY.to.use, lambdaX=lambdaX.to.use, 
#' eps=1e-5)
#' X.to.use = main3$X.select
#' Y.to.use = main3$Y.select
#' X.new = as.matrix(X1[,X.to.use])
#' Y.new = as.matrix(Y1[,Y.to.use])
#' colnames(Y.new) = colnames(Y1)
#' SPLS.biplot_no_labels(X.new, Y.new, 
#' algorithm=mod.SPLS, lambdaY=lambdaY.to.use, 
#' lambdaX=lambdaX.to.use, eps=1e-5, 
#' ax.tickvec.X=rep(3,ncol(X.new)), 
#' ax.tickvec.Y=rep(4,ncol(Y.new)))
#' 
#' #ash data  
#' if(require(chemometrics))
#' data(ash, package="chemometrics")
#' X1 = as.matrix(ash[,10:17], ncol=8)  
#' Y1 = as.matrix(ash$SOT)   
#' colnames(Y1) = paste("SOT")
#' #choosing a value for the penalty parameters lambdaY and lambdaX for this data 
#' main2 = opt.penalty.values(X=scale(X1), Y=scale(Y1), A=2, algorithm=mod.SPLS, eps=1e-5, 
#' from.value.X=0, to.value.X=500, from.value.Y=0, to.value.Y=0, lambdaY.len=1, lambdaX.len=100)
#' min.RMSEP.value = main2$min.RMSEP.value
#' lambdaY.to.use = main2$lambdaY.to.use
#' lambdaX.to.use = main2$lambdaX.to.use
#' list(lambdaY.to.use=lambdaY.to.use, lambdaX.to.use=lambdaX.to.use, min.RMSEP.value=min.RMSEP.value) 
#' #SPLS analysis
#' main3 = mod.SPLS(X=scale(X1), Y=scale(Y1), A=2, 
#' lambdaY=lambdaY.to.use, lambdaX=lambdaX.to.use, 
#' eps=1e-5)
#' X.to.use = main3$X.select
#' Y.to.use = main3$Y.select
#' X.new = as.matrix(X1[,X.to.use])
#' colnames(X.new)  #P=6
#' colnames(X1)  #P=8
#' Y.new = as.matrix(Y1[,Y.to.use])
#' colnames(Y.new) = colnames(Y1)
#' colnames(Y.new)
#' SPLS.biplot_no_labels(X.new, Y.new, 
#' algorithm=mod.SPLS, lambdaY=lambdaY.to.use,
#' lambdaX=lambdaX.to.use, eps=1e-5, 
#' ax.tickvec.X=rep(1,ncol(X.new)), 
#' ax.tickvec.Y=rep(5,ncol(Y.new)))
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export
SPLS.biplot_no_labels = function(X, Y, algorithm=NULL, eps, lambdaY=NULL, lambdaX=NULL, ax.tickvec.X=NULL, ax.tickvec.Y=NULL, ... ) 
{    
  options(digits=3)
  D = cbind(X,Y)
  X.mean = apply(X, 2, mean)  #column means of X or colMeans(X)
  Y.mean = apply(Y, 2, mean)  #column means of Y or colMeans(Y)
  X.std = apply(X, 2, sd)  #column standard deviations of X 
  Y.std = apply(Y, 2, sd)  #column standard deviations of Y   
  X.scal = scale(X, center=TRUE, scale=TRUE)  #scaled X
  Y.scal = scale(Y, center=TRUE, scale=TRUE)  #scaled Y  	
  
  #biplot	
  A = 2
  main = algorithm(X.scal,Y.scal,A,lambdaY,lambdaX,eps)
  Tmat = main$X.scores
  Pmat = main$X.loadings
  Qmat = main$Y.loadings
  Rmat = main$X.weights.trans		
  Bmat = Rmat %*% t(Qmat)  #(PxM) estimated SPLS coefficients matrix	 
  dimnames(Bmat) = list(colnames(X), colnames(Y))	    
  
  #overall quality of approximation		 
  X.hat = ((Tmat %*% t(Pmat)) * rep(X.std,each=nrow(X))) + rep(X.mean,each=nrow(X))
  dimnames(X.hat) = dimnames(X)
  Y.hat = ((Tmat %*% t(Qmat)) * rep(Y.std,each=nrow(Y))) + rep(Y.mean,each=nrow(Y))
  dimnames(Y.hat) = dimnames(Y)
  D.hat = cbind(X.hat, Y.hat)	   
  overall.quality = sum(diag(t(D.hat)%*%D.hat))  /  sum(diag(t(D)%*%D))
  
  #axes predictivity
  axis.pred = (diag(t(D.hat)%*%D.hat)) / (diag(t(D)%*%D))	    
  names(axis.pred) = colnames(D)		 
  
  #biplot points
  Z = rbind(Tmat, Rmat)  #((N+P)x2) biplot points
  dimnames(Z) = list(c(rownames(X),paste("b",1:ncol(X),sep="")), NULL)
  sample.style = biplot.sample.control(2, label=FALSE)
  sample.style$col = c("red", "purple")
  sample.style$pch = c(15,16)
  Gmat = cbind(c(rep(1,nrow(X)),rep(0,ncol(X))), c(rep(0,nrow(X)),rep(1,ncol(X))))
  dimnames(Gmat) = list(NULL, c("samples","coefficients"))
  classes = 1:2	 
  
  #axes direction    
  #X.hat=TP'
  X.axes.direction = (1 / (diag(Pmat%*%t(Pmat)))) * Pmat       	   
  #Y.hat=TQ'
  Y.axes.direction = (1 / (diag(Qmat%*%t(Qmat)))) * Qmat
  #B.spls=RQ'
  B.axes.direction = (1 / (diag(Qmat%*%t(Qmat)))) * Qmat	
  
  #calibration of axes	 	
  z.axes.X = lapply(1:ncol(X), calibrate.axis, X, X.mean, X.std, X.axes.direction, 
                    1:ncol(X), ax.tickvec.X, rep(0,ncol(X)), rep(0,ncol(X)), NULL)  
  z.axes.Y = lapply(1:ncol(Y), calibrate.axis, Y, Y.mean, Y.std, Y.axes.direction, 
                    1:ncol(Y), ax.tickvec.Y, rep(0,ncol(Y)), rep(0,ncol(Y)), NULL)
  z.axes.B = lapply(1:ncol(Bmat), calibrate.axis, Y.scal, rep(0,ncol(Y)), rep(1,ncol(Y)),
                    B.axes.direction, 1:ncol(Bmat), rep(4,ncol(Bmat)), rep(0,ncol(Bmat)), 
                    rep(0,ncol(Bmat)), NULL)	
  z.axes = vector("list", ncol(X)+ncol(Y)+ncol(Bmat))
  for (i in 1:ncol(X)) 
  {
    z.axes[[i]] = z.axes.X[[i]]
  }
  for (i in 1:ncol(Y)) 
  {
    z.axes[[i+ncol(X)]] = z.axes.Y[[i]]
  }
  for (i in 1:ncol(Bmat)) 
  {
    z.axes[[i+ncol(X)+ncol(Y)]] = z.axes.B[[i]]  
  }			  
  
  #biplot axes 
  ax.style = biplot.ax.control(ncol(X)+ncol(Y)+ncol(Bmat), c(colnames(X), colnames(Y),
                                                             colnames(Bmat)))
  ax.style$tick.col[1:ncol(X)] = "azure"
  ax.style$tick.label.col[1:ncol(X)] = "azure"
  ax.style$tick.col[ncol(X)+1:ncol(Y)] = "azure"
  ax.style$tick.label.col[ncol(X)+1:ncol(Y)] = "azure"
  ax.style$tick.col[ncol(X)+ncol(Y)+1:ncol(Bmat)] = "azure"
  ax.style$tick.label.col[ncol(X)+ncol(Y)+1:ncol(Bmat)] = "azure"
  #for X	
  ax.style$col[1:ncol(X)] = "blue"
  ax.style$label.col[1:ncol(X)] = "blue"
  #for Y	 
  ax.style$col[ncol(X)+1:ncol(Y)] = "black"
  ax.style$label.col[ncol(X)+1:ncol(Y)] = "black"
  #for B.spls
  ax.style$col[ncol(X)+ncol(Y)+1:ncol(Bmat)] = "black"
  ax.style$label.col[ncol(X)+ncol(Y)+1:ncol(Bmat)] = "black" 
  draw.biplot(Z, G=Gmat, classes=classes, sample.style=sample.style, z.axes=z.axes, 
              ax.style=ax.style, VV=rbind(X.axes.direction,Y.axes.direction))
  list(overall.quality=overall.quality, axis.pred=axis.pred, D.hat=D.hat, Bmat=Bmat)  	 
}

#' The Partial Least Squares (PLS) biplot for Sparse Partial Least Squares (SPLS), with the labels of the tick markers excluded
#' 
#' Takes in a set of predictor variables and a set of response variables and produces a PLS biplot for the SPLS with the labels of the tick markers excluded.
#' @param X A (NxP) predictor matrix 
#' @param Y A (NxM) response matrix
#' @param algorithm The SPLS algorithm 
#' @param eps Cut off value for convergence step 
#' @param lambdaY A value for the penalty parameters for the soft-thresholding penalization function for Y-weights 
#' @param lambdaX A value for the penalty parameters for the soft-thresholding penalization function for X-weights 
#' @param ax.tickvec.X tick marker length for each X-variable axis in the biplot
#' @param ax.tickvec.Y tick marker length for the Y-variable axis in the biplot
#' @param ... Other arguments. Currently ignored
#' @return The PLS biplot of a SPLS of D=[X Y] with some parameters
#' @examples
#' if(require(robustbase))
#' data(toxicity, package="robustbase")
#' Y1 = as.matrix(cbind(toxicity$toxicity))  
#' dimnames(Y1) = list(paste(1:nrow(Y1)), "toxicity")
#' X1 = as.matrix(cbind(toxicity[,2:10]))  
#' rownames(X1) = paste(1:nrow(X1))
#' #choosing a value for the penalty parameters lambdaY and lambdaX for this data 
#' main2 = opt.penalty.values(X=scale(X1), Y=scale(Y1), A=2, algorithm=mod.SPLS, eps=1e-5, 
#' from.value.X=0, to.value.X=10, from.value.Y=0, to.value.Y=0, lambdaY.len=1, lambdaX.len=100)
#' min.RMSEP.value = main2$min.RMSEP.value
#' lambdaY.to.use = main2$lambdaY.to.use
#' lambdaX.to.use = main2$lambdaX.to.use
#' list(lambdaY.to.use=lambdaY.to.use, lambdaX.to.use=lambdaX.to.use, min.RMSEP.value=min.RMSEP.value)
#' #SPLS analysis
#' main3 = mod.SPLS(X=scale(X1), Y=scale(Y1), A=2, 
#' lambdaY=lambdaY.to.use, lambdaX=lambdaX.to.use, 
#' eps=1e-5)
#' X.to.use = main3$X.select
#' Y.to.use = main3$Y.select
#' X.new = as.matrix(X1[,X.to.use])
#' Y.new = as.matrix(Y1[,Y.to.use])
#' colnames(Y.new) = colnames(Y1)
#' SPLS.biplot_no_ax.labels(X.new, Y.new, 
#' algorithm=mod.SPLS, lambdaY=lambdaY.to.use, 
#' lambdaX=lambdaX.to.use, 
#' eps=1e-5, ax.tickvec.X=rep(3,ncol(X.new)), 
#' ax.tickvec.Y=rep(4,ncol(Y.new)))
#' 
#' #ash data 
#' if(require(chemometrics)) 
#' data(ash, package="chemometrics")
#' X1 = as.matrix(ash[,10:17], ncol=8)  
#' Y1 = as.matrix(ash$SOT)   
#' colnames(Y1) = paste("SOT")
#' #choosing a value for the penalty parameters lambdaY and lambdaX for this data 
#' main2 = opt.penalty.values(X=scale(X1), Y=scale(Y1), A=2, algorithm=mod.SPLS, eps=1e-5, 
#' from.value.X=0, to.value.X=500, from.value.Y=0, to.value.Y=0, lambdaY.len=1, lambdaX.len=100)
#' min.RMSEP.value = main2$min.RMSEP.value
#' lambdaY.to.use = main2$lambdaY.to.use
#' lambdaX.to.use = main2$lambdaX.to.use
#' list(lambdaY.to.use=lambdaY.to.use, lambdaX.to.use=lambdaX.to.use, min.RMSEP.value=min.RMSEP.value) 
#' #SPLS analysis
#' main3 = mod.SPLS(X=scale(X1), Y=scale(Y1), A=2, 
#' lambdaY=lambdaY.to.use, lambdaX=lambdaX.to.use, 
#' eps=1e-5)
#' X.to.use = main3$X.select
#' Y.to.use = main3$Y.select
#' X.new = as.matrix(X1[,X.to.use])
#' colnames(X.new)  #P=6
#' colnames(X1)  #P=8
#' Y.new = as.matrix(Y1[,Y.to.use])
#' colnames(Y.new) = colnames(Y1)
#' colnames(Y.new)
#' SPLS.biplot_no_ax.labels(X.new, Y.new, 
#' algorithm=mod.SPLS, lambdaY=lambdaY.to.use, 
#' lambdaX=lambdaX.to.use, 
#' eps=1e-5, ax.tickvec.X=rep(1,ncol(X.new)), 
#' ax.tickvec.Y=rep(5,ncol(Y.new)))
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export
SPLS.biplot_no_ax.labels = function(X, Y, algorithm=NULL, eps, lambdaY=NULL, lambdaX=NULL, ax.tickvec.X=NULL, ax.tickvec.Y=NULL, ... ) 
{    
  options(digits=3)
  D = cbind(X,Y)
  X.mean = apply(X, 2, mean)  #column means of X or colMeans(X)
  Y.mean = apply(Y, 2, mean)  #column means of Y or colMeans(Y)
  X.std = apply(X, 2, sd)  #column standard deviations of X 
  Y.std = apply(Y, 2, sd)  #column standard deviations of Y   
  X.scal = scale(X, center=TRUE, scale=TRUE)  #scaled X
  Y.scal = scale(Y, center=TRUE, scale=TRUE)  #scaled Y  	
  
  #biplot	
  A = 2
  main = algorithm(X.scal,Y.scal,A,lambdaY,lambdaX,eps)
  Tmat = main$X.scores
  Pmat = main$X.loadings
  Qmat = main$Y.loadings
  Rmat = main$X.weights.trans		
  Bmat = Rmat %*% t(Qmat)  #(PxM) estimated SPLS coefficients matrix	 
  dimnames(Bmat) = list(colnames(X), colnames(Y))	    
  
  #overall quality of approximation		 
  X.hat = ((Tmat %*% t(Pmat)) * rep(X.std,each=nrow(X))) + rep(X.mean,each=nrow(X))
  dimnames(X.hat) = dimnames(X)
  Y.hat = ((Tmat %*% t(Qmat)) * rep(Y.std,each=nrow(Y))) + rep(Y.mean,each=nrow(Y))
  dimnames(Y.hat) = dimnames(Y)
  D.hat = cbind(X.hat, Y.hat)	   
  overall.quality = sum(diag(t(D.hat)%*%D.hat))  /  sum(diag(t(D)%*%D))
  
  #axes predictivity
  axis.pred = (diag(t(D.hat)%*%D.hat)) / (diag(t(D)%*%D))	    
  names(axis.pred) = colnames(D)		 
  
  #biplot points
  Z = rbind(Tmat, Rmat)  #((N+P)x2) biplot points
  dimnames(Z) = list(c(rownames(X),paste("b",1:ncol(X),sep="")), NULL)
  sample.style = biplot.sample.control(2, label=TRUE)
  sample.style$col = c("red", "purple")
  sample.style$pch = c(15,16)
  Gmat = cbind(c(rep(1,nrow(X)),rep(0,ncol(X))), c(rep(0,nrow(X)),rep(1,ncol(X))))
  dimnames(Gmat) = list(NULL, c("samples","coefficients"))
  classes = 1:2	 
  
  #axes direction    
  #X.hat=TP'
  X.axes.direction = (1 / (diag(Pmat%*%t(Pmat)))) * Pmat       	   
  #Y.hat=TQ'
  Y.axes.direction = (1 / (diag(Qmat%*%t(Qmat)))) * Qmat
  #B.spls=RQ'
  B.axes.direction = (1 / (diag(Qmat%*%t(Qmat)))) * Qmat	
  
  #calibration of axes	 	
  z.axes.X = lapply(1:ncol(X), calibrate.axis, X, X.mean, X.std, X.axes.direction, 
                    1:ncol(X), ax.tickvec.X, rep(0,ncol(X)), rep(0,ncol(X)), NULL)  
  z.axes.Y = lapply(1:ncol(Y), calibrate.axis, Y, Y.mean, Y.std, Y.axes.direction, 
                    1:ncol(Y), ax.tickvec.Y, rep(0,ncol(Y)), rep(0,ncol(Y)), NULL)
  z.axes.B = lapply(1:ncol(Bmat), calibrate.axis, Y.scal, rep(0,ncol(Y)), rep(1,ncol(Y)),
                    B.axes.direction, 1:ncol(Bmat), rep(4,ncol(Bmat)), rep(0,ncol(Bmat)), 
                    rep(0,ncol(Bmat)), NULL)	
  z.axes = vector("list", ncol(X)+ncol(Y)+ncol(Bmat))
  for (i in 1:ncol(X)) 
  {
    z.axes[[i]] = z.axes.X[[i]]
  }
  for (i in 1:ncol(Y)) 
  {
    z.axes[[i+ncol(X)]] = z.axes.Y[[i]]
  }
  for (i in 1:ncol(Bmat)) 
  {
    z.axes[[i+ncol(X)+ncol(Y)]] = z.axes.B[[i]]  
  }			  
  
  #biplot axes 
  ax.style = biplot.ax.control(ncol(X)+ncol(Y)+ncol(Bmat), c(colnames(X), colnames(Y),
                                                             colnames(Bmat)))
  ax.style$tick.col[1:ncol(X)] = "azure"
  ax.style$tick.label.col[1:ncol(X)] = "azure"
  ax.style$tick.col[ncol(X)+1:ncol(Y)] = "azure"
  ax.style$tick.label.col[ncol(X)+1:ncol(Y)] = "azure"
  ax.style$tick.col[ncol(X)+ncol(Y)+1:ncol(Bmat)] = "azure"
  ax.style$tick.label.col[ncol(X)+ncol(Y)+1:ncol(Bmat)] = "azure"
  #for X	
  ax.style$col[1:ncol(X)] = "blue"
  ax.style$label.col[1:ncol(X)] = "blue"
  #for Y	 
  ax.style$col[ncol(X)+1:ncol(Y)] = "black"
  ax.style$label.col[ncol(X)+1:ncol(Y)] = "black"
  #for B.spls
  ax.style$col[ncol(X)+ncol(Y)+1:ncol(Bmat)] = "black"
  ax.style$label.col[ncol(X)+ncol(Y)+1:ncol(Bmat)] = "black" 
  draw.biplot(Z, G=Gmat, classes=classes, sample.style=sample.style, z.axes=z.axes, 
              ax.style=ax.style, VV=rbind(X.axes.direction,Y.axes.direction))
  list(overall.quality=overall.quality, axis.pred=axis.pred, D.hat=D.hat, Bmat=Bmat)  	 
}

#' A zoomed-in display of the coefficient points in the Partial Least Squares (PLS) biplot for Sparse Partial Least Squares (SPLS)  
#' 
#' Takes in a set of predictor variables and a set of response variables and produces a zoomed-in display of the coefficient points in the PLS biplot for SPLS.
#' @param X A (NxP) predictor matrix 
#' @param Y A (NxM) response matrix
#' @param algorithm The SPLS algorithm 
#' @param eps Cut off value for convergence step 
#' @param lambdaY A value for the penalty parameters for the soft-thresholding penalization function for Y-weights 
#' @param lambdaX A value for the penalty parameters for the soft-thresholding penalization function for X-weights 
#' @param ax.tickvec.B (purple) tick marker length for the Y-variable axes in the biplot
#' @param ... Other arguments. Currently ignored
#' @return A zoomed-in display of the coefficient points in the PLS biplot of a SPLS-GLM of D=[X Y] with some parameters
#' @examples  
#' if(require(robustbase))
#' data(toxicity, package="robustbase")
#' Y1 = as.matrix(cbind(toxicity$toxicity))  
#' dimnames(Y1) = list(paste(1:nrow(Y1)), "toxicity")
#' X1 = as.matrix(cbind(toxicity[,2:10]))  
#' rownames(X1) = paste(1:nrow(X1))
#' #choosing a value for the penalty parameters lambdaY and lambdaX for this data 
#' main2 = opt.penalty.values(X=scale(X1), Y=scale(Y1), A=2, algorithm=mod.SPLS, eps=1e-5, 
#' from.value.X=0, to.value.X=10, from.value.Y=0, to.value.Y=0, lambdaY.len=1, lambdaX.len=100)
#' min.RMSEP.value = main2$min.RMSEP.value
#' lambdaY.to.use = main2$lambdaY.to.use
#' lambdaX.to.use = main2$lambdaX.to.use
#' list(lambdaY.to.use=lambdaY.to.use, lambdaX.to.use=lambdaX.to.use, min.RMSEP.value=min.RMSEP.value)
#' #SPLS analysis
#' main3 = mod.SPLS(X=scale(X1), Y=scale(Y1), A=2, 
#' lambdaY=lambdaY.to.use, lambdaX=lambdaX.to.use, 
#' eps=1e-5)
#' X.to.use = main3$X.select
#' Y.to.use = main3$Y.select
#' X.new = as.matrix(X1[,X.to.use])
#' Y.new = as.matrix(Y1[,Y.to.use])
#' colnames(Y.new) = colnames(Y1)
#' SPLS.biplot_Bmat(X.new, Y.new, algorithm=mod.SPLS, 
#' lambdaY=lambdaY.to.use, lambdaX=lambdaX.to.use, 
#' eps=1e-5, ax.tickvec.B=rep(5,ncol(Y.new)))
#'
#' #ash data  
#' if(require(chemometrics))
#' data(ash, package="chemometrics")
#' X1 = as.matrix(ash[,10:17], ncol=8)  
#' Y1 = as.matrix(ash$SOT)   
#' colnames(Y1) = paste("SOT")
#' #choosing a value for the penalty parameters lambdaY and lambdaX for this data 
#' main2 = opt.penalty.values(X=scale(X1), Y=scale(Y1), A=2, algorithm=mod.SPLS, eps=1e-5, 
#' from.value.X=0, to.value.X=500, from.value.Y=0, to.value.Y=0, lambdaY.len=1, lambdaX.len=100)
#' min.RMSEP.value = main2$min.RMSEP.value
#' lambdaY.to.use = main2$lambdaY.to.use
#' lambdaX.to.use = main2$lambdaX.to.use
#' list(lambdaY.to.use=lambdaY.to.use, lambdaX.to.use=lambdaX.to.use, min.RMSEP.value=min.RMSEP.value) 
#' #SPLS analysis
#' main3 = mod.SPLS(X=scale(X1), Y=scale(Y1), A=2, 
#' lambdaY=lambdaY.to.use, lambdaX=lambdaX.to.use, 
#' eps=1e-5)
#' X.to.use = main3$X.select
#' Y.to.use = main3$Y.select
#' X.new = as.matrix(X1[,X.to.use])
#' colnames(X.new)  #P=6
#' colnames(X1)  #P=8
#' Y.new = as.matrix(Y1[,Y.to.use])
#' colnames(Y.new) = colnames(Y1)
#' colnames(Y.new)
#' SPLS.biplot_Bmat(X.new, Y.new, algorithm=mod.SPLS,
#' lambdaY=lambdaY.to.use, lambdaX=lambdaX.to.use, 
#' eps=1e-5, ax.tickvec.B=5)
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export
SPLS.biplot_Bmat = function(X, Y, algorithm=NULL, eps, lambdaY=NULL, lambdaX=NULL, ax.tickvec.B=NULL, ... ) 
{    
  options(digits=3)      
  X.scal = scale(X, center=TRUE, scale=TRUE)  #scaled X
  Y.scal = scale(Y, center=TRUE, scale=TRUE)  #scaled Y  	
  
  #biplot	
  A = 2
  main = algorithm(X.scal,Y.scal,A,lambdaY,lambdaX,eps)
  Qmat = main$Y.loadings
  Rmat = main$X.weights.trans	
  Bmat = Rmat %*% t(Qmat)  #(PxM) estimated SPLS coefficients matrix	 
  dimnames(Bmat) = list(colnames(X), colnames(Y))	    
  
  #biplot points
  Z = Rmat  #
  dimnames(Z) = list(paste("b",1:ncol(X),sep=""), NULL)
  sample.style = biplot.sample.control(1, label=TRUE)
  sample.style$col = "purple"
  sample.style$pch = 16
  Gmat = cbind(rep(1,ncol(X)))
  dimnames(Gmat) = list(NULL, c("coefficients"))
  classes = 1	 
  
  #axes direction   
  #B.spls=RQ'
  B.axes.direction = (1 / (diag(Qmat%*%t(Qmat)))) * Qmat	
  
  #calibration of axes	
  z.axes.B = lapply(1:ncol(Bmat), calibrate.axis, Y.scal, rep(0,ncol(Y)), rep(1,ncol(Y)),
                    B.axes.direction, 1:ncol(Bmat), ax.tickvec.B, rep(0,ncol(Bmat)), 
                    rep(0,ncol(Bmat)), NULL)	
  z.axes = vector("list", ncol(Bmat))
  for (i in 1:ncol(Bmat)) 
  {
    z.axes[[i]] = z.axes.B[[i]]  
  }			  
  
  #biplot axes 
  ax.style = biplot.ax.control(ncol(Bmat), c(colnames(Bmat)))
  
  #for B.spls
  ax.style$col[1:ncol(Bmat)] = "black"
  ax.style$label.col[1:ncol(Bmat)] = "black" 
  ax.style$tick.col[1:ncol(Bmat)] = "purple"
  ax.style$tick.label.col[1:ncol(Bmat)] = "purple"   
  draw.biplot(Z, G=Gmat, classes=classes, sample.style=sample.style, z.axes=z.axes, 
              ax.style=ax.style, VV=B.axes.direction)
  list(Bmat=Bmat)  	 
}

#PLS biplot for SPLS-GLMs
#' The Partial Least Squares (PLS) biplot for Sparse Partial Least Squares-Generalized Linear Model (SPLS-GLM) 
#' 
#' Takes in a set of predictor variables and a set of response variables and produces a PLS biplot for the (univariate) SPLS-GLMs.
#' @param X A (NxP) predictor matrix 
#' @param y A (Nx1) response vector 
#' @param algorithm Any of the SPLS-GLM algorithm ("SPLS.GLM", "SPLS.binomial.GLM") 
#' @param eps Cut off value for convergence step 
#' @param lambdaY A value for the penalty parameters for the soft-thresholding penalization function for Y-weights 
#' @param lambdaX A value for the penalty parameters for the soft-thresholding penalization function for X-weights 
#' @param ax.tickvec.X tick marker length for each X-variable axis in the biplot
#' @param ax.tickvec.y tick marker length for the y-variable axis in the biplot
#' @param ax.tickvec.b (purple) tick marker length for the y-variable axis in the biplot
#' @param ... Other arguments. Currently ignored
#' @return The PLS biplot of a SPLS-GLM of D=[X y] with some parameters
#' @examples
#' if(require(robustbase))
#' possum.mat 
#' y = as.matrix(possum.mat[,1], ncol=1) 
#' dimnames(y) = list(paste("S", 1:nrow(possum.mat), seq=""), "Diversity")
#' X = as.matrix(possum.mat[,2:14], ncol=13) 
#' dimnames(X) = list(paste("S", 1:nrow(possum.mat), seq=""), colnames(possum.mat[,2:14])) 
#' #choosing a value for the penalty parameters lambdaY and lambdaX for this data 
#' main2B = opt.penalty.values(X=scale(X), Y=scale(y), A=2, algorithm=SPLS.GLM, eps=1e-3, 
#' from.value.X=1, to.value.X=4, from.value.Y=0, to.value.Y=0, lambdaY.len=1, lambdaX.len=100)
#' min.RMSEP.value = main2B$min.RMSEP.value
#' lambdaY.to.use = main2B$lambdaY.to.use
#' lambdaX.to.use = main2B$lambdaX.to.use
#' list(lambdaY.to.use=lambdaY.to.use, lambdaX.to.use=lambdaX.to.use, min.RMSEP.value=min.RMSEP.value)
#' #SPLS-GLM analysis
#' main3 = SPLS.GLM(scale(X), scale(y), A=2, lambdaY=lambdaY.to.use, lambdaX=lambdaX.to.use, eps=1e-3)
#' X.to.use = main3$X.select
#' X.new = as.matrix(X[,names(X.to.use)])
#' colnames(X.new)
#' main3$Y.select #note
#' SPLS.GLM.biplot(X.new, y, algorithm=SPLS.GLM, eps=1e-3, lambdaY=lambdaY.to.use, 
#' lambdaX=lambdaX.to.use, ax.tickvec.X=c(10,5,5,5,5,5,5,5,5,5,5,5,5), ax.tickvec.y=8, 
#' ax.tickvec.b=12)
#'
#' #Pima.tr data
#' if(require(MASS))
#' data(Pima.tr, package="MASS")
#' X = as.matrix(cbind(Pima.tr[,1:7]))  
#' dimnames(X) = list(1:nrow(X), colnames(X))
#' y = as.matrix(as.numeric(Pima.tr$type)-1, ncol=1)
#' #0=No and 1=Yes
#' dimnames(y) = list(1:nrow(y), paste("type"))
#' main2 = opt.penalty.values(X=scale(X), Y=scale(y), A=2, algorithm=SPLS.binomial.GLM, 
#' eps=1e-3, from.value.X=0, to.value.X=95, from.value.Y=0, to.value.Y=0, lambdaY.len=1, 
#' lambdaX.len=100)
#' min.RMSEP.value = main2$min.RMSEP.value
#' lambdaY.to.use = main2$lambdaY.to.use
#' lambdaX.to.use = main2$lambdaX.to.use
#' #SPLS-GLM analysis 
#' main3 = SPLS.binomial.GLM(scale(X), scale(y), A=2, lambdaY=lambdaY.to.use, 
#' lambdaX=lambdaX.to.use, eps=1e-3)
#' X.to.use = main3$X.select
#' X.new = as.matrix(X[,names(X.to.use)])
#' colnames(X.new) 
#' SPLS.GLM.biplot(X.new, y, algorithm=SPLS.binomial.GLM, eps=1e-3, lambdaY=lambdaY.to.use, 
#' lambdaX=lambdaX.to.use, ax.tickvec.X=c(3,3,3,3,3,3,1), ax.tickvec.y=3, ax.tickvec.b=3)
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export 
SPLS.GLM.biplot = function(X, y, algorithm=NULL, eps, lambdaY=NULL, lambdaX=NULL, ax.tickvec.X=NULL, ax.tickvec.y=NULL, ax.tickvec.b=NULL, ... ) 
{    
  options(digits=3)
  D = cbind(X,y)
  N = nrow(D)
  X.mean = apply(X, 2, mean)  #column means of X or colMeans(X)
  y.mean = mean(y)  #mean of y 
  X.std = apply(X, 2, sd)  #column standard deviations of X 
  y.std = sd(y)  #standard deviation of y   
  X.scal = scale(X, center=TRUE, scale=TRUE)  #scaled X
  y.scal = scale(y, center=TRUE, scale=TRUE)  #scaled y  
  
  #biplot 
  A = 2
  main = algorithm(X.scal,y.scal,A,lambdaY,lambdaX,eps)
  Tmat = main$X.scores
  Pmat = main$X.loadings
  qvec = main$y.loadings
  Rmat = main$X.weights.trans
  dimnames(Rmat) = list(colnames(X), paste("Comp", 1:A, seq=""))	 
  bvec = Rmat %*% t(qvec)  #(Px1) estimated SPLS-GLM coefficient vector 	 
  dimnames(bvec) = list(colnames(X), colnames(y))	
  
  #approximation		 
  X.hat = ( (Tmat %*% t(Pmat)) * rep(X.std,each=N) ) + rep(X.mean,each=N)
  dimnames(X.hat) = dimnames(X)
  y.hat = ( (Tmat %*% t(qvec)) * rep(y.std,each=N) ) + rep(y.mean,each=N)
  dimnames(y.hat) = list(rownames(y), paste("Expected_",colnames(y),sep=""))
  D.hat = cbind(X.hat, y.hat)   	    
  
  #biplot points
  Z = rbind(Tmat, Rmat)  #((N+P)x2) biplot points
  dimnames(Z) = list(c(rownames(X),paste("b",1:ncol(X),sep="")), NULL)
  sample.style = biplot.sample.control(2, label=TRUE)
  sample.style$col = c("red", "purple")
  sample.style$pch = c(15,16)
  Gmat = cbind(c(rep(1,nrow(X)),rep(0,ncol(X))), c(rep(0,nrow(X)),rep(1,ncol(X))))
  dimnames(Gmat) = list(NULL, c("samples","coefficients"))
  classes = 1:2	 
  
  #axes direction    
  #X.hat=TP'
  X.axes.direction = (1 / (diag(Pmat%*%t(Pmat)))) * Pmat       	   
  #y.hat=Tq
  qvec = matrix(qvec, nrow=1)  #since for M=1, j = 1:1
  y.axes.direction = (1 / (diag(qvec%*%t(qvec)))) * qvec
  #b.spls-glm=RQ'
  b.axes.direction = (1 / (diag(qvec%*%t(qvec)))) * qvec	
  
  #calibration of axes	 	
  z.axes.X = lapply(1:ncol(X), calibrate.axis, X, X.mean, X.std, X.axes.direction, 
                    1:ncol(X), ax.tickvec.X, rep(0,ncol(X)), rep(0,ncol(X)), NULL)  
  z.axes.y = lapply(1, calibrate.axis, y, y.mean, y.std, y.axes.direction, 
                    1, ax.tickvec.y, rep(0,1), rep(0,1), NULL)
  z.axes.b = lapply(1, calibrate.axis, y.scal, rep(0,1), rep(1,1), b.axes.direction, 
                    1, ax.tickvec.b, rep(0,1), rep(0,1), NULL)		
  z.axes = vector("list", ncol(X)+ncol(y)+ncol(bvec))
  for (i in 1:ncol(X)) 
  {
    z.axes[[i]] = z.axes.X[[i]]
  }
  for (i in 1:ncol(y)) 
  {
    z.axes[[i+ncol(X)]] = z.axes.y[[i]]
  }
  for (i in 1:ncol(bvec)) 
  {
    z.axes[[i+ncol(X)+ncol(y)]] = z.axes.b[[i]]  
  }			  
  
  #biplot axes 
  ax.style = biplot.ax.control(ncol(X)+ncol(y)+ncol(bvec), c(colnames(X), colnames(y), colnames(bvec)))
  #for X	 
  ax.style$col[1:ncol(X)] = "blue"
  ax.style$label.col[1:ncol(X)] = "blue"
  ax.style$tick.col[1:ncol(X)] = "blue"
  ax.style$tick.label.col[1:ncol(X)] = "blue"  
  #for y	 
  ax.style$col[ncol(X)+1:ncol(y)] = "black"
  ax.style$label.col[ncol(X)+1:ncol(y)] = "black"
  ax.style$tick.col[ncol(X)+1:ncol(y)] = "black"
  ax.style$tick.label.col[ncol(X)+1:ncol(y)] = "black"  
  #for b.spls-glm
  ax.style$col[ncol(X)+ncol(y)+1:ncol(bvec)] = "black"
  ax.style$label.col[ncol(X)+ncol(y)+1:ncol(bvec)] = "black" 
  ax.style$tick.col[ncol(X)+ncol(y)+1:ncol(bvec)] = "purple"
  ax.style$tick.label.col[ncol(X)+ncol(y)+1:ncol(bvec)] = "purple"   
  draw.biplot(Z, G=Gmat, classes=classes, sample.style=sample.style, z.axes=z.axes, ax.style=ax.style, VV=rbind(X.axes.direction,y.axes.direction)) 
  
  list(D.hat=D.hat, bvec=bvec)  	 
}

#' The Partial Least Squares (PLS) biplot for Sparse Partial Least Squares-Generalized Linear Model (SPLS-GLM), with the labels of the sample points excluded
#' 
#' Takes in a set of predictor variables and a set of response variable and produces a PLS biplot for the SPLS-GLM with the labels of the sample points excluded.
#' @param X A (NxP) predictor matrix 
#' @param y A (Nx1) response vector 
#' @param algorithm Any of the SPLS-GLM algorithm ("SPLS.GLM", "SPLS.binomial.GLM") 
#' @param eps Cut off value for convergence step 
#' @param lambdaY A value for the penalty parameters for the soft-thresholding penalization function for Y-weights 
#' @param lambdaX A value for the penalty parameters for the soft-thresholding penalization function for X-weights 
#' @param ax.tickvec.X tick marker length for each X-variable axis in the biplot
#' @param ax.tickvec.y tick marker length for the y-variable axis in the biplot
#' @param ax.tickvec.b (purple) tick marker length for the y-variable axis in the biplot
#' @param ... Other arguments. Currently ignored
#' @return The PLS biplot of a SPLS-GLM of D=[X y] with some parameters
#' @examples
#' if(require(robustbase))
#' possum.mat 
#' y = as.matrix(possum.mat[,1], ncol=1) 
#' dimnames(y) = list(paste("S", 1:nrow(possum.mat), seq=""), "Diversity")
#' X = as.matrix(possum.mat[,2:14], ncol=13) 
#' dimnames(X) = list(paste("S", 1:nrow(possum.mat), seq=""), colnames(possum.mat[,2:14])) 
#' #choosing a value for the penalty parameters lambdaY and lambdaX for this data 
#' main2B = opt.penalty.values(X=scale(X), Y=scale(y), A=2, algorithm=SPLS.GLM, 
#' eps=1e-3, from.value.X=1, to.value.X=4, from.value.Y=0, to.value.Y=0, lambdaY.len=1,
#' lambdaX.len=100)
#' min.RMSEP.value = main2B$min.RMSEP.value
#' lambdaY.to.use = main2B$lambdaY.to.use
#' lambdaX.to.use = main2B$lambdaX.to.use
#' list(lambdaY.to.use=lambdaY.to.use, lambdaX.to.use=lambdaX.to.use, min.RMSEP.value=min.RMSEP.value)
#' #SPLS-GLM analysis
#' main3 = SPLS.GLM(scale(X), scale(y), A=2, lambdaY=lambdaY.to.use, lambdaX=lambdaX.to.use, 
#' eps=1e-3)
#' X.to.use = main3$X.select
#' X.new = as.matrix(X[,names(X.to.use)])
#' colnames(X.new)
#' main3$Y.select #note
#' SPLS.GLM.biplot_no.SN(X.new, y, algorithm=SPLS.GLM, eps=1e-3, lambdaY=lambdaY.to.use, 
#' lambdaX=lambdaX.to.use, ax.tickvec.X=c(10,5,5,5,5,5,5,5,5,5,5,5,5), ax.tickvec.y=8, 
#' ax.tickvec.b=12) 
#'
#' #Pima.tr data
#' if(require(MASS))
#' data(Pima.tr, package="MASS")
#' X = as.matrix(cbind(Pima.tr[,1:7]))  
#' dimnames(X) = list(1:nrow(X), colnames(X))
#' y = as.matrix(as.numeric(Pima.tr$type)-1, ncol=1)
#' #0=No and 1=Yes
#' dimnames(y) = list(1:nrow(y), paste("type"))
#' main2 = opt.penalty.values(X=scale(X), Y=scale(y), A=2, algorithm=SPLS.binomial.GLM, 
#' eps=1e-3, from.value.X=0, to.value.X=95, from.value.Y=0, to.value.Y=0, lambdaY.len=1, 
#' lambdaX.len=100)
#' min.RMSEP.value = main2$min.RMSEP.value
#' lambdaY.to.use = main2$lambdaY.to.use
#' lambdaX.to.use = main2$lambdaX.to.use
#' #SPLS-GLM analysis 
#' main3 = SPLS.binomial.GLM(scale(X), scale(y), A=2, lambdaY=lambdaY.to.use, 
#' lambdaX=lambdaX.to.use, eps=1e-3)
#' X.to.use = main3$X.select
#' X.new = as.matrix(X[,names(X.to.use)])
#' colnames(X.new) 
#' SPLS.GLM.biplot_no.SN(X.new, y, algorithm=SPLS.binomial.GLM, eps=1e-3, 
#' lambdaY=lambdaY.to.use, lambdaX=lambdaX.to.use, 
#' ax.tickvec.X=c(3,3,3,3,3,3,1), ax.tickvec.y=3, ax.tickvec.b=3)
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export
SPLS.GLM.biplot_no.SN = function(X, y, algorithm=NULL, eps, lambdaY=NULL, lambdaX=NULL, ax.tickvec.X=NULL, ax.tickvec.y=NULL, ax.tickvec.b=NULL, ... ) 
{    
  options(digits=3)
  D = cbind(X,y)
  N = nrow(D)
  X.mean = apply(X, 2, mean)  #column means of X or colMeans(X)
  y.mean = mean(y)  #mean of y 
  X.std = apply(X, 2, sd)  #column standard deviations of X 
  y.std = sd(y)  #standard deviation of y   
  X.scal = scale(X, center=TRUE, scale=TRUE)  #scaled X
  y.scal = scale(y, center=TRUE, scale=TRUE)  #scaled y  
  
  #biplot 
  A = 2
  main = algorithm(X.scal,y.scal,A,lambdaY,lambdaX,eps)
  Tmat = main$X.scores
  Pmat = main$X.loadings
  qvec = main$y.loadings
  Rmat = main$X.weights.trans
  dimnames(Rmat) = list(colnames(X), paste("Comp", 1:A, seq=""))	 
  bvec = Rmat %*% t(qvec)  #(Px1) estimated SPLS-GLM coefficient vector 	 
  dimnames(bvec) = list(colnames(X), colnames(y))	
  
  #approximation		 
  X.hat = ( (Tmat %*% t(Pmat)) * rep(X.std,each=N) ) + rep(X.mean,each=N)
  dimnames(X.hat) = dimnames(X)
  y.hat = ( (Tmat %*% t(qvec)) * rep(y.std,each=N) ) + rep(y.mean,each=N)
  dimnames(y.hat) = list(rownames(y), paste("Expected_",colnames(y),sep=""))
  D.hat = cbind(X.hat, y.hat)   	    
  
  #biplot points
  Z = rbind(Tmat, Rmat)  #((N+P)x2) biplot points
  dimnames(Z) = list(c(rownames(X),paste("b",1:ncol(X),sep="")), NULL)
  sample.style = biplot.sample.control(2, label=c(FALSE,TRUE))
  sample.style$col = c("red", "purple")
  sample.style$pch = c(15,16)
  Gmat = cbind(c(rep(1,nrow(X)),rep(0,ncol(X))), c(rep(0,nrow(X)),rep(1,ncol(X))))
  dimnames(Gmat) = list(NULL, c("samples","coefficients"))
  classes = 1:2	 
  
  #axes direction    
  #X.hat=TP'
  X.axes.direction = (1 / (diag(Pmat%*%t(Pmat)))) * Pmat       	   
  #y.hat=Tq
  qvec = matrix(qvec, nrow=1)  #since for M=1, j = 1:1
  y.axes.direction = (1 / (diag(qvec%*%t(qvec)))) * qvec
  #b.spls-glm=RQ'
  b.axes.direction = (1 / (diag(qvec%*%t(qvec)))) * qvec	
  
  #calibration of axes	 	
  z.axes.X = lapply(1:ncol(X), calibrate.axis, X, X.mean, X.std, X.axes.direction, 
                    1:ncol(X), ax.tickvec.X, rep(0,ncol(X)), rep(0,ncol(X)), NULL)  
  z.axes.y = lapply(1, calibrate.axis, y, y.mean, y.std, y.axes.direction, 
                    1, ax.tickvec.y, rep(0,1), rep(0,1), NULL)
  z.axes.b = lapply(1, calibrate.axis, y.scal, rep(0,1), rep(1,1), b.axes.direction, 
                    1, ax.tickvec.b, rep(0,1), rep(0,1), NULL)		
  z.axes = vector("list", ncol(X)+ncol(y)+ncol(bvec))
  for (i in 1:ncol(X)) 
  {
    z.axes[[i]] = z.axes.X[[i]]
  }
  for (i in 1:ncol(y)) 
  {
    z.axes[[i+ncol(X)]] = z.axes.y[[i]]
  }
  for (i in 1:ncol(bvec)) 
  {
    z.axes[[i+ncol(X)+ncol(y)]] = z.axes.b[[i]]  
  }			  
  
  #biplot axes 
  ax.style = biplot.ax.control(ncol(X)+ncol(y)+ncol(bvec), c(colnames(X), colnames(y), colnames(bvec)))
  #for X	 
  ax.style$col[1:ncol(X)] = "blue"
  ax.style$label.col[1:ncol(X)] = "blue"
  ax.style$tick.col[1:ncol(X)] = "blue"
  ax.style$tick.label.col[1:ncol(X)] = "blue"  
  #for y	 
  ax.style$col[ncol(X)+1:ncol(y)] = "black"
  ax.style$label.col[ncol(X)+1:ncol(y)] = "black"
  ax.style$tick.col[ncol(X)+1:ncol(y)] = "black"
  ax.style$tick.label.col[ncol(X)+1:ncol(y)] = "black"  
  #for b.spls-glm
  ax.style$col[ncol(X)+ncol(y)+1:ncol(bvec)] = "black"
  ax.style$label.col[ncol(X)+ncol(y)+1:ncol(bvec)] = "black" 
  ax.style$tick.col[ncol(X)+ncol(y)+1:ncol(bvec)] = "purple"
  ax.style$tick.label.col[ncol(X)+ncol(y)+1:ncol(bvec)] = "purple"   
  draw.biplot(Z, G=Gmat, classes=classes, sample.style=sample.style, z.axes=z.axes, ax.style=ax.style, VV=rbind(X.axes.direction,y.axes.direction)) 
  
  list(D.hat=D.hat, bvec=bvec)  	 
}

#' A zoomed-in display of the coefficient points in the Partial Least Squares (PLS) biplot for Sparse Partial Least Squares-Generalized Linear Model (SPLS-GLM)  
#' 
#' Takes in a set of predictor variables and a set of response variables and produces a zoomed-in display of the coefficient points in the PLS biplot for the (univariate) SPLS-GLMs.
#' @param X A (NxP) predictor matrix 
#' @param y A (Nx1) response vector 
#' @param algorithm Any of the SPLS-GLM algorithm ("SPLS.GLM", "SPLS.binomial.GLM") 
#' @param eps Cut off value for convergence step 
#' @param lambdaY A value for the penalty parameters for the soft-thresholding penalization function for Y-weights 
#' @param lambdaX A value for the penalty parameters for the soft-thresholding penalization function for X-weights 
#' @param ax.tickvec.b (purple) tick marker length for the y-variable axis in the biplot
#' @param ... Other arguments. Currently ignored
#' @return A zoomed-in display of the coefficient points in the PLS biplot of a SPLS-GLM of D=[X y] with some parameters
#' @examples  
#' if(require(robustbase))
#' possum.mat 
#' y = as.matrix(possum.mat[,1], ncol=1) 
#' dimnames(y) = list(paste("S", 1:nrow(possum.mat), seq=""), "Diversity")
#' X = as.matrix(possum.mat[,2:14], ncol=13) 
#' dimnames(X) = list(paste("S", 1:nrow(possum.mat), seq=""), colnames(possum.mat[,2:14])) 
#' #choosing a value for the penalty parameters lambdaY and lambdaX for this data 
#' main2B = opt.penalty.values(X=scale(X), Y=scale(y), A=2, algorithm=SPLS.GLM, eps=1e-3, 
#' from.value.X=1, to.value.X=4, from.value.Y=0, to.value.Y=0, lambdaY.len=1, lambdaX.len=100)
#' min.RMSEP.value = main2B$min.RMSEP.value
#' lambdaY.to.use = main2B$lambdaY.to.use
#' lambdaX.to.use = main2B$lambdaX.to.use
#' list(lambdaY.to.use=lambdaY.to.use, lambdaX.to.use=lambdaX.to.use, min.RMSEP.value=min.RMSEP.value)
#' #SPLS-GLM analysis
#' main3 = SPLS.GLM(scale(X), scale(y), A=2, lambdaY=lambdaY.to.use, lambdaX=lambdaX.to.use, eps=1e-3)
#' X.to.use = main3$X.select
#' X.new = as.matrix(X[,names(X.to.use)])
#' colnames(X.new)
#' main3$Y.select #note
#' SPLS.GLM.biplot_bvec(X.new, y, algorithm=SPLS.GLM, eps=1e-3, lambdaY=lambdaY.to.use,
#' lambdaX=lambdaX.to.use, ax.tickvec.b=30)   
#'
#' #Pima.tr data
#' if(require(MASS))
#' data(Pima.tr, package="MASS")
#' X = as.matrix(cbind(Pima.tr[,1:7]))  
#' dimnames(X) = list(1:nrow(X), colnames(X))
#' y = as.matrix(as.numeric(Pima.tr$type)-1, ncol=1)
#' #0=No and 1=Yes
#' dimnames(y) = list(1:nrow(y), paste("type"))
#' main2 = opt.penalty.values(X=scale(X), Y=scale(y), A=2, algorithm=SPLS.binomial.GLM, eps=1e-3, 
#' from.value.X=0, to.value.X=95, from.value.Y=0, to.value.Y=0, lambdaY.len=1, lambdaX.len=100)
#' min.RMSEP.value = main2$min.RMSEP.value
#' lambdaY.to.use = main2$lambdaY.to.use
#' lambdaX.to.use = main2$lambdaX.to.use
#' #SPLS-GLM analysis 
#' main3 = SPLS.binomial.GLM(scale(X), scale(y), A=2, lambdaY=lambdaY.to.use, lambdaX=lambdaX.to.use,
#' eps=1e-3)
#' X.to.use = main3$X.select
#' X.new = as.matrix(X[,names(X.to.use)])
#' colnames(X.new) 
#' SPLS.GLM.biplot_bvec(X.new, y, algorithm=SPLS.binomial.GLM, eps=1e-3, lambdaY=lambdaY.to.use,
#' lambdaX=lambdaX.to.use, ax.tickvec.b=30) 
#' @author Opeoluwa F. Oyedele and Sugnet Gardner-Lubbe 
#' @export
SPLS.GLM.biplot_bvec = function(X, y, algorithm=NULL, eps, lambdaY=NULL, lambdaX=NULL, ax.tickvec.b=NULL, ... ) 
{    
  options(digits=3)
  X.scal = scale(X, center=TRUE, scale=TRUE)  #scaled X
  y.scal = scale(y, center=TRUE, scale=TRUE)  #scaled y  
  
  #biplot 
  A = 2
  main = algorithm(X.scal,y.scal,A,lambdaY,lambdaX,eps)
  qvec = main$y.loadings
  Rmat = main$X.weights.trans
  dimnames(Rmat) = list(colnames(X), paste("Comp", 1:A, seq=""))	 
  bvec = Rmat %*% t(qvec)  #(Px1) estimated SPLS-GLM coefficient vector 	 
  dimnames(bvec) = list(colnames(X), colnames(y))	
  
  #biplot points
  Z = Rmat  #
  dimnames(Z) = list(paste("b",1:ncol(X),sep=""), NULL)
  sample.style = biplot.sample.control(1, label=TRUE)
  sample.style$col = "purple"
  sample.style$pch = 16
  Gmat = cbind(rep(1,ncol(X)))
  dimnames(Gmat) = list(NULL, c("coefficients"))
  classes = 1	
  
  #axes direction  
  qvec = matrix(qvec, nrow=1)  #since for M=1, j = 1:1
  #b.spls-glm=RQ'
  b.axes.direction = (1 / (diag(qvec%*%t(qvec)))) * qvec	
  
  #calibration of axes	
  z.axes.b = lapply(1, calibrate.axis, y.scal, rep(0,1), rep(1,1), b.axes.direction, 
                    1, ax.tickvec.b, rep(0,1), rep(0,1), NULL)	
  z.axes = vector("list", ncol(bvec))
  for (i in 1:ncol(bvec )) 
  {
    z.axes[[i]] = z.axes.b[[i]]  
  }			  
  
  #biplot axes 
  ax.style = biplot.ax.control(ncol(bvec), c(colnames(bvec)))
  
  #for B.spls
  ax.style$col[1:ncol(bvec)] = "black"
  ax.style$label.col[1:ncol(bvec)] = "black" 
  ax.style$tick.col[1:ncol(bvec)] = "purple"
  ax.style$tick.label.col[1:ncol(bvec)] = "purple"   
  draw.biplot(Z, G=Gmat, classes=classes, sample.style=sample.style, z.axes=z.axes, 
              ax.style=ax.style, VV=b.axes.direction)
  list(bvec=bvec)  	 	    	 
}

# --------------------------------------------------------------------------------------------------------------------
