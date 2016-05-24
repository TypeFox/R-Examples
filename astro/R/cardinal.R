cardinal = function(rotate = 0, pos = "bottomright", inset = 0.5, len = 2, gap = 1, col = "white", linecol = col, bg = "grey25", lwd = 1, invert = TRUE, south = TRUE, west = TRUE, textrot = TRUE, ...){
    
    # origins
    xlo2 = par("usr")[1]; xhi2 = par("usr")[2]; ylo2 = par("usr")[3]; yhi2 = par("usr")[4]
    xdim=par("pin")[1]; ydim=par("pin")[2]
    if(length(inset)==2){xlo = xlo2+(((xhi2-xlo2)/xdim)*inset[1]); xhi = xhi2-(((xhi2-xlo2)/xdim)*inset[1]); ylo = ylo2+(((yhi2-ylo2)/ydim)*inset[2]); yhi = yhi2-(((yhi2-ylo2)/ydim)*inset[2])}else{xlo = xlo2+(((xhi2-xlo2)/xdim)*inset); xhi = xhi2-(((xhi2-xlo2)/xdim)*inset); ylo = ylo2+(((yhi2-ylo2)/ydim)*inset); yhi = yhi2-(((yhi2-ylo2)/ydim)*inset)}
    xmid = ((xhi2-xlo2)/2)+xlo2; ymid = ((yhi2-ylo2)/2)+ylo2
    if(pos=="bottomleft"){x = xlo; y = ylo}
    if(pos=="topleft"){x = xlo; y = yhi}
    if(pos=="topright"){x = xhi; y = yhi}
    if(pos=="bottomright"){x = xhi; y = ylo}
    
    # length
    cd11 = cos(rotate*(pi/180))
    cd12 = -sin(rotate*(pi/180))
    cd21 = sin(rotate*(pi/180))
    cd22 = cos(rotate*(pi/180))
    if(invert){
        cd11 = -cd11
        cd12 = -cd12
    }
    rule = strheight("lee")
    scale = len*rule
    tscale = (len*rule) + (gap*rule)
    nscale = scale/sqrt((cd21^2) + (cd22^2))
    escale = scale/sqrt((cd11^2) + (cd12^2))
    ntscale = tscale/sqrt((cd21^2) + (cd22^2))
    etscale = tscale/sqrt((cd11^2) + (cd12^2))
    
    # north-south
    xn = x + (nscale*cd21)
    yn = y + (nscale*cd22)
    xs = x - (nscale*cd21)
    ys = y - (nscale*cd22)
    xnt = x + (ntscale*cd21)
    ynt = y + (ntscale*cd22)
    xst = x - (ntscale*cd21)
    yst = y - (ntscale*cd22)
    
    # east-west
    xe = x + (escale*cd11)
    ye = y + (escale*cd12)
    xw = x - (escale*cd11)
    yw = y - (escale*cd12)
    xet = x + (etscale*cd11)
    yet = y + (etscale*cd12)
    xwt = x - (etscale*cd11)
    ywt = y - (etscale*cd12)
    
    # arrows
    if(!south){
        xs = x
        ys = y
    }
    if(!west){
        xw = x
        yw = y
    }
    arrows(xs,ys,xn,yn, length=0.05, col=bg, lwd=3*lwd)
    arrows(xw,yw,xe,ye, length=0.05, col=bg, lwd=3*lwd)
    arrows(xs,ys,xn,yn, length=0.05, col=linecol, lwd=lwd)
    arrows(xw,yw,xe,ye, length=0.05, col=linecol, lwd=lwd)
    if(!textrot){
        rotate = 0
    }
    shadowtext(xnt,ynt,labels="N",col=col,srt=-rotate,...)
    shadowtext(xet,yet,labels="E",col=col,srt=-rotate,...)
    
}

