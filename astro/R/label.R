label = function(pos = "topleft", lab = "label", txt = NULL, inset = 0.1, whitespace = 0.08, col = "black", bgcol = "white", bty = "n", bordercol = "black", lwd = 1, cex = 1, align = "center"){
    #pos = "topleft"; lab = "label"; txt = NULL; inset = 0.1; whitespace = 0.08; col = "black"; bgcol = "white"; bty = "n"; bordercol = "black"; lwd = 1; cex = 1; align = "center"
    # LABEL Function
    # adds a label to a plot in the plotting region
    if(length(txt)>0){
        lab=txt     # backwards compatability for txt, eventually remove
    }
    xlo2 = par("usr")[1]; xhi2 = par("usr")[2]; ylo2 = par("usr")[3]; yhi2 = par("usr")[4]
    #if(round((axTicks(1)[length(axTicks(1))]-axTicks(1)[length(axTicks(1))-1])-(axTicks(1)[length(axTicks(1))-1]-axTicks(1)[length(axTicks(1))-2]), digits=30)!=0){logx=TRUE}else{logx=FALSE}
    #if(round((axTicks(2)[length(axTicks(2))]-axTicks(2)[length(axTicks(2))-1])-(axTicks(2)[length(axTicks(2))-1]-axTicks(2)[length(axTicks(2))-2]), digits=30)!=0){logy=TRUE}else{logy=FALSE}
    xdim=par("pin")[1]; ydim=par("pin")[2]
    if(length(inset)==2){xlo = xlo2+(((xhi2-xlo2)/xdim)*inset[1]); xhi = xhi2-(((xhi2-xlo2)/xdim)*inset[1]); ylo = ylo2+(((yhi2-ylo2)/ydim)*inset[2]); yhi = yhi2-(((yhi2-ylo2)/ydim)*inset[2])}else{xlo = xlo2+(((xhi2-xlo2)/xdim)*inset); xhi = xhi2-(((xhi2-xlo2)/xdim)*inset); ylo = ylo2+(((yhi2-ylo2)/ydim)*inset); yhi = yhi2-(((yhi2-ylo2)/ydim)*inset)}
    xmid = ((xhi2-xlo2)/2)+xlo2; ymid = ((yhi2-ylo2)/2)+ylo2
    textwidth = strwidth(lab)*cex; textheight = strheight(lab)*cex
    rectwidth = textwidth+(((xhi2-xlo2)/xdim)*whitespace); rectheight = textheight+(((yhi2-ylo2)/ydim)*whitespace)
    if(pos=="topleft"){xleft=xlo; xright=xlo+rectwidth; ybottom=yhi-rectheight; ytop=yhi
    }else if(pos=="top"){xleft=xmid-rectwidth/2; xright=xmid+rectwidth/2; ybottom=yhi-rectheight; ytop=yhi
    }else if(pos=="topright"){xleft=xhi-rectwidth; xright=xhi; ybottom=yhi-rectheight; ytop=yhi
    }else if(pos=="right"){xleft=xhi-rectwidth; xright=xhi; ybottom=ymid-rectheight/2; ytop=ymid+rectheight/2
    }else if(pos=="bottomright"){xleft=xhi-rectwidth; xright=xhi; ybottom=ylo; ytop=ylo+rectheight
    }else if(pos=="bottom"){xleft=xmid-rectwidth/2; xright=xmid+rectwidth/2; ybottom=ylo; ytop=ylo+rectheight
    }else if(pos=="bottomleft"){xleft=xlo; xright=xlo+rectwidth; ybottom=ylo; ytop=ylo+rectheight
    }else if(pos=="left"){xleft=xlo; xright=xlo+rectwidth; ybottom=ymid-rectheight/2; ytop=ymid+rectheight/2
    }else if(pos=="centre"){xleft=xmid-rectwidth/2; xright=xmid+rectwidth/2; ybottom=ymid-rectheight/2; ytop=ymid+rectheight/2}
    if(align=="left"){tadj=c(0,0.5); xtext = xleft+(((xhi2-xlo2)/xdim)*whitespace)/2
    }else if(align=="right"){tadj=c(1,0.5); xtext = xleft+rectwidth-(((xhi2-xlo2)/xdim)*whitespace)/2
    }else{tadj=c(0.5,0.5); xtext = xleft+rectwidth/2}
    ytext=ybottom+rectheight/2
    if(par("xlog")){xleft=10^xleft; xright=10^xright; xtext=10^xtext}
    if(par("ylog")){ybottom=10^ybottom; ytop=10^ytop; ytext=10^ytext}
    if (bty!="n") {rect(xleft, ybottom, xright, ytop, col=bgcol, border=bordercol, lwd=lwd)}
    text(x=xtext, y=ytext, labels=lab, col=col, cex=cex, adj=tadj)
}

