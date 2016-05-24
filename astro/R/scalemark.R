scalemark = function(len = "AUTO", txt = "AUTO", pixsize = NA, col = "white", linecol = col, bg = "grey25", pos = "bottomleft", inset = 0.1, cex = 0.9, lwd = 1.2){
    # adds a scale mark to a plot in the plotting region
    xlo2 = par("usr")[1]; xhi2 = par("usr")[2]; ylo2 = par("usr")[3]; yhi2 = par("usr")[4]
    if (round((axTicks(1)[length(axTicks(1))]-axTicks(1)[length(axTicks(1))-1])-(axTicks(1)[length(axTicks(1))-1]-axTicks(1)[length(axTicks(1))-2]), digits=10)!=0) {logx=TRUE}else{logx=FALSE}
    if (round((axTicks(2)[length(axTicks(2))]-axTicks(2)[length(axTicks(2))-1])-(axTicks(2)[length(axTicks(2))-1]-axTicks(2)[length(axTicks(2))-2]), digits=10)!=0) {logy=TRUE}else{logy=FALSE}
    xdim=par("pin")[1]; ydim=par("pin")[2]
    if(length(inset)==2){xlo = xlo2+(((xhi2-xlo2)/xdim)*inset[1]); xhi = xhi2-(((xhi2-xlo2)/xdim)*inset[1]); ylo = ylo2+(((yhi2-ylo2)/ydim)*inset[2]); yhi = yhi2-(((yhi2-ylo2)/ydim)*inset[2])}else{xlo = xlo2+(((xhi2-xlo2)/xdim)*inset); xhi = xhi2-(((xhi2-xlo2)/xdim)*inset); ylo = ylo2+(((yhi2-ylo2)/ydim)*inset); yhi = yhi2-(((yhi2-ylo2)/ydim)*inset)}
    xmid = ((xhi2-xlo2)/2)+xlo2; ymid = ((yhi2-ylo2)/2)+ylo2
    if(len=="AUTO"){
        if(is.na(pixsize)){pixsize=1}
        possibles=c(0.1,0.5,1,5,10,50,100,500,1000)
        len=possibles[which.min(abs(((possibles/pixsize)*3)-(xhi2-xlo2)))]/pixsize
        txt=paste(len*pixsize,"\'\'",sep="")
    }
    rectwidth = len; rectheight = strheight("lee")*0.5
    if (pos=="topleft") {xleft=xlo; xright=xlo+rectwidth; ybottom=yhi-rectheight; ytop=yhi}
    else if (pos=="top") {xleft=xmid-rectwidth/2; xright=xmid+rectwidth/2; ybottom=yhi-rectheight; ytop=yhi}
    else if (pos=="topright") {xleft=xhi-rectwidth; xright=xhi; ybottom=yhi-rectheight; ytop=yhi} 
    else if (pos=="right") {xleft=xhi-rectwidth; xright=xhi; ybottom=ymid-rectheight/2; ytop=ymid+rectheight/2} 
    else if (pos=="bottomright") {xleft=xhi-rectwidth; xright=xhi; ybottom=ylo; ytop=ylo+rectheight} 
    else if (pos=="bottom") {xleft=xmid-rectwidth/2; xright=xmid+rectwidth/2; ybottom=ylo; ytop=ylo+rectheight} 
    else if (pos=="bottomleft") {xleft=xlo; xright=xlo+rectwidth; ybottom=ylo; ytop=ylo+rectheight} 
    else if (pos=="left") {xleft=xlo; xright=xlo+rectwidth; ybottom=ymid-rectheight/2; ytop=ymid+rectheight/2}
    xtext = xleft+((xhi2-xlo2)/xdim)*0.12; ytext=ybottom+((yhi2-ylo2)/ydim)*0
    if (logx) {xleft=10^xleft; xright=10^xright; xtext=10^xtext}
    if (logy) {ybottom=10^ybottom; ytop=10^ytop; ytext=10^ytext}
    lines(c(xleft,xleft),c(ytop,ybottom),lwd=3*lwd,col=bg)
    lines(c(xleft,xright),c(ybottom,ybottom),lwd=3*lwd,col=bg)
    lines(c(xright,xright),c(ybottom,ytop),lwd=3*lwd,col=bg)
    lines(c(xleft,xleft),c(ytop,ybottom),lwd=lwd,col=linecol)
    lines(c(xleft,xright),c(ybottom,ybottom),lwd=lwd,col=linecol)
    lines(c(xright,xright),c(ybottom,ytop),lwd=lwd,col=linecol)
    shadowtext(x=xtext, y=ytext, labels=txt, col=col, cex=cex, pos=3)
}

