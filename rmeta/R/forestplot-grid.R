

forestplot<-function(labeltext,mean,lower,upper,align=NULL, is.summary=FALSE, clip=c(-Inf,Inf), 
                     xlab="", zero= 0,graphwidth=unit(2,"inches"),col=meta.colors(),xlog=FALSE,
                     xticks=NULL, boxsize=NULL, ...){

  require("grid") || stop("`grid' package not found")
  require("rmeta") || stop("`rmeta' package not found")

  ## Function to draw a non-summary rect-plus-CI
  drawNormalCI <- function(LL, OR, UL, size) {
    size=0.75*size

    clipupper<-convertX(unit(UL, "native"), "npc", valueOnly=TRUE) > 1
    cliplower<-convertX(unit(LL, "native"), "npc", valueOnly=TRUE) <0
    box<- convertX(unit(OR, "native"), "npc", valueOnly=TRUE)
    clipbox <- box<0 || box>1
    
    
    ## Draw arrow if exceed col range
    ## convertX() used to convert between coordinate systems
    if (clipupper || cliplower){
      ends<-"both"
      lims<-unit(c(0, 1), c("npc", "npc"))
      if (!clipupper) {
        ends<-"first"
        lims<-unit(c(0, UL), c("npc","native"))
      }
      if (!cliplower) {
        ends<-"last"
        lims<-unit(c(LL, 1), c("native", "npc"))
      }
      grid.lines(x=lims, y=0.5,arrow=arrow(ends=ends,length=unit(0.05, "inches")),
                 gp=gpar(col=col$lines))

      if (!clipbox)
          grid.rect(x=unit(OR, "native"),
                    width=unit(size, "snpc"), height=unit(size, "snpc"),
                    gp=gpar(fill=col$box,col=col$box))
      
      } else   {
      ## Draw line white if totally inside rect
      grid.lines(x=unit(c(LL, UL), "native"), y=0.5,
                 gp=gpar(col=col$lines))
      grid.rect(x=unit(OR, "native"),
                width=unit(size, "snpc"), height=unit(size, "snpc"),
                gp=gpar(fill=col$box,col=col$box))
      if ((convertX(unit(OR, "native") + unit(0.5*size, "lines"), "native", valueOnly=TRUE) > UL) &&
          (convertX(unit(OR, "native") - unit(0.5*size, "lines"), "native", valueOnly=TRUE) < LL))
        grid.lines(x=unit(c(LL, UL), "native"), y=0.5, gp=gpar(col=col$lines))
    }
  }
  
  ## Function to draw a summary "diamond"
  drawSummaryCI <- function(LL, OR, UL, size) {
    grid.polygon(x=unit(c(LL, OR, UL, OR), "native"),
                 y=unit(0.5 + c(0, 0.5*size, 0, -0.5*size), "npc"),gp=gpar(fill=col$summary,col=col$summary))
  }
  
  plot.new()
  ## calculate width based on labels with something in every column
  widthcolumn<-!apply(is.na(labeltext),1,any)
  
  nc<-NCOL(labeltext)
  labels<-vector("list",nc)

  
  if (is.null(align))
    align<-c("l",rep("r",nc-1))
  else
    align<-rep(align,length=nc)
  
  nr<-NROW(labeltext)
  is.summary<-rep(is.summary,length=nr)
  
  for(j in 1:nc){
    labels[[j]]<-vector("list", nr)
    for(i in 1:nr){
      if (is.na(labeltext[i,j]))
        next
      x<-switch(align[j],l=0,r=1,c=0.5)
      just<-switch(align[j],l="left",r="right",c="center")
      labels[[j]][[i]]<-textGrob(labeltext[i,j], x=x,just=just,
                                 gp=gpar(fontface=if(is.summary[i]) "bold" else "plain",
                                                                    col=rep(col$text,length=nr)[i]) )
    }
  }
  
  colgap<-unit(3,"mm")
  
  colwidths<-unit.c(max(unit(rep(1,sum(widthcolumn)),"grobwidth",labels[[1]][widthcolumn])),colgap)
  if (nc>1){
    for(i in 2:nc)
      colwidths<-unit.c(colwidths, max(unit(rep(1,sum(widthcolumn)),"grobwidth",labels[[i]][widthcolumn])),colgap)
    
  }
  colwidths<-unit.c(colwidths,graphwidth)
  
  pushViewport(viewport(layout=grid.layout(nr+1,nc*2+1,
                          widths=colwidths,
                          heights=unit(c(rep(1, nr),0.5), "lines"))))
  
  
  
  cwidth<-(upper-lower)

  xrange<-c(max(min(lower,na.rm=TRUE),clip[1]), min(max(upper,na.rm=TRUE),clip[2]))
  
  info<-1/cwidth
  info<-info/max(info[!is.summary], na.rm=TRUE)
  info[is.summary]<-1

  if (!is.null(boxsize)) info<-rep(boxsize,length=length(info))
  
  for(j in 1:nc){
    for(i in 1:nr){
      if (!is.null(labels[[j]][[i]])){
        pushViewport(viewport(layout.pos.row=i,layout.pos.col=2*j-1))
        grid.draw(labels[[j]][[i]])
          popViewport()
      }
    }
  }
  
  pushViewport(viewport(layout.pos.col=2*nc+1, xscale=xrange))
  grid.lines(x=unit(zero, "native"), y=0:1,gp=gpar(col=col$zero))
  if (xlog){
      if(is.null(xticks)){
          ticks<-pretty(exp(xrange))
          ticks<-ticks[ticks>0]
      } else{
          ticks<-xticks
      }
      if (length(ticks)){ 
          if (min(lower,na.rm=TRUE)<clip[1]) ticks<-c(exp(clip[1]),ticks)
          if (max(upper,na.rm=TRUE)>clip[2]) ticks<-c(ticks,exp(clip[2]))
          xax<-xaxisGrob(gp=gpar(cex=0.6,col=col$axes),at=log(ticks),name="xax")
          xax1<-editGrob(xax, gPath("labels"), label=format(ticks,digits=2))
          grid.draw(xax1)
      }
  } else{
      if (is.null(xticks)){
          grid.xaxis(gp=gpar(cex=0.6,col=col$axes))
      } else if(length(xticks)) {
          grid.xaxis(at=xticks, gp=gpar(cex=0.6,col=col$axes))
      }
      
  }
  grid.text(xlab, y=unit(-2, "lines"),gp=gpar(col=col$axes))
  popViewport()
  for (i in 1:nr) {
    if (is.na(mean[i])) next
    pushViewport(viewport(layout.pos.row=i, layout.pos.col=2*nc+1,
                          xscale=xrange))
    if (is.summary[i])
      drawSummaryCI(lower[i], mean[i], upper[i], info[i])
    else
      drawNormalCI(lower[i], mean[i], upper[i], info[i])
    
    popViewport()
  }
  
  
  popViewport()
}





