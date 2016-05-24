
#plots a box...

pbox <- function(sdobj,xdim,ydim,boxnum,fromtype="oldbox",lwd=3,gborder="blue",mdborder="red",col=NA,lims=c(-2^1000,-2^1000,2^1000,2^1000)){

#Currently:
#'fromtype' may take values:  supgem, cart, or oldbox
#indicating whether the box was 

#for prim you can specify individual boxes, for cart you accept a complexity level
#when prim=TRUE, allows you to set boxnum=NA and get all boxes plotted

if(fromtype=="cart"){
  for (i in 1:length(sdobj$cddate.bsets[[boxnum]])){
  
    o <- sdobj$cddate.bsets[[boxnum]][[i]][[2]]
    xdrn <- which(sdobj$cddate.bsets[[boxnum]][[i]][[1]]==xdim)
    ydrn <- which(sdobj$cddate.bsets[[boxnum]][[i]][[1]]==ydim)
  
    xleft <- o[xdrn,1]
    ybottom <- o[ydrn,1]
    xright <- o[xdrn,2]
    ytop  <- o[ydrn,2]
  
    if(length(sdobj$cddate.bsets[[boxnum]][[i]][[1]])>2){border <- mdborder}
    else{border <- gborder}
  
    rect(max(lims[1],xleft), max(lims[2],ybottom), min(xright,lims[3]), min(ytop,lims[4]),lwd=lwd  ,border=border,col=col)
  
    }
  
  }
  
else if(fromtype=="supgem"){

  if (is.na(boxnum)){boxnum <- 1:length(sdobj$cddate.bsets)}

  for (i in boxnum){
  
    o <- sdobj$cddate.bsets[[i]][[2]]
    xdrn <- which(sdobj$cddate.bsets[[i]][[1]]==xdim)
    ydrn <- which(sdobj$cddate.bsets[[i]][[1]]==ydim)
  
    xleft <- o[xdrn,1]
    ybottom <- o[ydrn,1]
    xright <- o[xdrn,2]
    ytop  <- o[ydrn,2]
  
    if(length(sdobj$cddate.bsets[[i]][[1]])>2){border <- mdborder}
    else{border <- gborder}
  
    rect(max(lims[1],xleft), max(lims[2],ybottom), min(xright,lims[3]), min(ytop,lims[4]),lwd=lwd  ,border=border,col=col)
  
    }
     
  }

  else if(fromtype=="oldbox"){
   
    o <- sdobj[[2]]
    xdrn <- which(sdobj[[1]]==xdim)
    ydrn <- which(sdobj[[1]]==ydim)
  
    xleft <- o[xdrn,1]
    ybottom <- o[ydrn,1]
    xright <- o[xdrn,2]
    ytop  <- o[ydrn,2]
  
    if(length(sdobj[[1]])>2){border <- mdborder}
    else{border <- gborder}
  
    #Old bad way that didn't really do limits properly
    #rect(max(lims[1],xleft), max(lims[2],ybottom), min(xright,lims[3]), min(ytop,lims[4]),lwd=lwd  ,border=border,col=col)
    
    pb <- par("usr")
    
    rect(max(pb[1], xleft), max(pb[3], ybottom), min(pb[2],xright), min(ytop, pb[4]), lwd = lwd, border = border, 
    col = col)
  
  }

}
