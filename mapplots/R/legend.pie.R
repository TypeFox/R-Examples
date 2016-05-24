legend.pie <-
function(x,y=NULL,z=NULL,labels,radius=1,bty='o',mab=1.2,bg=NULL,inset=0,...){
  if(is.null(z)) z <- rep(1,length.out=length(labels))
  box <- legend.box(x,y,radius,mab,inset)
  if(bty=='o') rect(box[1],box[2],box[3],box[4],col=bg)
  x <- (box[1]+box[3])/2
  y <- box[4]+mab*radius
  add.pie(z,x,y,labels,radius,...)
}

