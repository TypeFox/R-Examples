"plotellipse" <- function(mat,alpha=0.05,coord=c(1,2),eig,cex=1,color=NULL,title=NULL){

if (length(color)==0) color = c("black","red","green3","blue",
  "cyan","magenta","darkgray","darkgoldenrod","darkgreen","violet","turquoise","orange","lightpink","lavender","yellow","lightgreen","lightgrey",
  "lightblue","darkkhaki", "darkmagenta","darkolivegreen","lightcyan", "darkorange",
  "darkorchid","darkred","darksalmon","darkseagreen","darkslateblue","darkslategray","darkslategrey",
  "darkturquoise","darkviolet", "lightgray","lightsalmon","lightyellow", "maroon")
 plotellipseinter(mat,alpha=alpha,coord=coord,nbgroup=1,moy=TRUE,eig=eig,color=color,cex=cex,title=title)
if (length(mat$partiel)!=0) {
dev.new()
  nbgroup=length(levels(mat$partiel$simul[,ncol(mat$partiel$simul)])) / length(levels(mat$moy$simul[,ncol(mat$moy$simul)]))
  plotellipseinter(mat,alpha=alpha,coord=coord,nbgroup=nbgroup,moy=FALSE,eig=eig,color=color,cex=cex,title=title)
}
}
