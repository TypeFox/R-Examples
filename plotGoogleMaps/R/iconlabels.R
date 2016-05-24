iconlabels<-function(attribute,
                     colPalette=NULL,
                     at=NULL,
                     height=10,
                     icon=FALSE,
                     scale=0.6
                     
  ){
  
  if(!is.null(colPalette)){
    rgbc<-col2rgb(colPalette)
    colPalette<-apply(rgbc,2,function(x) rgb(x[1],x[2],x[3],maxColorValue=255))}
  
  if(length(colPalette)==length(attribute)) {
    color <- sub('\\#','',as.character(colPalette))
} else{
  x=PolyCol(attribute,colPalette,at)
  color <- sub('\\#','',as.character(x$cols))
  }
  
  vals <- paste("http://chart.apis.google.com/chart?chst=d_text_outline&chld=", color, "|",height,"|h|000000|b|", attribute, sep="")
  
  if(icon){
    vals <- paste("http://chart.apis.google.com/chart?chst=d_map_spin&chld=",scale,"|0|",color,"|",height,"|b|", attribute, sep="")
  }
  
  
  
  return(vals)
}

