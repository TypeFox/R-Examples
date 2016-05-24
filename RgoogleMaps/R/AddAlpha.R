AddAlpha <- structure(function#add alpha level to color that lacks one
###add alpha level to color that lacks one
(
  plotclr, ##<< color to be modified
  alpha=0.5, ##<< alpha level
  verbose=0 ##<< level of verbosity
){
  tmp <- col2rgb(plotclr, alpha=alpha)
  tmp[4,] = round(alpha*255)
  for (i in 1:ncol(tmp)){
    plotclr[i] = rgb(tmp[1,i], tmp[2,i], tmp[3,i], tmp[4,i], maxColorValue = 255)
  }
  return(plotclr)
### modified color with alpha value
}, ex = function(){

  #example: 
  #require(RColorBrewer)
  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    plotclr <- RColorBrewer::brewer.pal(8,"YlOrRd")
    plotclr = AddAlpha(plotclr,0.5)
  } else {
    print("package RColorBrewer must be installed for this example")
  }
  
})
