col2html = function(col, opacity=1, alpha){
  

  colL=as.list(as.data.frame(t(col2rgb(col))))
  colL$maxColorValue=255
  colL$names = names(col)
  
  
  if(missing(alpha)){
    alpha = pmin(255,round(opacity*255))
  } else {
    if(is.character(alpha))
      alpha = as.hexmode(alpha)
    if(class(alpha)=='hexmode')
      alpha = as.integer(alpha)
  }
  if(any(alpha<255))
    colL$alpha = rep_len(alpha,length(col))
  
  do.call(rgb,colL)
  
}