#' Determine whether {f1>=c1} a subset of {f2>=c2}.
#'
#' @param x x-Coordinates of the grid on which the functions are observed.
#' @param y y-Coordinates of the grid on which the functions are observed.
#' @param f1 A matrix of dimension c(length(x),length(y)) defining the first 
#'           function.
#' @param f2 A matrix of dimension c(length(x),length(y)) defining the second 
#'           function.
#' @param c1 The first level.
#' @param c2 The second level.
#' @return A boolean indicating whether {f1>=c1} a subset of {f2>=c2}.
SubsetContour = function(x,y,f1,c1,f2,c2){
  cont1 = contourLines(x,y,f1,levels=c1,nlevels=1)
  cont2 = contourLines(x,y,f2,levels=c2,nlevels=1)
  
  if(all(f1<c1)) return(TRUE)
  if(all(f1>c1) & all(f2>=c2)) return(TRUE)
  if(all(f1>c1) & !all(f2>=c2)) return(FALSE)
  if(all(f2>c2)) return(TRUE)
  if(all(f2<c2)) return(FALSE)
  
  cont1_x = c()
  cont1_y = c()
  for(l in 1:length(cont1)) {cont1_x = c(cont1_x,cont1[[l]]$x); cont1_y = c(cont1_y,cont1[[l]]$y)}
  cont1 = cbind(cont1_x,cont1_y)
  
  cont2_x = c()
  cont2_y = c()
  for(l in 1:length(cont2)) {cont2_x = c(cont2_x,cont2[[l]]$x); cont2_y = c(cont2_y,cont2[[l]]$y)}
  cont2 = cbind(cont2_x,cont2_y)
  
  min(fields::interp.surface(list(x=x,y=y,z=f2),cont1))>=c2 & max(fields::interp.surface(list(x=x,y=y,z=f1),cont2))<=c1
}
