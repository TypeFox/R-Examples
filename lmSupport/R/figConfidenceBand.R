figConfidenceBand <- function(X, Y, CILo, CIHi, Color = 'gray'){
  
  #Color = col2rgb(Color, TRUE)
  #Color[4] = .15  #make color semi transparent
  polygon(c(X, rev(X)), c(CILo, rev(Y)),col = Color, border = NA)
  polygon(c(X, rev(X)), c(CIHi, rev(Y)),col = Color, border = NA)
 
}