plot.blca.vb <-
function(x, which=1L, main="", ...){
  #class(x)<- "blca"
  #print("NextMethodUsed")
  show<- rep(FALSE, 5)
  show[which]<- TRUE
  which<- c(1:2,5:6,9)[show]
  NextMethod("plot")
}
