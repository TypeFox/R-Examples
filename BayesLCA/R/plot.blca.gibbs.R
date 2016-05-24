plot.blca.gibbs <-
function(x, which=1L, main="", ...){
  #class(x)<- "blca"
  #print("NextMethodUsed")
  show<- rep(FALSE, 5)
  show[which]<- TRUE
  which<- c(1:4,8)[show]
  NextMethod("plot")
}
