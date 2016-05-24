plot.blca.boot <-
function(x, which=1L, main="", ...){
  #class(x)<- "blca"
  #print("NextMethodUsed")
  show<- rep(FALSE, 5)
  show[which]<- TRUE
  if(show[5]){ 
  	warning("No diagnostic plot for bootstrapping method")
  	show[5]<- FALSE
  	}
  which<- c(1:4)[show]
  NextMethod("plot")
}
