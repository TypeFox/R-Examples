plot.blca.em <-
function(x, which=1L, main="", ...){
  #class(x)<- "blca"
  #print("NextMethodUsed")
  show<- rep(FALSE, 5)
  show[which]<- TRUE
  if(any(show[3:4])){ 
  	warning("Density estimates not supported for em estimates.")
  	show[3:4]<- FALSE
  	 }
  which<- c(1:4,7)[show]
  NextMethod("plot")
}
