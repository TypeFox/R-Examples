#' @export
dis.cost<-function(cost){
d.cost <- cost
d.cost$time<-ceiling(cost$time/365)
return(d.cost)
}