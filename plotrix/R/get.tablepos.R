get.tablepos<-function(x) {
 plotlim<-par("usr")
 tablepos<-list()
 if(x == "bottomleft") {
  tablepos$x<-plotlim[1]
  tablepos$y<-plotlim[3]
  tablepos$xjust<-0
  tablepos$yjust<-1
 }
 if(x == "bottom") {
  tablepos$x<-(plotlim[2]+plotlim[1])/2
  tablepos$y<-plotlim[3]
  tablepos$xjust<-0.5
  tablepos$yjust<-1
 }
 if(x == "bottomright") {
  tablepos$x<-plotlim[2]
  tablepos$y<-plotlim[3]
  tablepos$xjust<-1
  tablepos$yjust<-1
 }
 if(x == "left") {
  tablepos$x<-plotlim[1]
  tablepos$y<-(plotlim[3]+plotlim[4])/2
  tablepos$xjust<-0
  tablepos$yjust<-0.5
 }
 if(x == "right") {
  tablepos$x<-plotlim[2]
  tablepos$y<-(plotlim[3]+plotlim[4])/2
  tablepos$xjust<-1
  tablepos$yjust<-0.5
 }
 if(x == "topleft") {
  tablepos$x<-plotlim[1]
  tablepos$y<-plotlim[4]
  tablepos$xjust<-0
  tablepos$yjust<-0
 }
 if(x == "top") {
  tablepos$x<-(plotlim[2]+plotlim[1])/2
  tablepos$y<-plotlim[4]
  tablepos$xjust<-0.5
  tablepos$yjust<-0
 }
 if(x == "topright") {
  tablepos$x<-plotlim[2]
  tablepos$y<-plotlim[4]
  tablepos$xjust<-1
  tablepos$yjust<-0
 }
 # if no recognizable position was passed, put it in the center
 if(x == "center" || length(tablepos)==0) {
  tablepos$x<-(plotlim[1]+plotlim[2])/2
  tablepos$y<-(plotlim[3]+plotlim[4])/2
  tablepos$xjust<-0.5
  tablepos$yjust<-0.5
 }
 return(tablepos)
}
