winmean <- function(x, tr = 0.2, na.rm = FALSE){
if(na.rm)x=elimna(x)
winmean<-mean(winval(x,tr))
winmean
}