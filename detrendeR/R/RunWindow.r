RunWindow=function(rwl, start=FirstYear(rwl), winLength=50){
a<-as.integer(rownames(rwl))
b<-a>=start & a<=start+winLength-1
win<-subset(rwl, subset=b)
win<-na.omit(t(win))
win<-as.data.frame(t(win))
return(win)
}
#RunWindow(rwl, start=1901)