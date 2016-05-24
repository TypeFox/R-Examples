`join.freq` <-
function(histogram,join){
classes<-histogram$breaks
frec<-histogram$counts
frec[join[1]]<-sum(frec[join])
join<-join[-1]
classes<-classes[-join]
frec<-frec[-join]
h<-graph.freq(classes,counts=frec,plot=FALSE)
return(h)
}

