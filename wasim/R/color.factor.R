color.factor<-function(color, value, max){
    l.color<-length(color)
    t.rgb=rep(col2rgb(color),length(value))
    t.rgb[is.na(t.rgb)]<-255
    dim(t.rgb)<-c(3,l.color*length(value))
    value[is.na(value)]<-0
    t.rgb<-255-((255-t.rgb)*rep(value, each=3)/max)
    return(rgb(red=t.rgb[1,], green=t.rgb[2,], blue=t.rgb[3,], max=255))

}

