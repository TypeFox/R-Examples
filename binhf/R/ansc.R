"ansc" <-
function(x,binsize){
y<-asin(sqrt((x+3/8)/(binsize+3/4)))
y <- y*sqrt(4*(binsize+.5))
y
}

