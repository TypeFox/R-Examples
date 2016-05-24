phuso <-
function(s, t.bar=NA, f, d){
# s = probability that a dead bat is still there after 24 hours.
# t.bar = mean persistence time in days
# f = searchers efficiency, probability that a dead bat that is present is detected by a searcher
# d = (average) number of days between two searches
if(is.na(t.bar)) t.bar<-1/(-log(s)) # expected number of days until dissapearance (for exponential distribution of laying time) 
d.tilde<- -log(0.01)*t.bar
d.star <- min(d,d.tilde)
r<-t.bar * (1-exp(-d.star/t.bar))/d.star
k<-min(1, d.tilde/d)
phuso<-r*f*k
phuso
}

