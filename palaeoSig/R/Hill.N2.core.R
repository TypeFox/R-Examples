Hill.N2.core <-
function(spp){
 hn2<-Hill.N2(spp, margin=1)
 return(quantile(hn2,prob=c(0,.25,.5)))
}

