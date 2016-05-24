evanor<-function(x){

d<-length(x) 
eta<-sum(x^2)     #vektorin x pituuden nelio
normvakio<-(sqrt(2*pi))^{-d}
tulos<-exp(-eta/2)*normvakio
return(tulos)
}
