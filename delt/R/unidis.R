unidis<-function(ulot){
#tasainen jakauma vec:n elementeille
#ulot=dim(vec)
#tasainen jakauma joukossa {1,2,...,ulot}
#
arpa<-runif(1)
ele<-ceiling(ulot*arpa)
return(ele)
}
