detsi<-function(pysty,lehtilkm){
#Hakee pysty-vektorista "pysty" sen indeksin, jonka kohdalla
#esiintyy luku "lehtilkm"
#
lkm<-length(pysty)
i<-1
while ((pysty[i]!=lehtilkm) && (i<=lkm)) i<-i+1
return(i)
}
