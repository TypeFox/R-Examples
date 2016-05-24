cfrekv<-function(levels,arvot){
#laskee tasojoukon osien frekvenssit
#arvo on reaaliluku
#kumu on kork*n-matriisi, n saraketta, kuvaa kork kpl:tta tasojoukon osia
#1 jos vastaava data-matriisin rivin indikoima pallo kuuluu tasojouon osaan
#muodostetaan matriisi, jonka 1. sarakkeessa "arvo", 
#2. sarakkeessa kunkin tasojoukon osan frekvenssi  
#ts. laskettu kuinka monesta pallosta tasojoukko on yhdistetty
#
tasolkm<-length(levels[,1])     #levels:n rivien maara
frek<-matrix(0,tasolkm,1)
a<-1
while (a<=tasolkm){
   frek[a]<-sum(levels[a,]*arvot)
   a<-a+1 
}
return(t(frek))
}







