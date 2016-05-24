t3pval<-function(cmat,tmeans,v,h){
alph<-c(1:99)/100
for(i in 1:99){
irem<-i
chkit<-johan(cmat,tmeans,v,h,alph[i])
if(chkit$teststat>chkit$crit)break
}
p.value <- irem/100
        if(p.value <= 0.1) {
                iup <- (irem + 1)/100
                alph <- seq(0.001, iup, 0.001)
                for(i in 1:length(alph)) {
                        p.value <- alph[i]
                        chkit<-johan(cmat,tmeans,v,h,alph[i])
if(chkit$teststat>chkit$crit)break
                }
        }
  if(p.value <= 0.001) {
                alph <- seq(0.0001, 0.001, 0.0001)
                for(i in 1:length(alph)) {
                        p.value <- alph[i]
chkit<-johan(cmat,tmeans,v,h,alph[i])
if(chkit$teststat>chkit$crit)break
                }
        }
p.value
}
