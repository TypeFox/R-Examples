
library(lsa)
demo(lsa_landauer)

plot(landauerSpace$dk[,1:2]*landauerSpace$sk[1:2], pch=17, col="darkgreen")
points(landauerSpace$tk[,1:2]*landauerSpace$sk[1:2], pch=23, col="darkred")

text(landauerSpace$tk[,1:2]*landauerSpace$sk[1:2],rownames(landauerSpace$tk), col="darkred")
text(landauerSpace$dk[,1:2]*landauerSpace$sk[1:2],rownames(landauerSpace$dk), col="darkgreen")

