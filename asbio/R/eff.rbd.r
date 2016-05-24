eff.rbd <- function(lm){
alm<-anova(lm)
block <- alm[rownames(alm)=="block",]
trt <- alm[rownames(alm)!="block"&rownames(alm)!="Residuals",] 
error <- alm[rownames(alm)=="Residuals",] 
nb <- as.numeric(block[1])+1
r <- as.numeric(trt[1])+1
MSB <- as.numeric(block[3])+1
MSE <- as.numeric(error[3])+1

(((nb - 1)* MSB)+(nb * (r - 1) * MSE))/((nb * r) - 1 * MSE)
}

