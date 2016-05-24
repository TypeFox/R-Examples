compareImages <-function(Im1,Im2, testsize=min(nrow(Im1),nrow(Im2)),alpha=0.05,...){

dimIm1 <- dim(Im1)
dimIm2 <- dim(Im2)

if(dimIm1[1] > dimIm2[1]){
	Im1 <- Im1[1:dimIm2[1],1:dimIm2[1]]
}

if(dimIm2[1] > dimIm1[1]){
	Im2 <- Im2[1:dimIm1[1],1:dimIm2[1]]
}

CombIm1Im2 <- cbind(Im1,Im2)

cropCombIm <- cropimage(CombIm1Im2, newsize=testsize)

montagep <- TOS2D(cropCombIm,...)$p.value

montageres <- (montagep > alpha)

return(montageres)

}

