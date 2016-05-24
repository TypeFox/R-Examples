allelicOR <-
function(xtable){
	#a 2x3 table
	xtable2 = cbind(xtable[,1]*2 + xtable[,2], xtable[,2] + xtable[,3]*2)
	xtable2[1,1]*xtable2[2,2]/(xtable2[1,2]*xtable2[2,1])
}

