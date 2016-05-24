cdcsis <- function(x, y, z, thres){
	p <- dim(x)[2]
	DW.c <- sapply(1:p, function(j) cdcor.ada(x[,j], y, z))
	DCW<-unlist(DW.c[1,])
          	DCWord<-order(abs(DCW),decreasing = T)
	CDCSISind <- DCWord[1:thres]

	return(list(CDCSISind =CDCSISind, thres=thres, DC=DCW, DCord=DCWord ))
}
