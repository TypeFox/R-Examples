scoreMDscan <-function(input,k,matriu,direction){
	
	Score_MDscan<-vector("numeric",length=(k-ncol(matriu)+1))
	for (i in c(length(input):1)) {
		for (j in c((nrow(input[[i]]$Secuencia)):1)) {
   			if (( input[[i]]$Secuencia[j,2]=="f" | input[[i]]$Secuencia[j,2]=="r") & is.na(input[[i]]$Secuencia[j,3])==FALSE & as.numeric(input[[i]]$Secuencia[j,3])<=(k-ncol(matriu)+1)) {
				Score_MDscan[as.numeric(input[[i]]$Secuencia[j,3] )]<-as.numeric(input[[i]]$Score)}
			}
		}
	return(Score_MDscan)
}

