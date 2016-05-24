plotMAPE <-
function(MAPE.obj,cutoff,MAPE.method=c('MAPE_I','MAPE_P','MAPE_G')){

MAPE.method=match.arg(MAPE.method)

if(MAPE.method=='MAPE_I'){
	MAPE.mtx=as.matrix(subset(MAPE.obj$qvalue,MAPE.obj$qvalue$MAPE_I<cutoff))
	} else if(MAPE.method=='MAPE_P') {
	MAPE.mtx=as.matrix(subset(MAPE.obj$qvalue,MAPE.obj$qvalue$MAPE_P<cutoff))
	keep.idx=c(1:(ncol(MAPE.mtx)-3),which(colnames(MAPE.mtx)==MAPE.method))
	MAPE.mtx=MAPE.mtx[,keep.idx]
	} else if(MAPE.method=='MAPE_G') {
	MAPE.mtx=as.matrix(subset(MAPE.obj$qvalue,MAPE.obj$qvalue$MAPE_G<cutoff))
	keep.idx=c(1:(ncol(MAPE.mtx)-3),which(colnames(MAPE.mtx)==MAPE.method))
	MAPE.mtx=MAPE.mtx[,keep.idx]
} else {stop('Please check: the MAPE.method should be one of MAPE_I,MAPE_P or MAPE_G') }


MAPE.bin=MAPE.mtx[,(ncol(MAPE.mtx)-2):ncol(MAPE.mtx)]
MAPE.bin=ifelse(MAPE.bin>cutoff,0,1)

if(is.null(nrow(MAPE.bin))) {
	cat('No enriched pathways were identified using the current cutoff') 
	} else {
	cut4plot=1e-5
	MAPE.mtx[MAPE.mtx<=cut4plot]=cut4plot
	MAPE.mtx=log10(MAPE.mtx)	
	sort.idx=order(MAPE.mtx[,ncol(MAPE.mtx)])
	MAPE.mtx=MAPE.mtx[sort.idx,]
	heatmap(MAPE.mtx, col = gray(0:128/128),Colv=NA,Rowv=NA,scale='none',cexCol=1,cexRow=1,main='Heatmap of enriched pathways')		
} 

}
