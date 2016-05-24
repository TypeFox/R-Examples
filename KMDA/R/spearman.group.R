spearman.group <-
function(data,threshold){
	p=nrow(data)
	cormatrix=matrix(NA,p,p)
	for(i in 1:p){
	cormatrix[i,i]=1
	}
	for(i in 1:(p-1)){
	for(j in (i+1):p){
		cormatrix[i,j]=cormatrix[j,i]=cor(data[i,],data[j,],method='spearman')
			}
	}
	label=c(1:p)
	for(i in 2:p){
		for(j in 1:(i-1)){
			if(cormatrix[i,j]>threshold) label[i]=label[j]
		}
	}
	return(label)
	## Rows with the same label are in the same group.
}
