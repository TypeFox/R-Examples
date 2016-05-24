


pch=function(data.matrix){
partial.capture.histories=matrix(NA,nrow=nrow(data.matrix),ncol=ncol(data.matrix))

for(i in 1:nrow(data.matrix)){
for(j in 1:(ncol(data.matrix)-1)){


partial.capture.histories[i,1]=""	
partial.capture.histories[i,j+1]=paste(data.matrix[i,1:j],collapse="")
	
	
	}}
	
out=partial.capture.histories
return(out)	
	
}
