## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

writet <- function(f,ni,xt){ 
for(i in 1:size(ni,1) ){
	for(j in 1:ni[i]){
		fprintf(f,'%g ',xt[i,j])
	}
	fprintf(f,'\n')
}
return( ) 
}
