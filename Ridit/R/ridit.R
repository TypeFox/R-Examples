ridit <-
function(x,g,ref=NULL){
if(length(dim(x))==2) {
	# table, data.frame, matrix, array with 2 dimensions
	if(length(g)!=1 || g<1 || g>2)	
	stop("when x is a crosstabualtion, g must be 1 or 2 that specify dimension")
	m=ict(x) 
	if(g==1) output=ridit.raw(x=m$d2,g=m$d1,ref) else output=ridit.raw(x=m$d1,g=m$d2,ref)
}	# end if(length(dim(x))==2)
else if(length(dim(x))==0){
	
	output=ridit.raw(x,g,ref)  # factor, vector, array with 1 dimension, list, formula and other objects!!!
} # end if(length(dim(x))==0)
else stop("dimension of x is incorrect")  #  array with 3 or more dimensions
output
}
