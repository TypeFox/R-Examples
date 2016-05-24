getCosine <-
function(v1,v2){
	# Find the cosine value of angle of vectors v1 and v2
	if(sum(v1^2)*sum(v2^2)==0)return(0)
	return(sum(v1*v2)/sqrt(sum(v1^2)*sum(v2^2)))
    
}
