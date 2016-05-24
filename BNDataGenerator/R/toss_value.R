# Written by Jae-seong Yoo 20141208
# developed of tosscoin in prob package

# Toss Value

toss_value = function (times, num_of_cases, makespace=FALSE) 
{
	mat_values = merge("Value", c(1:num_of_cases))
	
	temp = list()
	for (i in 1:times) {
		temp[[i]] = paste(mat_values[,1], mat_values[,2], sep="")
	}
	res = expand.grid(temp, KEEP.OUT.ATTRS = FALSE)
	names(res) = c(paste(rep("toss", times), 1:times, sep = ""))
	if (makespace) 
		res$probs = rep(1, 2^times)/2^times
		
	return(res)
}
