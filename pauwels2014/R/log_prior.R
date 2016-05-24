log_prior <-
function(theta){
	## Prior independant of the dynamical charaxcteristics
	- 1/2*sum((theta - 50)^2/ 10000 + log(100) ) / length(theta)
}
