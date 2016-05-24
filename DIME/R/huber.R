`huber` <-
function(input, co = -1.345, shape = c('full','lower','upper'))
{
	input <- unlist(input);
	len <- length(input);
  shape <- match.arg(shape)
	input <- (input - mean(input))/sd(input)
	change <- switch(shape,
		full = which(abs(input) > abs(co)),
		lower = which(input <= co),
		upper = which(input >= abs(co))
	)
  	if (length(change)<1) input <- rep(1,len) ;
    input[change] <- abs(co)/abs(input[change])
  	input[-change] <- 1
	return(input)
}