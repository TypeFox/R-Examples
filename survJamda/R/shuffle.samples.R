shuffle.samples <-
function(n, censor, train.nb)
{
	for (i in 1:5){
        	o <- sample(1:n)
        	train.ind = o[1:train.nb]
		test.ind = o[(train.nb+1):n]
			
		if(sum(censor[train.ind]) > 1)
			break
	}
	if (sum(censor[train.ind]) <= 1)
		stop("\rToo few non-censored patients. Choose a higher number of samples for the training set., call. = F")

	return (list(train.ind = train.ind, test.ind = test.ind))
}

