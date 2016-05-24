G<-
function(z, link) {
	if (link=="logit") {
			exp(z) / (1 + exp(z))
	} else if (link=="cloglog") {
			1-exp(-exp(z))
	}
}

