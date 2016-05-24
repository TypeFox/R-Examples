getNearest <- function(obs,gridx,gridy,WIN,anypoint = FALSE){
	if(any(is.na(obs))) return(NA)
	if(anypoint){
		diffx <- gridx-obs[1]
		diffy <- gridy-obs[2]
		hyps <- sqrt(diffx^2+diffy^2)
		return(which(hyps==min(hyps))[1])

	} else {
		gridin <- which(inside.owin(gridx,gridy,WIN))
		gridinx <- gridx[gridin]
		gridiny <- gridy[gridin]
		diffx <- gridinx-obs[1]
		diffy <- gridiny-obs[2]
		hyps <- sqrt(diffx^2+diffy^2)
		return(gridin[which(hyps==min(hyps))[1]])
	}
}