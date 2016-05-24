f.post.contrasts <- function(test.type, n.res, n.pars){


.vis <- F

if(test.type == "interaction"){


	f.vis(.id <- diag(n.pars), vis = .vis)
	f.vis(.A1 <- matrix(rep(1, n.res - 1)), vis = .vis)

	f.vis(.A2 <- -diag(n.res - 1), vis = .vis)
	f.vis(cbind(.A1, .A2), vis = .vis)
	f.vis(.A <- cbind(.A1, .A2) %x% .id, vis = .vis)
	
}

return(.A)

}
