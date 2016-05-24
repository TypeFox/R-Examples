SimMLDS <- function(Trials, Scale, Sigma, n = 1){
	p <- dim(Trials)[2]
	cnt <- if (p == 4) c(1, -1, -1, 1) else c(1, -2, 1)
	Del <- drop(matrix(Scale[Trials], ncol = p) %*% 
		cnt/Sigma) + rnorm(nrow(Trials) * n)
   Resp <- Del > 0
#   d <- data.frame(Resp = Resp[1:nrow(Trials)], S = Trials)
#   if(p==4) as.mlds.df(d) else as.mlbs.df(d) 
   if(n == 1) {
   		d <- data.frame(Resp = Resp[1:nrow(Trials)], 
   							S = Trials)
   		if(p==4) as.mlds.df(d) else as.mlbs.df(d) 
   	} else {
   		lapply(seq_len(n), function(nn) {
   			d <- data.frame(Resp = 
   			Resp[seq(((nn - 1) * nrow(Trials) + 1), nn * nrow(Trials))], S = Trials)
   			if(p==4) as.mlds.df(d) else as.mlbs.df(d)
   			})
   	}
}
