sginitc <-
function(G,sg,pig,p,n,x){
	sgc <- matrix(0,p,p)
	for(g in 1:G){
		sgc <- sgc + pig[g]*sg[,,g]
	}
	sgc
}
