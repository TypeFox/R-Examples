permute <- function(nn,mm,N,L,Z,tmp,NDIi){
	if (nn == 1){
		NDIi = rep(0,N*Z)
		dim(NDIi) <- c(N,Z)
	}
	for (elem in 1:L[nn]){
		tmp[nn] = elem
       if (nn == N){
    		NDIi[,mm] = tmp
			mm = mm + 1
       } else {
			res = permute(nn+1,mm,N,L,Z,tmp,NDIi)
			mm = res$listmm
			NDIi = res$listNDIi
       }
	}
	res = list(listmm = mm, listNDIi = NDIi)
	return(res)
}