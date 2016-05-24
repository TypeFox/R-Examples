fcd.start <-
function(A, K = 2, nlambda = 1e+3, lambda.min.ratio = 1e-5, alpha = 0.8, scale = FALSE){
	
	### preparation 
	iso.A = isolate(A)
    iso.seq = iso.A$isolate
    noniso.seq = iso.A$nonisolate
    A.noniso = A[noniso.seq, noniso.seq]
    
    ### Laplacian and things related
    L = laplacian(A.noniso)
	eig.L = eigen(L, symmetric = T)
	n = dim(L)[2]
	Y = matrix(0, nrow = n, ncol = K)
	for(i in 1: K){
		Y[, i] = eig.L$vectors[, (n - i + 1)]
	}

	### edge incidence matrix and things related
	B = fcd.trans(A.noniso)
	B.inv = ginv(B)
	
	### glmnet and things related
	fit.glmnet = list()
	beta.mat = list()
	length.vec = c()
	lambda.list = list()
	for(i in 1: K){
		fit.glmnet[[i]] = glmnet(x = t(B), y = Y[, i], nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, alpha = alpha)
		beta.mat[[i]] = apply(fit.glmnet[[i]]$beta, 2, function(x){B.inv %*% x})
		length.vec[i] = dim(beta.mat[[i]])[2]
		lambda.list[[i]] = fit.glmnet[[i]]$lambda
	}
	min.length = min(length.vec)

	### combine columns
	beta.combind = list()
	for(i in 1: min.length){
		beta.combind[[i]] = beta.mat[[1]][, i]
		for(j in 2: K){
			beta.combind[[i]] = cbind(beta.combind[[i]], beta.mat[[j]][, i])
		}
		if(scale) beta.combind[[i]] = scale(beta.combind[[i]], center = F)	
	}
	
	### return
	return(list(beta.combind = beta.combind, iso.seq = iso.seq, lambda.list = lambda.list))

}
