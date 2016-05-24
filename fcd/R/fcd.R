fcd <-
function(A, K = 2, nlambda = 1e+3, lambda.min.ratio = 1e-5, alpha = 0.8, scale = FALSE){
	
	### get fcd.start object first
	t1 = proc.time()
	obj = fcd.start(A = A, K = K, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, alpha = alpha, scale = scale)
	t2 = proc.time() - t1
	beta.combind = obj$beta.combind
	iso.seq = obj$iso.seq
	lambda.list = obj$lambda.list
	
	### clustering
	t3 = proc.time()
	temp.cluster = fcd.cluster(obj = beta.combind, K = K)
	t4 = proc.time() - t3
	
	### criteria
	t5 = proc.time()
	temp.criteria = fcd.criteria(A = A, fcd.cluster = temp.cluster, K = K, iso.seq = iso.seq)
	t6 = proc.time() - t5
	
	### get cluster
	t7 = proc.time()
	temp = get.cluster(A = A, iso.seq = iso.seq, criteria.list = temp.criteria, clusters.list = temp.cluster)
	t8 = proc.time() - t7
	
	time = list(start = t2[1], cluster = t4[1], criteria = t6[1], getcluster = t8[1])
	
	### return
	return(list(beta.combind = beta.combind, iso.seq = iso.seq, cluster.list = temp.cluster, criteria.list = temp.criteria, final.ratio.cluster = temp$final.ratio.cluster, ratio.location = temp$ratio.location, final.normalised.cluster = temp$final.normalised.cluster, normalised.location = temp$normalised.location, lambda.list = lambda.list, time = time))
}
