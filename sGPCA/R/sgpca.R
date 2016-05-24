sgpca = function(X,Q,R,K=1,lamu=0,lamvs=0,posu=FALSE,posv=FALSE,threshold=1e-7,maxit=1000,full.path = FALSE){
	
	
## error checking
	
	if(class(X) != "matrix" && class(X) != "dgCMatrix"){
		stop("X must be a matrix of class 'matrix'  or 'dgCMatrix' ")
	}
	
	if(class(Q) != "matrix" && class(Q) != "dgCMatrix"){
		stop("Q must be a matrix of class 'matrix'  or 'dgCMatrix' ")
	}
	
	
	if(class(R) != "matrix" && class(R) != "dgCMatrix"){
		stop("R must be a matrix of class 'matrix'  or 'dgCMatrix' ")
	}
	
	
	
	n = dim(X)[1]
	p = dim(X)[2]
	
	qn = dim(Q)[1]
	pn = dim(Q)[2]
	
	rn = dim(R)[1]
	rp = dim(R)[2]
	
	if(qn != pn){
		stop("Q must be a square matrix")
	}
	
	if(qn != n){
		stop("Q must be an n x n matrix")
	}
	
	if(rn != rp){
		stop("R must be a square matrix")
	}
	
	if(rn != p){
		stop("R must be a p x p matrix")
	}
	
	
	K = floor(K)
	if(K <1){
		stop("k must be an integer greater than 0")
	}
	K = as.integer(K)
	
	if(lamu < 0){
		stop("lamu must be greater than 0")
	}
	
	minlamv = min(lamvs)
	if(minlamv < 0){
		stop("lamvs cannot be less than 0")
	}
	lamvs = as.matrix(sort(lamvs))
	
	if(threshold < 0){
		stop("threshold must be greater than 0")
	}
	
	if(maxit < 0){
		stop("maxit must be greater than 0")
	}
	
## Initialize  column compressed format
	ir_q = 0
	jc_q = 0
	pr_q = Q
	
	ir_r = 0
	jc_r = 0
	pr_r = R
	
	QisSparse = 0
	RisSparse = 0
	
	
	if(class(Q)=="dgCMatrix"){
		
		ir_q = Q@i
		jc_q = Q@p
		pr_q = Q@x
		QisSparse = 1
	}
	
	if(class(R) == "dgCMatrix"){
		ir_r = R@i
		jc_r = R@p
		pr_r = R@x
		RisSparse = 1
	}
	
	QisSparse = as.integer(QisSparse)
	RisSparse = as.integer(RisSparse)
	
	jc_q= as.integer(jc_q)
	ir_q = as.integer(ir_q)
	jc_r = as.integer(jc_r)
	ir_r= as.integer(ir_r)
	
	posu = as.integer(posu)
	posv = as.integer(posv)
	maxit = as.integer(maxit)
	
## Manage posu posv option	
	pass.posu = as.integer(0)
	if(posu){
		pass.posu = as.integer(1)
	}
	
	pass.posv = as.integer(0)
	if(posv){
		pass.posv = as.integer(1)
	}
## Run Analysis	
	
	
	r = length(lamvs)
	Us = matrix(nrow = n,ncol = r*K)
	Vs = matrix(nrow = p,ncol = r*K)
	Ds = rep(0,r*K)
	
	Xhat = X
	row.index = 1
	bic.index = rep(0,K)
	bics = rep(0,K)
	optlams = rep(0,K)
	
## computation
	
	
	
	for (i in 1:K){
		for (j in 1:r){
			
			# index of component
			old.index = row.index
			row.index = (i-1)*r + j
			
			#calculate component
			if(j==1){
				start = gmdLA(Xhat,Q,R,1,n,p)
				gpmf.results = .Call("gpmfSparse",Xhat,pr_q,pr_r,jc_q,ir_q,jc_r,ir_r,RisSparse,QisSparse,lamu,lamvs[j],posu,posv,threshold,maxit,as.matrix(start$u),as.matrix(start$v))
			}else{
				gpmf.results = .Call("gpmfSparse",Xhat,pr_q,pr_r,jc_q,ir_q,jc_r,ir_r,RisSparse,QisSparse,lamu,lamvs[j],posu,posv,threshold,maxit,as.matrix(Us[,old.index]),as.matrix(Vs[,old.index]))
			}
			
			#Store component
			Us[,row.index]=gpmf.results[[1]]
			Vs[,row.index]=gpmf.results[[2]]
			Ds[row.index] =gpmf.results[[3]]
			
			#calculate bic for component
			df = sum(as.numeric(Vs[,row.index]!=0))
			temp1 = Xhat - Ds[row.index]*Us[,row.index]%*%t(Vs[,row.index])
			temp2 = as.matrix(Q %*% temp1 %*% R %*% t(temp1))
			bic = log(sum(diag(temp2))/(n*p)) + (log(n*p)/(n*p))*df

			# get minimum bic and index of minimum
			if(j == 1){
				minbic = bic
				bic.index[i] = row.index
				bics[i] = bic
				optlams[i] = lamvs[j]
			} else if(bic < minbic){
				minbic = bic
				bics[i] = bic
				bic.index[i] = row.index
				optlams[i] = lamvs[j]
			}
			
		}
		
		#deflate data matrix
		Xhat = Xhat - Ds[bic.index[i]]*Us[,bic.index[i]]%*%t(Vs[,bic.index[i]])
	}
	
	
	#calculate proportion of variance
	var = calculateVariance(X,as.matrix(Q),as.matrix(R),Us[,bic.index],Vs[,bic.index],K)
	
## Return results	
	
	if(full.path){
		return(list("U"=Us,"V" = Vs,"D" = Ds,"cumulative.prop.var" = var,"bics" = bics, "optlams" = optlams))
	}else{
		return(list("U"=Us[,bic.index],"V" = Vs[,bic.index],"D" = Ds[bic.index],"cumulative.prop.var" = 
var,"bics" = bics, "optlams" = optlams))
	}

		
	
}

