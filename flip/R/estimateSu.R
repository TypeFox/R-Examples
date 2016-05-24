###############################################33
.estimateSuMultiILS = function(Y,Z,S,tol=0.001,max.int=10){


	tr =function(x){return(sum(diag(as.matrix(x))))}
	## si: subject's standard errors
	n = nrow(Y) ; p = ncol(Y)
	Su = array(0,dim=c(p,p))
	TR = array(0,dim=c(max.int,1))
	bb=2 ; ID = 0
	si = array(0,dim=c(n,1))

	while(bb <= max.int & ID == 0){
		for(i in 1:n){#print(tr(S[i,,]))
                        si[i] = tr(S[i,,])+tr(Su)
		}
    
    #print(si)
		V = diag(as.vector(1/sqrt(si))) ; # W =  diag(as.vector(sqrt(si)))
		ys = V%*%Y ; Zs = V%*%Z
		Hs = diag(n)-Zs%*%solve(t(Zs)%*%Zs)%*%t(Zs)
    #force the symmetry (errors on approx)
    Hs = (Hs+t(Hs))/2
		Es=eigen(Hs)
		Ls = Es$vectors%*%diag(Es$values)
		Rs = t(Ls)%*%ys  				#; Rr = W%*%Rs
		G = t(Ls)%*%V
		w = diag(t(G)%*%G)
		SS = array(0,dim=c(p,p))

		for(i in 1:n){SS = SS + w[i]*S[i,,]}
		#Rs = W%*%Rs
		SSu = (t(Rs)%*%Rs-SS)/sum(w) 
		
		########## simmetrizza per non avere PARTI IMMAGINARIE dovute a scarsa approssimazione
    E = eigen((SSu+t(SSu))/2)
		E$values[E$values<0]=0
		if(length(E$values)==1) E$values=matrix(E$values)
		
		A = diag(E$values) 
		TR[bb] = sum(A)
		Su = E$vectors%*%A%*%t(E$vectors)  ; 
		if(abs(TR[bb]-TR[bb-1])<tol){ID = 1}
		bb = bb+1
	}
  colnames(Su)=rownames(Su)=colnames(Y)
  attr(Su,"n.iter")=bb-1
	attr(Su,"TR")=TR
	return(Su)
}