Cond.prob.MSAR <-
function(data,theta,yrange=NULL,covar.emis=NULL,covar.trans=NULL){
	if (is.null(yrange)) {
		yrange = seq(min(data,na.rm=TRUE),max(data,na.rm=TRUE),length.out=100)
	}
    lyr = length(yrange)
	T = dim(data)[1]
	N.samples = dim(data)[2]
	M = attributes(theta)$NbRegimes
	order=attributes(theta)$order
	scale = matrix(1,1,T)
	loglik = 0
	alpha = matrix(0,M,T)
	prob = array(0,c(lyr,T,N.samples))
	prior = theta$prior
	if (substr(attributes(theta)$label,1,1)=="H") {transition = array(theta$transmat,c(M,M,T))
	} else if ( substr(attributes(theta)$label,1,1)=="N") {transition = attributes(theta)$nh.transitions(covar.trans,theta$par.trans,theta$transmat)}
#	transmat = theta$transmat
	Yhat = matrix(0,T,N.samples)
	a = array(0,c(M,lyr,T,N.samples))
	for (ex in 1:N.samples) {
		p1 = matrix(0,2,1)
		py = matrix(0,lyr,M)
		for (o in (1:order)) {
			for (m in 1:M) {
				p1[m] = p1[m]+pdf.norm(matrix(data[1,ex,1],1,1),theta$A0[m],matrix(theta$sigma[[m]]))*prior[m]
				for (k in 1:lyr) {
					py[k,m] = py[k,m]+pdf.norm(matrix(yrange[k],1,1),theta$A0[m],matrix(theta$sigma[[m]]))*prior[m]
				}
			}
		}
		alpha[,1] = p1
		alpha[,1]  = normalise(alpha[,1])
		obslik.y = matrix(0,M,T)
		for (o in 1:order) {
			y = data[,ex,]
			g = emisprob.MSAR(y,theta,covar.emis[,ex,])
			obslik.y[,(order+1):dim(obslik.y)[2]] = g
    			if (o>1) {tmp1 = (t(transition[,,o]) %*% alpha[,o-1])
				tmp = tmp1 * obslik.y[,o]
				alpha[,o] = normalise(tmp)}
			
			for (k in 1:lyr) {
				y[o] = yrange[k] # On remplace la valeur au temps t par le y qui nous interesse
    			g = emisprob.MSAR(y,theta,covar.emis[,ex,]) # a optimiser
    			#
    			if (o>1) {
    				obslik.y[,1:(o-1)] = p1
    				}
    			obslik.y[,o:order] = py[k,]	
    				a = normalise(py[k,])
				beta = matrix(1,M,1)
				# BACKWARD .........................................................
				for (tb in seq(T-1,o,-1)) { # On s'arrete au temps t? 
					b = beta * obslik.y[,tb+1];
					#beta = normalise((transmat2 %*% b)) #P(Y_t=y|Y_{t+1:T}=y_{t+1:T}) 					
					beta = normalise((transition[,,tb+1] %*% b)) #P(Y_t=y|Y_{t+1:T}=y_{t+1:T}) 

				}
				prob[k,o,ex] = sum(a * beta)
			}
			y = prob[,o,ex] 
			prob[,o,ex] =  prob[,o,ex]/sum(diff(yrange)*(y[-length(y)]+y[-1])/2)
			y = yrange*prob[,o,ex]
			Yhat[o,ex] =  sum(diff(yrange)*(y[-length(y)]+y[-1])/2)
					prob[,o,ex] =  prob[,o,ex]/sum(prob)

		}
        g = emisprob.MSAR(data[,ex,],theta,covar.emis[,ex,]) # pour calculer les alpha
        obslik=matrix(0,M,T)
    	obslik[,(order+1):T] = g
		for (t in max((order+1),2):(T-1)) {
			y = data[,ex,]
			tmp1 = (t(transition[,,t]) %*% alpha[,t-1])
			tmp = tmp1 * obslik[,t]
			alpha[,t] = normalise(tmp)
    		for (k in 1:lyr) {
    			y[t] = yrange[k] # On remplace la valeur au temps t par le y qui nous interesse
    			g = emisprob.MSAR(y,theta,covar.emis[,ex,]) # a optimiser
    			obslik.y[,(order+1):dim(obslik.y)[2]] = g
    			# FORWARD ..............................................................
				tmp = tmp1 * obslik.y[,t]
				a = normalise(tmp1)*obslik.y[,t] #P(Y_t=y|Y_{1:t-1}=y_{1:t-1},S) 
#a[k,t,ex] = sum(tmp) # On doit trouver la meme chose dans le code forecast	
				beta = matrix(0,M,T) 
				beta[,T] = matrix(1,M,1)
				# BACKWARD .........................................................
				for (tb in seq(T-1,t,-1)) { # On s'arrete au temps t? 
					b = beta[,tb+1] * obslik.y[,tb+1];
					beta[,tb] = normalise((transition[,,tb+1] %*% b)) #P(Y_t=y|Y_{t+1:T}=y_{t+1:T}) 
				}
				prob[k,t,ex] = sum(a * beta[,t])
			}
			y = prob[,t,ex] 
			prob[,t,ex] =  prob[,t,ex]/sum(diff(yrange)*(y[-length(y)]+y[-1])/2)
			y = yrange*prob[,t,ex]
			Yhat[t,ex] =  sum(diff(yrange)*(y[-length(y)]+y[-1])/2) 
					prob[,t,ex] =  prob[,t,ex]/sum(prob)

		}	
		t = T
		tmp1 = (transition[,,T] %*% alpha[,t-1])
		y = data[,ex,]
		for (k in 1:lyr) {
			y[t] = yrange[k]
			g = emisprob.MSAR(y[(t-order):t],theta,covar.emis[(t-order):t,ex,])
			a = normalise(tmp1)*g[,dim(g)[2]]
			prob[k,t,ex] = sum(a)
		}
		y = prob[,t,ex] 
		prob[,t,ex] =  prob[,t,ex]/sum(diff(yrange)*(y[-length(y)]+y[-1])/2)
		y = yrange*prob[,t,ex]
		Yhat[t,ex] =  sum(diff(yrange)*(y[-length(y)]+y[-1])/2)
		prob[,t,ex] =  prob[,t,ex]/sum(prob)
	}
    list(yrange=yrange,prob=prob,Yhat=Yhat)
}
