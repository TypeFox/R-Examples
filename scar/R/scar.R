scar <- function(x,y,shape=rep("l",d), family=gaussian(), weights=rep(1,length(y)), epsilon=1e-8){

        ## Some checking on inputs
        if(is.vector(x)==TRUE) x=matrix(x,ncol=1)
        y = as.vector(y)
        shape = as.vector(shape)
	if(is.matrix(x)==FALSE||is.numeric(x)==FALSE||is.vector(y)==FALSE||is.numeric(y)==FALSE) stop("Input error!")
        if((length(y)!=nrow(x))||(length(weights)!=nrow(x))) stop("Input error: x, y, d mismatch.")
        d=NCOL(x)
        if(length(shape)!=d) stop("Input error: x and shape mismatch.")
        N=NROW(x)

        ## delete some boundary points if needed
        for(j in 1:d){
            if(family$family == "poisson"){
                if((shape[j]=="ccv")||(shape[j]=="ccvin")||(shape[j]=="in")){
                    while(TRUE){
                      bound = which.min(x[,j])
                      if (y[bound]==0) {x = matrix(x[-bound,],ncol=d); y = y[-bound]; weights = weights[-bound]}
                      else break
                    }
                }
            } else if(family$family == "binomial"){
                if((shape[j]=="cvx")||(shape[j]=="cvxin")||(shape[j]=="in")){
                    while(TRUE){
                      bound = which.max(x[,j])
                      if (y[bound]==1) {x = matrix(x[-bound,],ncol=d); y = y[-bound]; weights = weights[-bound]}
                      else break
                    }
                }
                if((shape[j]=="ccv")||(shape[j]=="ccvde")||(shape[j]=="de")){
                    while(TRUE){
                      bound = which.max(x[,j])
                      if (y[bound]==0) {x = matrix(x[-bound,],ncol=d); y = y[-bound]; weights = weights[-bound]}
                      else break
                    }
                }
                if((shape[j]=="cvx")||(shape[j]=="cvxde")||(shape[j]=="de")){
                    while(TRUE){
                      bound = which.min(x[,j])
                      if (y[bound]==1) {x = matrix(x[-bound,],ncol=d); y = y[-bound]; weights = weights[-bound]}
                      else break
                    }
                }
                if((shape[j]=="ccv")||(shape[j]=="ccvin")||(shape[j]=="in")){
                    while(TRUE){
                      bound = which.min(x[,j])
                      if (y[bound]==0) {x = matrix(x[-bound,],ncol=d); y = y[-bound]; weights = weights[-bound]}
                      else break
                    }
                }
            } else {# must be gaussian / gamma family so no need to do anything
            }
        }

        n=nrow(x)
	if(n+1<d) stop("Input error: too few observations")

  
        ## run scar via active set algorithm

        # Setup basis functions
        xbasis=matrix(0, nrow=n, ncol=d)
        xbasisindex=matrix(0, nrow=n, ncol=d)
        xdiff=matrix(0, nrow=n-1, ncol=d)
        deriv=matrix(0, nrow=n, ncol=d)
        NumberOfBasis=rep(0,ncol=d)
        for (j in 1:d){
                xsort = sort.int(x[,j],index.return=TRUE)
                xbasis[,j]=xsort$x
                xbasisindex[,j]=xsort$ix
                xdiff[,j]=c(diff(xbasis[,j]))
        } 

        # initialization
        # build the constant and boundary matrix which are always active
        inactive = 0
	xinactive<-matrix(1,nrow=n,ncol=1)
        Ninitialinactive = 1
	for (j in 1:d){
           if ((shape[j]=="l")||(shape[j]=="cvx")||(shape[j]=="ccv")) {
               inactive = c(inactive, (1+(j-1)*n))
               xtemp = x[,j]-xbasis[1,j]
               xinactive = cbind(xinactive,xtemp,deparse.level = 0)
               Ninitialinactive = Ninitialinactive + 1
           }
        }
  
	# initial run 
        obj = glm.fit(xinactive,y,family=family,control=list(epsilon=epsilon),weights=weights,intercept = FALSE)
        iter = 1

        while(TRUE){               
            # find the derivatives and direction under 9 different settings
            score = obj$residuals*obj$weights
            score[is.na(score)] = 0
            for (j in 1:d){
                if(shape[j]=="cvxde"){
                    colxdiff = c(0,xdiff[,j])
                    colscore = score[xbasisindex[,j]]
                    colscoresum = c(0,aggsum(colscore)[-n])
                    deriv[,j]= aggsum(colxdiff*colscoresum)
                } else if ((shape[j]=="cvx")||(shape[j]=="cvxin")) {
                    colxdiff = c(xdiff[,j],0)
                    colscore = score[xbasisindex[,j]]
                    colscoresum = c(aggsum(colscore,reverse=TRUE)[-1],0)
                    deriv[,j]= aggsum(colxdiff*colscoresum,reverse=TRUE)
                } else if(shape[j]=="ccvin"){
                    colxdiff = c(0,xdiff[,j])
                    colscore = score[xbasisindex[,j]]
                    colscoresum = c(0,aggsum(colscore)[-n])
                    deriv[,j]= -aggsum(colxdiff*colscoresum)
                } else if (shape[j]=="ccv"||shape[j]=="ccvde") {
                    colxdiff = c(xdiff[,j],0)
                    colscore = score[xbasisindex[,j]]
                    colscoresum = c(aggsum(colscore,reverse=TRUE)[-1],0)
                    deriv[,j]= -aggsum(colxdiff*colscoresum,reverse=TRUE)
                } else if (shape[j]=="in") {
                    colscore = score[xbasisindex[,j]]
                    deriv[,j]= c(aggsum(colscore,reverse=TRUE)[-1],0)
                } else if (shape[j]=="de") {
                    colscore = score[xbasisindex[,j]]
                    deriv[,j]= c(0,aggsum(colscore)[-n])
                } else {
                    # shape == "l"
                    deriv[,j] = 0
                }
            }
            maxindex = which.max(deriv)
            direction = c(maxindex-n*(as.integer((maxindex-1)/n)),as.integer((maxindex-1)/n)+1)

            # exit the loop if no good direction can be found
            if(deriv[maxindex]<epsilon) break
            # exit the loop if the direction added in is already there...
            if(length(unique(sort(c(inactive,maxindex))))!=length(inactive)+1) break

           # add new element into inactive set
            inactive = c(inactive,maxindex)
            if (shape[direction[2]]=="cvxde"){
                xtemp = xbasis[direction[1],direction[2]] - x[,direction[2]]
                xinactive = cbind(xinactive,xtemp* as.integer(xtemp>0),deparse.level = 0)                 
            } else if ((shape[direction[2]]=="cvxin")||(shape[direction[2]]=="cvx")){
                xtemp = x[,direction[2]] - xbasis[direction[1],direction[2]]
                xinactive = cbind(xinactive,xtemp* as.integer(xtemp>0),deparse.level = 0)
            } else if ((shape[direction[2]]=="ccvin")){
                xtemp = x[,direction[2]] - xbasis[direction[1],direction[2]]
                xinactive = cbind(xinactive,xtemp* as.integer(xtemp<0),deparse.level = 0)
            } else if ((shape[direction[2]]=="ccvde")||(shape[direction[2]]=="ccv")){
                xtemp = xbasis[direction[1],direction[2]] - x[,direction[2]]
                xinactive = cbind(xinactive,xtemp* as.integer(xtemp<0),deparse.level = 0)
            } else if (shape[direction[2]]=="in"){
                xtemp = x[,direction[2]] - xbasis[direction[1],direction[2]]
                xinactive = cbind(xinactive,as.integer(xtemp>0),deparse.level = 0)
            } else if (shape[direction[2]]=="de"){
                xtemp = xbasis[direction[1],direction[2]] - x[,direction[2]]
                xinactive = cbind(xinactive,as.integer(xtemp>0),deparse.level = 0)
            } else {
                # must be linear "l"
                break
            }

     	    # try a initial run with new elements added
            wtsstart = c(obj$coefficients,0)
            objnew = glm.fit(xinactive,y,family=family,control=list(epsilon=epsilon),start=wtsstart,weights=weights,intercept = FALSE)

            # make sure all the inactive coefficients > 0
            while(TRUE){
                wtsnew = objnew$coefficients
                ninact = length(wtsnew)
                if(ninact <= Ninitialinactive) break else
                if (prod(wtsnew[(Ninitialinactive+1):ninact]>0)==1) break else 
                {
                    wtdenom = wtsstart - wtsnew
                    wtdenom = wtdenom * (abs(wtdenom) > epsilon) - epsilon * (wtdenom >= -epsilon) * ( wtdenom <=0) + epsilon * (wtdenom <= epsilon) * (wtdenom >0)
                    ratio = wtsstart/wtdenom
                    ratio = (ratio > 0)*(wtsnew <= 0)*ratio + 1*(wtsnew>0) + 1*(ratio<=0) 
                    ratio = min(ratio[(Ninitialinactive+1):ninact])
                    wtsstart = wtsstart * (1-ratio) + ratio*wtsnew
                    retain  = c(rep(TRUE,Ninitialinactive),(wtsstart>epsilon)[(Ninitialinactive+1):ninact])
                    xinactive = xinactive[,retain]
                    wtsstart = wtsstart[retain]
                    inactive = inactive[retain]
                    # check whether any element has been dropped
                    if(length(wtsstart)==ninact) break else
                    objnew = glm.fit(xinactive,y,family=family,control=list(epsilon=epsilon),start=wtsstart,weights=weights,intercept = FALSE)
                }
            }
            
           if (obj$deviance - objnew$deviance < epsilon){obj = objnew; break} else    
           {
               obj = objnew
               iter = iter + 1
           }
        }

        # output
        componentfit = matrix(0,ncol=d,nrow=n)
        constant = obj$coefficients[1]
        for(j in 1:d){
           comwts = matrix((inactive <= j*n)*(inactive > (j-1)*n)*obj$coefficients,ncol=1)
           componentfit[,j]=xinactive %*% comwts 
 	   if (min(x[,j]) <=0 && max(x[,j]) >=0) {
		constantj = approx(x[,j],componentfit[,j],0,rule=2)$y
		componentfit[,j] = componentfit[,j] - constantj
		constant = constant + constantj
  	   } else {
		constantj = approx(x[,j],componentfit[,j],mean(componentfit[,j]),rule=2)$y
		componentfit[,j] = componentfit[,j] - constantj
		constant = constant + constantj
   	   }   
        }
        deviance = obj$deviance
        nulldeviance = obj$nulldeviance
	result = list(x=x,y=y,shape=shape, weights=weights, family=family, componentfit=componentfit, constant=constant, deviance=deviance, nulldeviance=nulldeviance, iter=iter)
	class(result) <- "scar"
        return(result)
}
