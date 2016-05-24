iwls_mh <-
function(curr.y,curr.X,curr.beta,iprior.var){

curr.LP<-as.vector(curr.X%*%matrix(curr.beta,ncol=1)) 		## current linear predictor

curr.w<-exp(curr.LP)  						## current elements of weight matrix
icurr.C<-iprior.var+crossprod(x=curr.X*curr.w,y=curr.X)  	## inverse proposal variance
#curr.C<-solve(icurr.C)
curr.C<-chol2inv(chol(icurr.C))						## proposal variance
curr.z<-curr.LP+(curr.y-curr.w)/curr.w 				## current value of working vector
curr.m<-as.vector(tcrossprod(x=curr.C,y=curr.X*curr.w)%*%matrix(curr.z,ncol=1)) ## proposal mean
prop.beta<-as.vector(rmvnorm(n=1,mean=curr.m,sigma=curr.C))	## proposal

prop.LP<-as.vector(curr.X%*%matrix(prop.beta,ncol=1))		## proposed linear predictor

prop.w<-exp(prop.LP)						## proposed elememts of weight matrix
iprop.C<-iprior.var+crossprod(x=curr.X*prop.w,y=curr.X)		## inverse current variance
#prop.C<-solve(iprop.C)						## current variance
prop.C<-chol2inv(chol(iprop.C))						## current variance
prop.z<-prop.LP+(curr.y-prop.w)/prop.w				## proposed working vector
prop.m<-as.vector(tcrossprod(x=prop.C,y=curr.X*prop.w)%*%matrix(prop.z,ncol=1))##current mean

top<-sum(curr.y*prop.LP)-sum(prop.w)-0.5*as.vector(matrix(prop.beta[-1],nrow=1)%*%iprior.var[-1,-1]%*%matrix(prop.beta[-1],ncol=1))+dmvnorm(x=curr.beta,mean=prop.m,sigma=prop.C,log=TRUE) ## log numerator of acceptance probability
bot<-sum(curr.y*curr.LP)-sum(curr.w)-0.5*as.vector(matrix(curr.beta[-1],nrow=1)%*%iprior.var[-1,-1]%*%matrix(curr.beta[-1],ncol=1))+dmvnorm(x=prop.beta,mean=curr.m,sigma=curr.C,log=TRUE) ## log denominatrot of acceptance probability

prob<-exp(top-bot)						## acceptance probability

if(prob>=runif(1)){						## accept or reject
new.beta<-prop.beta} else{
new.beta<-curr.beta}

new.beta}
