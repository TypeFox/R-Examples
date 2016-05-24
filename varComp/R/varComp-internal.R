refugeEnvironment=new.env()
evalq(
	alt.fit <- control <- denom <- diag.1.n <- eigK <- k <- LI <- LIkLI <- LIy <- n <- negGrad <- negHess <- nK <- null <- null.fit <- nums <- nums2 <- numsPart <- tau.idx <- test <- tr1 <- tr2 <- y <- NULL
  , refugeEnvironment)
evalq(
	infoFunc <- function()NULL
  , refugeEnvironment)
  


V=evalq(function(tau){ ## diag.1.n, k
diag.1.n + Reduce('+', mapply('*', tau, k, SIMPLIFY=FALSE))
}, refugeEnvironment)

{
#  E=function(tau){
#    eigen(V(tau))
#  }
#  VIhalf=function(tau){
#    e=E(tau)
#    ans=tcrossprod( sweep( e$vec, 2, 1/sqrt(sqrt(e$val)), '*'))
#    attr(ans, 'values')=e$val
#    ans
#  }
}

updateLI=evalq(function(tau){
	if(isTRUE( attr(LI, 'tau')==tau )) return(LI)
	if(nK==1L){ 
		tmp= crossprod(eigK$tvec, 1/sqrt(eigK$val*tau+1) * eigK$tvec)
	}else
		tmp = t(backsolve(chol(V(tau)), diag.1.n))
	attr(tmp,'tau')=tau
	LI<<- tmp
}, refugeEnvironment)

updateLIkLI=evalq(function(){ #LI, k
	if(isTRUE( attr(LIkLI, 'tau')==attr(LI, 'tau'))) return(LIkLI)
	tmp=LIkLI
	for(j in 1:nK) tmp[[j]] = tcrossprod(LI%*%k[[j]], LI)
	attr(tmp, 'tau')=attr(LI, 'tau')
	LIkLI <<- tmp
}, refugeEnvironment)

# updateLIkLI=function(){ #LI, k
	# for(j in 1:nK) LIkLI[[j]]<<- tcrossprod(LI%*%k[[j]], LI)
# }

updateLIy=evalq(function(){
	if(isTRUE( attr(LIy, 'tau')==attr(LI, 'tau'))) return(LIy)
	LIy<<- structure( LI%*%y, tau=attr(LI, 'tau'))
}, refugeEnvironment)

PREML=evalq(function(){ # LIy, LI, n
	drop(
		if(nK==1L){
			.5*(
				-n*log(crossprod(LIy)) - sum(log(eigK$val*attr(LI,'tau')+1)) -n -n*log(2*pi)-n*log(n)
			)
		}else
			.5*(
				-n*log(crossprod(LIy)) +sum(log(diag(LI)))*2 -n -n*log(2*pi)-n*log(n)
			)
	)
}, refugeEnvironment)

preprocPREML=evalq(function(tau){
	updateLI(tau)
	updateLIy(); 
}  , refugeEnvironment)

obj=evalq(function(tau){
	updateLI(tau)
	updateLIy()
	PREML()
}, refugeEnvironment)

obj2=evalq(function(ltau)obj(exp(ltau)), refugeEnvironment)

updateNumsPart=evalq(function(){ # LIkLI, LIy
	if(isTRUE( attr(numsPart, 'tau') == attr(LIkLI, 'tau') )) return( numsPart )
	numsPart<<- sapply(LIkLI, '%*%', LIy)
}, refugeEnvironment)

updateNums=evalq(function(){ # LIy, numsPart
	if(isTRUE( attr(nums, 'tau') == attr(numsPart, 'tau'))) return( nums )
	nums<<-drop(crossprod(LIy, numsPart))
}, refugeEnvironment)

updateDenom=evalq(function(){ # LIy
	if(isTRUE( attr(denom, 'tau') == attr(LIy, 'tau'))) return( denom )
	denom<<-drop(crossprod(LIy))
}, refugeEnvironment)

updateTr1=evalq(function(){ # LIkLI
	if(isTRUE( attr(tr1, 'tau') == attr(LIkLI, 'tau'))) return( tr1 )
	tr1<<-sapply(LIkLI, function(z)sum(diag(z)))
}, refugeEnvironment)

updateTr2=evalq(function(){ # LIkLI
	if(isTRUE( attr(tr2, 'tau') == attr(LIkLI, 'tau'))) return( tr2 )
	for(i in 1:nK) {
		for(j in i:nK){
			kij=LIkLI[[i]]%*%LIkLI[[j]]
			tr2[i,j]<<-tr2[j,i]<<-sum(diag(kij))
		}
	}
}, refugeEnvironment)

updateNums2=evalq(function(){ # numsPart
	if(isTRUE( attr(nums2, 'tau') == attr(numsPart, 'tau'))) return( nums2 )
	for(i in 1:nK) {
		for(j in i:nK){
			nums2[i,j]<<-nums2[j,i]<<-2*crossprod(numsPart[,i], numsPart[,j])
		}
	}
}, refugeEnvironment)

score=evalq(function(){ # LIy, LIkLI, n, denoms, nums
	.5*(n*nums/denom - tr1)
}, refugeEnvironment)

updateNegGrad=evalq(function(){
	negGrad <<- -score()
}, refugeEnvironment)

gradFunc=evalq(function(tau){ # this agrees with numerical gradient :)
	updateLI(tau)
	updateLIkLI()
	updateLIy(); 
	updateNumsPart();   
	updateNums(); updateDenom(); updateTr1()
	score()
}, refugeEnvironment)

OI=evalq(function(){ # nums2, nums, denom, n, tr2
	.5*(
		(nums2-outer(drop(nums), drop(nums))/denom ) /denom*n - tr2
	)
}, refugeEnvironment)

hess=evalq(function(tau){  # this agrees with numerical hessian :)
	updateLI(tau)
	updateLIkLI(); updateLIy(); updateNumsPart(); 
	updateNums(); updateDenom()
	updateNums2()
	updateTr2()
	-OI()
}, refugeEnvironment)

EI=evalq(function(){ # tr2, tr1, n
	.5/(n+2)*(n*tr2-outer(tr1,tr1))
}, refugeEnvironment)

AOI=evalq(function(){ # nums2, nums, denom
	.5*(
		(nums2/2-outer(nums, nums)/denom )/denom*n
	)
}, refugeEnvironment)

AEI=evalq(function(){ # nums2, tr1, denom
	.5/(n+2)*(
		n*n*nums2/denom/2 - outer(tr1, tr1)
	)
}, refugeEnvironment)

WAI=evalq(function(){ # nums2, denom, nums, tr1, n
	.5*.5/(n+1)*(
		n*n*nums2/denom - n*n*outer(nums, nums)/denom/denom - outer(tr1, tr1)
	)
}, refugeEnvironment)

updateNegHess=evalq(function(){
	negHess <<- infoFunc()
}, refugeEnvironment)

ToDoList=c(
"update/add1/drop1",
"anova interface with add1/drop1",
"anova / fixef : printing significance code",
"logLik",
"fitted",
"residual",
"BLUP",
"kernel prediction",
"lme conversion",
"bartlett"
)