hiv.mortmod <-
function(prev=NULL, e0=NULL, child.mort=NULL, adult.mort=NULL, model=1, region=1, sex=1, lt=FALSE, nax=NULL, opt=TRUE){
#data(HIV-MLTs-obs)
# e0
if(model==1){
	mx.out <- mortmod.e0(e0=e0, prev=prev, region=region, sex=sex, opt=opt)
}
# 5q0
if(model==2){
	mx.out <- mortmod.5q0(child.mort=child.mort, prev=prev, region=region, sex=sex, opt=opt)
}
# 5q0 and 45q15
if(model==3){
	mx.out <- mortmod.45q15(child.mort=child.mort, adult.mort=adult.mort, prev=prev, region=region, sex=sex, opt=opt)
}
age <- c(0,1,seq(5,100,5))
nmx <- mx.out
## make a life table
if(lt==TRUE){
	if(is.null(nax)){

		nax <- rep(2.5,length(age))

		if (sex=="male") { # TRUE for male, FALSE for female
			if (nmx[1] >= 0.107) {
				nax[1] <- 0.33
				nax[2] <- 1.352
			} else {
				nax[1] <- 0.045+2.684*nmx[1]
				nax[2] <- 1.651-2.816*nmx[1]
				}
			} 
		if (sex=="female") {
			if (nmx[1] >= 0.107) {
				nax[1] <- 0.35
				nax[2] <- 1.361
			} else {
				nax[1] <- 0.053+2.8*nmx[1]
				nax[2] <- 1.522-1.518*nmx[1]
				}
			}
	} else {nax=nax}
    
    n <- c(diff(age), 999)
    nqx <- (n*nmx)/(1+(n-nax)*nmx)
    nqx <- c(nqx[-(length(nqx))], 1)
    for(i in 1:length(nqx)){
    	if(nqx[i] > 1) nqx[i] <- 1
    	}
    nage <- length(age)

    #nqx <- round(nqx, 4)

    npx <- 1 - nqx
    max.age <- min(which(npx==0))
    l0=100000
    lx <- round(cumprod(c(l0, npx)))
    #lx <- (cumprod(c(l0, npx)))
    ndx <- -diff(lx)
    lxpn <- lx[-1]
    nLx <- n * lxpn + ndx * nax
    lx <- lx[1:length(age)]
    nLx[max.age] <- lx[max.age]/nmx[max.age]
    Tx <- rev(cumsum(rev(nLx)))
    ex <- Tx/lx
    lt <- cbind(Age = age, nax = nax, nmx = nmx, nqx = nqx, npx = npx, ndx = ndx, lx = lx, nLx = round(nLx), Tx = round(Tx), 
        ex = round(ex, 2))
     lt <- lt[lt[,6]!=0,]
     e0 <- lt[1,10]
     lt.45q15 <- 1-(lx[14]/lx[5])
     lt.5q0 <- 1-(lx[3]/lx[1])
     return(list(e0=e0, lt.5q0=lt.5q0, lt.45q15=lt.45q15, lt=lt))
} else {
return(mx.out)}
}
