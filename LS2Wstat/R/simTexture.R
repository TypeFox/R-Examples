simTexture <-function (n = 256, sd = 1, K = 150, imtype = "S1", ...){
 
opt.args <- c("order", "fn", "start", "end", "a", "prop", 
        "sd2", "prop2", "pos1", "pos2","nu","rho","type")

# set default optional values:  

order <- 5
fn <- lincurve
start <- 1
end <- 2
a <- 1
prop <- 0.25
sd2 <- 1
prop2 <- 0.25
pos1 <- c(n * (0.5 - prop)/2, n * (0.5 - prop)/2)
pos2 <- c(n/4, 3 * n/4)
nu <- 1
rho <- 0.9
type <- "queen"

argl <- as.list(match.call())
oargind <- match(names(argl), opt.args)
argind <- which(!is.na(oargind))

if (length(argind) > 0) {
        called.args <- names(argl)[argind]
        for (i in 1:length(called.args)) {
        	assign(called.args[i], eval(argl[[argind[i]]]))
        }
}

images <- list()
switch(imtype, 
	S1 = {
        	for (i in 1:K) {
            		images[[i]] <- matrix(rnorm(n^2, mean = 0, sd = sd), nrow = n, ncol = n, byrow = TRUE)
        	}
	},
	S2 = {
	
		sma.transform <- function(x,w,rho){

		wx <- lag.listw(nb2listw(w),x)
		smax <- rho * wx + x
		smax

		}

		wm<-cell2nb(n,n,type=type)

		for(i in 1:K){
			x<-rnorm(n^2,sd=sd)

			y<-sma.transform(x,wm,rho=rho)

			y<-matrix(y,nrow=n)

			images[[i]]<-y
		}

 	},
	S3={
		s<-seq(0,1,length=n)

		for(i in 1:K){
			x<-GaussRF(s,s,model="matern",param=c(mean=0,variance=sd^2,nugget=0,scale=.5,nu=nu),grid=T)
			images[[i]]<-x
		}
	}, 
	S4={
		data<-grf(n^2,grid="reg",nugget=0.1,cov.model="exponential", cov.pars=c(sd^2,2))$data

		data<-array(data,dim=c(n,n,K))

		for(i in 1:K){
			images[[i]]<-data[,,i]
		}
	},
	S5 = {
        	for (i in 1:K) {
            		images[[i]] <- Haar2MA.diag(n, sd, ...)
        	}
    	},
	NS1 = {
        	for (i in 1:K) {
            		a <- matrix(rnorm((n^2)/2, mean = 0, sd = 1), nrow = n, ncol = n/2, byrow = TRUE)
  		        b <- matrix(rnorm((n^2)/2, mean = 0, sd = sd), nrow = n, ncol = n/2, byrow = TRUE)
            		images[[i]] <- cbind(a, b)
        	}
    	}, 
	NS2 = {
		s<-seq(0,1,length=n)

		for(i in 1:K){
			b1 <-GaussRF(s,s,model="matern",param=c(mean=0,variance=sd^2,nugget=0,scale=.5,nu=nu),grid=T)
            		b <- matrix(rnorm(n^2, mean = 0, sd = sd), nrow = n, ncol = n, byrow = TRUE)
			x<-cropimage(cbind(b,b1),newsize=n)
			images[[i]]<-x
		}
	},
	NS3 = {
        	for (i in 1:K) {
            		images[[i]] <- HaarMontage(n, sd = sd, ...)
        	}
	}, 
	NS4 = {
        	x <- seq(0, 1, length = n)
        	for (i in 1:K) {
            		images[[i]] <- matrix(rnorm(n^2, 0, fn(x, start = start, end = end, a = a)), n, n)
        	}
    	}, 
	NS5 = {
        	inner <- list()
        	S1 <- matrix(rnorm(n^2, mean = 0, sd = 1), nrow = n, ncol = n, byrow = TRUE)
        	for (i in 1:K) {
            		inner[[i]] <- matrix(rnorm(n^2, mean = 0, sd = sd), nrow = n, ncol = n, byrow = TRUE)
        	}
        	images <- lapply(inner, mix2images, imageB = S1, prop = prop)
    	}, 
	NS6 = {
        	inner <- list()
        	S1 <- matrix(rnorm(n^2, mean = 0, sd = 1), nrow = n, ncol = n, byrow = TRUE)
        	for (i in 1:K) {
            		inner[[i]] <- matrix(rnorm(n^2, mean = 0, sd = sd), nrow = n, ncol = n, byrow = TRUE)
        	}
        	images <- lapply(inner, mix2images, imageB = S1, prop = prop, 
				pos = c(n * (0.5 - prop)/2, n * (0.5 - prop)/2))
    	}, 
	NS7 = {
        	inner <- list()
        	inner2 <- list()
        	S1 <- matrix(rnorm(n^2, mean = 0, sd = 1), nrow = n, ncol = n, byrow = TRUE)
        	for (i in 1:K) {
            		inner[[i]] <- matrix(rnorm(n^2, mean = 0, sd = sd), nrow = n, ncol = n, byrow = TRUE)
            		inner2[[i]] <- matrix(rnorm(n^2, mean = 0, sd = sd2), nrow = n, ncol = n, byrow = TRUE)
        	}
        	images <- lapply(inner, mix2images, imageB = S1, prop = prop, pos = pos1)
        	for (i in 1:K) {
            		images[[i]] <- mix2images(inner2[[i]], images[[i]], prop = prop2, pos = pos2)
        	}
    	})

return(images)

}
