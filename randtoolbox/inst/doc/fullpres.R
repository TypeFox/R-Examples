### R code from vignette source 'fullpres.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: wh.predict
###################################################
wh.predict <- function(x)
{
    M1 <- 30269
    M2 <- 30307
    M3 <- 30323
    y <- round(M1*M2*M3*x)
    s1 <- y %% M1
    s2 <- y %% M2
    s3 <- y %% M3
    s1 <- (171*26478*s1) %% M1
    s2 <- (172*26070*s2) %% M2
    s3 <- (170*8037*s3) %% M3
    (s1/M1 + s2/M2 + s3/M3) %% 1
}

RNGkind("Wichmann-Hill")
xnew <- runif(1)
maxerr <- 0
for (i in 1:1000) {
    xold <- xnew
    xnew <- runif(1)
	err <- abs(wh.predict(xold) - xnew)
    maxerr <- max(err, maxerr)
}
print(maxerr)


###################################################
### code chunk number 2: congrurand (eval = FALSE)
###################################################
## congruRand(10)


###################################################
### code chunk number 3: congrurand2
###################################################
library(randtoolbox)
options(width = 40)
congruRand(10)


###################################################
### code chunk number 4: congrurandseed1 (eval = FALSE)
###################################################
## setSeed(1)
## congruRand(10)


###################################################
### code chunk number 5: congrurandseed2
###################################################
options( width =40)
setSeed(1)
congruRand(10)


###################################################
### code chunk number 6: congrurandseed3 (eval = FALSE)
###################################################
## setSeed(1)
## congruRand(10, echo=TRUE)


###################################################
### code chunk number 7: congrurandseed4
###################################################
options( width =40)
setSeed(1)
congruRand(10, echo=TRUE)


###################################################
### code chunk number 8: congrurandseed5 (eval = FALSE)
###################################################
## setSeed(1614852353)
## congruRand(5, echo=TRUE)


###################################################
### code chunk number 9: congrurandseed6
###################################################
options( width =40)
setSeed(1614852353)
congruRand(5, echo=TRUE)


###################################################
### code chunk number 10: congrurandseed5 (eval = FALSE)
###################################################
## setSeed(12)
## congruRand(5, mod = 2^8, mult = 25, incr = 16, echo= TRUE)


###################################################
### code chunk number 11: congrurandseed6
###################################################
options( width =30)
setSeed(12)
congruRand(5, mod = 2^8, mult = 25, incr = 16, echo= TRUE)


###################################################
### code chunk number 12: sfmt (eval = FALSE)
###################################################
## SFMT(10)
## SFMT(5, 2) #bi dimensional variates


###################################################
### code chunk number 13: sfmt2
###################################################
options( width =40)
SFMT(10)
SFMT(5, 2)


###################################################
### code chunk number 14: sfmt3 (eval = FALSE)
###################################################
## SFMT(10, mexp = 607)


###################################################
### code chunk number 15: sfmt4
###################################################
options( width =40)
SFMT(10, mexp = 607)


###################################################
### code chunk number 16: halton (eval = FALSE)
###################################################
## halton(10)
## halton(10, 2)


###################################################
### code chunk number 17: halton2
###################################################
options( width =40)
halton(10)
halton(10, 2)


###################################################
### code chunk number 18: halton3 (eval = FALSE)
###################################################
## halton(5)
## halton(5, init=FALSE)


###################################################
### code chunk number 19: halton4
###################################################
options( width =40)
halton(5)
halton(5, init=FALSE)


###################################################
### code chunk number 20: sobol (eval = FALSE)
###################################################
## sobol(10)
## sobol(10, scramb=3)


###################################################
### code chunk number 21: sobol2
###################################################
options( width =40)
sobol(10)
sobol(10, scramb=3)


###################################################
### code chunk number 22: unitsquare1 (eval = FALSE)
###################################################
## par(mfrow = c(2,1))
## plot(sobol(1000, 2))
## plot(sobol(10^3, 2, scram=1))


###################################################
### code chunk number 23: unitsquare2
###################################################
par(mfrow = c(2,1))
plot(sobol(1000, 2), xlab ="u", ylab="v", main="Sobol (no scrambling)")
plot(sobol(10^3, 2, scram=1), xlab ="u", ylab="v", main="Sobol (Owen)")


###################################################
### code chunk number 24: torus (eval = FALSE)
###################################################
## torus(10)


###################################################
### code chunk number 25: torus2
###################################################
options( width =40)
torus(10)


###################################################
### code chunk number 26: torus3 (eval = FALSE)
###################################################
## torus(5, use =TRUE)


###################################################
### code chunk number 27: torus4
###################################################
options( width =40)
torus(5, use =TRUE)


###################################################
### code chunk number 28: torus5 (eval = FALSE)
###################################################
## torus(5, p =7)


###################################################
### code chunk number 29: torus6
###################################################
options( width =40)
torus(5, p =7)


###################################################
### code chunk number 30: torus7 (eval = FALSE)
###################################################
## torus(5,  mixed =TRUE)


###################################################
### code chunk number 31: torus8
###################################################
options( width =40)
torus(5,  mixed =TRUE)


###################################################
### code chunk number 32: torus9 (eval = FALSE)
###################################################
## par(mfrow = c(2,1))
## acf(torus(10^5))
## acf(torus(10^5, mix=TRUE))


###################################################
### code chunk number 33: torusacf
###################################################
par(mfrow = c(2,1))
acf(torus(10^5))
acf(torus(10^5, mix=TRUE))


###################################################
### code chunk number 34: unitsquare3 (eval = FALSE)
###################################################
## par(mfrow = c(2,1))
## plot(SFMT(1000, 2))
## plot(torus(10^3, 2))


###################################################
### code chunk number 35: unitsquare4
###################################################
par(mfrow = c(2,1))
plot(SFMT(1000, 2), xlab ="u", ylab="v", main="SFMT")
plot(torus(1000, 2), xlab ="u", ylab="v", main="Torus")


###################################################
### code chunk number 36: unitsquare5 (eval = FALSE)
###################################################
## par(mfrow = c(2,1))
## plot(WELL(1000, 2))
## plot(sobol(10^3, 2, scram=2))


###################################################
### code chunk number 37: unitsquare6
###################################################
par(mfrow = c(2,1))
plot(WELL(1000, 2), xlab ="u", ylab="v", main="WELL 512a")
plot(sobol(10^3, 2, scram=2), xlab ="u", ylab="v", main="Sobol (Faure-Tezuka)")


###################################################
### code chunk number 38: integralcos (eval = FALSE)
###################################################
## I25 <- -1356914
## nb <- c(1200, 14500, 214000)
## ans <- NULL
## for(i in 1:3)
## {
## 	tij <- sobol(nb[i], dim=25, scramb=2, norm=TRUE )
## 	Icos <- mean(cos(sqrt( apply( tij^2/2, 1, sum ) ))) * pi^(25/2)
## 	ans <- rbind(ans, c(n=nb[i], I25=Icos, Delta=(Icos-I25)/I25 ))
## }
## data.frame(ans)


###################################################
### code chunk number 39: integralcos2
###################################################
options( width =40)
I25 <- -1356914
nb <- c(1200, 14500, 214000)
ans <- NULL
for(i in 1:3)
{
	tij <- sobol(nb[i], dim=25, scramb=2, norm=TRUE )
	Icos <- mean(cos(sqrt( apply( tij^2/2, 1, sum ) ))) * pi^(25/2)
	ans <- rbind(ans, c(n=nb[i], I25=Icos, Delta=(Icos-I25)/I25 ))
}
data.frame(ans)


###################################################
### code chunk number 40: hist (eval = FALSE)
###################################################
## par(mfrow = c(2,1))
## hist(SFMT(10^3), 100)
## hist(torus(10^3), 100)


###################################################
### code chunk number 41: hist
###################################################
par(mfrow = c(2,1))
hist(SFMT(10^3), 100)
hist(torus(10^3), 100)


###################################################
### code chunk number 42: gap1 (eval = FALSE)
###################################################
## gap.test(runif(1000))


###################################################
### code chunk number 43: gap2
###################################################
options( width =40)
res <- gap.test(runif(1000), echo=FALSE)
	stat <- res$statistic
	pvalue <- res$p.value
	df <- res$parameter
	obsnum <- res$observed
	expnum <- res$expected
	
	cat("\n\t Gap test\n")
	cat("\nchisq stat = ", stat, ", df = ",df, "\n, p-value = ", pvalue, "\n", sep="")
	cat("\n(sample size : ",1000,")\n\n", sep="")
       cat("length observed freq theoretical freq\n")
       for(i in 1:(df+1))
                cat(i,"\t", obsnum[i],"\t", expnum[i],"\n")


###################################################
### code chunk number 44: gap3 (eval = FALSE)
###################################################
## gap.test(SFMT(1000), 1/3, 2/3)


###################################################
### code chunk number 45: order1 (eval = FALSE)
###################################################
## order.test(runif(4000), d=4)


###################################################
### code chunk number 46: order2
###################################################
options( width =40)
res <- order.test(runif(4000), d=4, echo=FALSE)
	stat <- res$statistic
	pvalue <- res$p.value
	df <- res$parameter
	obsnum <- res$observed
	expnum <- res$expected
	
	cat("\n\t Order test\n")
         cat("\nchisq stat = ", stat, ", df = ",df, "\n, p-value = ", pvalue, "\n", sep="")
         cat("\n (sample size : ", 1000,")\n\n", sep="")
         cat("observed number ",obsnum[1:6],"\n",obsnum[7:18],"\n", obsnum[19:24],"\n")
         cat("expected number ",expnum,"\n")  


###################################################
### code chunk number 47: freq (eval = FALSE)
###################################################
## freq.test(runif(1000), 1:4)


###################################################
### code chunk number 48: freq2
###################################################
options( width =40)
res <- freq.test(runif(1000), 1:4, echo=FALSE)
	stat <- res$statistic
	pvalue <- res$p.value
	df <- res$parameter
	obsnum <- res$observed
	expnum <- res$expected
	         
         cat("\n\t Frequency test\n")
         cat("\nchisq stat = ", stat, ", df = ",df, "\n, p-value = ", pvalue, "\n", sep="")
         cat("\n (sample size : ",1000,")\n\n", sep="")
         cat("observed number ",obsnum,"\n")
         cat("expected number ",expnum,"\n")    


###################################################
### code chunk number 49: serial (eval = FALSE)
###################################################
## serial.test(runif(3000), 3)


###################################################
### code chunk number 50: serial2
###################################################
options( width =40)
res <- serial.test(runif(3000), 3, echo=FALSE)
	stat <- res$statistic
	pvalue <- res$p.value
	df <- res$parameter
	obsnum <- res$observed
	expnum <- res$expected
	         
         cat("\n\t Serial test\n")
         cat("\nchisq stat = ", stat, ", df = ",df, "\n, p-value = ", pvalue, "\n", sep="")
         cat("\n (sample size : ",3000,")\n\n", sep="")
         cat("observed number ",obsnum[1:4],"\n", obsnum[5:9])
         cat("expected number ",expnum,"\n")    


###################################################
### code chunk number 51: coll1 (eval = FALSE)
###################################################
## coll.test(runif, 2^7, 2^10, 1)


###################################################
### code chunk number 52: coll2
###################################################
options( width =40)
res <- coll.test(runif, 2^7, 2^10, 1, echo=FALSE)
	stat <- res$statistic
	pvalue <- res$p.value
	df <- res$parameter
	obsnum <- res$observed
	expnum <- res$expected
	         
         cat("\n\t Collision test\n")
         cat("\nchisq stat = ", stat, ", df = ", df, "\n, p-value = ", pvalue, "\n", sep="")
         cat("\n exact distribution \n(sample number : ", 1000,"/sample size : ", 128,"\n / cell number : ", 1024,")\n\n", sep="")
         cat("collision observed expected\n", "number    count    count\n", sep="")
         for(i in 1:(df + 1) )
              cat(" ", i,"    ", obsnum[i],"    ", expnum[i],"\n")


###################################################
### code chunk number 53: coll3 (eval = FALSE)
###################################################
## coll.test(congruRand, 2^8, 2^14, 1)


###################################################
### code chunk number 54: coll4
###################################################
options( width =40)
res <- coll.test(congruRand, 2^8, 2^14, 1, echo=FALSE)
	stat <- res$statistic
	pvalue <- res$p.value
	df <- res$parameter
	obsnum <- res$observed
	expnum <- res$expected
	         
         cat("\n\t Collision test\n")
         cat("\nchisq stat = ", stat, ", df = ", df, "\n, p-value = ", pvalue, "\n", sep="")
         cat("\n Poisson approximation \n(sample number : ", 1000,"/sample size : ", 256,"\n / cell number : ", 16384,")\n\n", sep="")
         cat("collision observed expected\n", "number    count    count\n", sep="")
         for(i in 1:(df + 1) )
              cat(" ", i-1,"    ", obsnum[i],"    ", expnum[i],"\n")


###################################################
### code chunk number 55: poker (eval = FALSE)
###################################################
## poker.test(SFMT(10000))


###################################################
### code chunk number 56: poker2
###################################################
options( width =40)
res <- poker.test(SFMT(10000), echo=FALSE)
	stat <- res$statistic
	pvalue <- res$p.value
	df <- res$parameter
	obsnum <- res$observed
	expnum <- res$expected
	         
            cat("\n\t Poker test\n")
            cat("\nchisq stat = ", stat, ", df = ", df, "\n, p-value = ", pvalue, "\n", sep="")
            cat("\n (sample size : ", 10000,")\n\n", sep="")
            cat("observed number ", obsnum,"\n")
            cat("expected number ", expnum,"\n") 


###################################################
### code chunk number 57: wh1
###################################################
wh.predict <- function(x)
{
    M1 <- 30269
    M2 <- 30307
    M3 <- 30323
    y <- round(M1*M2*M3*x)
    s1 <- y %% M1
    s2 <- y %% M2
    s3 <- y %% M3
    s1 <- (171*26478*s1) %% M1
    s2 <- (172*26070*s2) %% M2
    s3 <- (170*8037*s3) %% M3
    (s1/M1 + s2/M2 + s3/M3) %% 1
}

RNGkind("Wichmann-Hill")
xnew <- runif(1)
err <- 0
for (i in 1:1000) {
    xold <- xnew
    xnew <- runif(1)
    err <- max(err, abs(wh.predict(xold) - xnew))
}
print(err)


