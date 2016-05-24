### R code from vignette source 'tree-dust.Rnw'

###################################################
### code chunk number 1: tree-dust.Rnw:13-14 (eval = FALSE)
###################################################
## options(width=80)


###################################################
### code chunk number 2: tree-dust.Rnw:19-21 (eval = FALSE)
###################################################
## library(catdata)
## data(dust)


###################################################
### code chunk number 3: tree-dust.Rnw:26-27 (eval = FALSE)
###################################################
## library(rpart)


###################################################
### code chunk number 4: tree-dust.Rnw:32-37 (eval = FALSE)
###################################################
## tree1 <-rpart(as.factor(bronch) ~ years, data = dust, 	method = "class",
##               control = rpart.control(cp = 0.001,  parms=list(split='information'),
##                                       maxdepth = 4))
## plot(tree1, xpd=TRUE)
## text(tree1)


###################################################
### code chunk number 5: tree-dust.Rnw:42-70 (eval = FALSE)
###################################################
##  	pred <- predict(tree1)
##  	year<- dust$years
##  	year [dust$years<15.5] <- 1
##  	year [dust$years>15.5 & dust$years<36.5] <- 2
##  	year [dust$years>36.5 & dust$years<47.5] <- 3
##  	year [dust$years>47.5 & dust$years<50.5] <- 4 	
##  	year [dust$years>50.5] <- 5
##  	
## pre5 <- unique( pred[,2][year==5])
## pre4 <- unique( pred[,2][year==4])
## pre3 <- unique( pred[,2][year==3])
## pre2 <- unique( pred[,2][year==2])
## pre1 <- unique( pred[,2][year==1])
##  	      
## meanyear <- c()
## 
## for (i in min(dust$years):max(dust$years)){
## meanyear[i] <- sum(dust$bronch[dust$year==i])
## if(sum(dust$bronch[dust$year==i])!=0){
## meanyear[i] <- mean(dust$bronch[dust$year==i])
## }
## }
## 
## dust$means<- rep(2, nrow(dust)) 
##  
##  for (k in 1:nrow(dust)){
##  dust$means[k] <- meanyear[dust$years[k]]
##  }


###################################################
### code chunk number 6: tree-dust.Rnw:73-83 (eval = FALSE)
###################################################
##  plot(dust$years, dust$means, xlab="years",ylab="")
##  segments(x0=3,x1=15.5,y0=pre1)
##  segments(x0=15.5,x1=15.5,y0=pre1,y1=pre2)
##  segments(x0=15.5,x1=36.5,y0=pre2)
##  segments(x0=36.5,x1=36.5,y0=pre2,y1=pre3)
##  segments(x0=36.5,x1=47.5,y0=pre3)
##   segments(x0=47.5,x1=47.5,y0=pre3,y1=pre4)
##  segments(x0=47.5,x1=50.5,y0=pre4)
##   segments(x0=50.5,x1=50.5,y0=pre4,y1=pre5)
##  segments(x0=50.5,x1=66,y0=pre5)


###################################################
### code chunk number 7: tree-dust.Rnw:88-89 (eval = FALSE)
###################################################
## library(party)


###################################################
### code chunk number 8: tree-dust.Rnw:94-96 (eval = FALSE)
###################################################
## treeP1 <-ctree(as.factor(bronch) ~ years, data = dust)
## plot(treeP1)


###################################################
### code chunk number 9: tree-dust.Rnw:98-108 (eval = FALSE)
###################################################
##  	year<- dust$years
##  	year [dust$years<7.5] <- 1
##  	year [dust$years>7.5 & dust$years<15.5] <- 2
##  	year [dust$years>15.5 & dust$years<36.5] <- 3
##  	year [dust$years>36.5] <- 4
##  	
## pre4 <- mean(dust$bronch[year==4])
## pre3 <- mean(dust$bronch[year==3])
## pre2 <- mean(dust$bronch[year==2])
## pre1 <- mean(dust$bronch[year==1])


###################################################
### code chunk number 10: tree-dust.Rnw:111-120 (eval = FALSE)
###################################################
##  plot(dust$years, dust$means, xlab="years",ylab="")
##  segments(x0=3,x1=7.5,y0=pre1)
##  segments(x0=7.5,x1=7.5,y0=pre1,y1=pre2)
##  segments(x0=7.5,x1=15.5,y0=pre2)
##  segments(x0=15.5,x1=15.5,y0=pre2,y1=pre3)
##  segments(x0=15.5,x1=36.5,y0=pre3)
##   segments(x0=36.5,x1=36.5,y0=pre3,y1=pre4)
##  segments(x0=36.5,x1=66,y0=pre4)
## 


###################################################
### code chunk number 11: tree-dust.Rnw:126-128 (eval = FALSE)
###################################################
## treeP2 <-ctree(as.factor(bronch) ~ smoke + years + dust, data = dust)
## plot(treeP2)


###################################################
### code chunk number 12: tree-dust.Rnw:132-134 (eval = FALSE)
###################################################
## detach(package:rpart)
## detach(package:party)


