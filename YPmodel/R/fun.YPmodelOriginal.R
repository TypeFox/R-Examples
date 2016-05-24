fun.YPmodelOriginal <-
function(Data, Parameters)
########################################################
#fun.YPmodelOriginal(Data,Parameters)
#######################################################
# version 0.1
# May 19, 2012
# Junlong Sun
# return [best,r]
#######################################################
# May 19, 2012 - v0.1 Create
# Jul 19, 2012 - v0.2 Add output
#######################################################
{
#-----------------------------------------------------------------#
## loading data
#-----------------------------------------------------------------#
	nm <- Parameters$nm 
	maxIteration1<- Parameters$maxIter1
	maxIteration2<- Parameters$maxIter2


Z <- Data$Z
Delta <- Data$Delta
n <- Data$length

nb <- array(c(0,0),c(1,2))
ob <- nb + 0.1
bc <- t(nb)
jh <- 0.01

##return [po,s]
#data2 <- fun.hcvxitr(bc,jh,Z,Delta,n)
#po <- data2$po
#s <- data2$s

data2 <- fun.hcvxitr1(bc,Data,Parameters)
po <- data2$po
s <- data2$s

of <- s[1]
nm <- 100

numIteration <- 0

#-----------------------------------------------------------------#
while ((sum(abs(round(1000*ob) - round(1000*nb)))>0) & (numIteration < maxIteration1) ) {
    
    numIteration <- numIteration +1
    t <- 1
    bc <- t(nb) + t*po
    idb <- sum((bc > log(nm)) + (bc < -log(nm)))
    
	numIterationInner<- 0
    
   #-----------------------------------------------------------------#
   while ((idb>0)&(numIterationInner<maxIteration2)){
        numIterationInner<- numIterationInner+1
        t <- t*0.5
        bc <- t(nb)+ t*po
        idb <- sum(( bc > log(nm)) + (bc < -log(nm)))
    }
   #-----------------------------------------------------------------#

    m1 <- 1
    b <- bc
    
   ##return [s,ru]
   #data4 <- fun.ntitr(b,Z,Delta,n,m1)
   #s <- data4$s
   #ru <- data4$ru

	data4 <- fun.oldp2(b,1,Data)
	s <- data4$s	
	ru <- data4$ru

    if (of==0){
        break
    } else if (abs((s/of-1)*1.e+8)<1){
        break
    }
    
    ob <- nb
    nb <- t(bc)
    of <- s

   ##return [po,s]
   #data3 <- fun.hcvxitr(bc,jh,Z,Delta,n)
   #po <- data3$po
   #s <- data3$s

   data3 <- fun.hcvxitr1(bc,Data,Parameters)
   po <- data3$po
   s <- data3$s

}
#-----------------------------------------------------------------#

b <- cbind(t(nb), t(ob))
m <- 2

##return [s,ru]
#data5 <- fun.ntitr(b,Z,Delta,n,m)
#s <- data5$s
#ru <- data5$ru

data5 <- fun.oldp2(b,m,Data)
s <- data5$s
ru <- data5$ru
        
ots <- min(s)
jt <- which.min(s)
best <- t(b[,jt])
r <- ru[,jt]

#-----------------------------------------------------------------#
## Output Resuts 
#-----------------------------------------------------------------#
    output<- list(best=best,r=r)
    return(output)
#-----------------------------------------------------------------#

}
