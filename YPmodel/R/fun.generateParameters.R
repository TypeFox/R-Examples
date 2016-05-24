fun.generateParameters <-
function(Data)
########################################################
#fun.generateParameters(y,d,z,n)
#Culculate [fb1,fb2]
#######################################################
# version 0.2
# May 6, 2012
# Junlong Sun
# return [fb1,fb2]
#######################################################
# May 10, 2012 - v0.1 Create
# May 19, 2012 - v0.1 change input and output
#######################################################
{
#-----------------------------------------------------------------#
## loading data
#-----------------------------------------------------------------#
	n <- Data$length
	y <- Data$X
	d <- Data$Delta
	z <- Data$Z

#-----------------------------------------------------------------#
## main function
#-----------------------------------------------------------------#

	temp <- matrix(c(1:n), nrow = n, ncol = 1)
	zd1 <- temp[z==0]
	zd2 <- temp[z==1]
	
	y1 <- y[zd1]
	y2 <- y[zd2]
	d1 <- d[zd1]
	d2 <- d[zd2]

	n1 <- length(y1)
	n2 <- length(y2)

	k1 <- matrix(c(n1:1), nrow = n1, ncol = 1 )
	k2 <- matrix(c(n2:1), nrow = n2, ncol = 1 )
	dla1 <- d1/k1
	dla2 <- d2/k2

	la2 <- fun.cumsum(dla2)
	s1 <- exp(-fun.cumsum(dla1))
	s2 <- exp(-fun.cumsum(dla2))	

	k1all <- matrix(1, nrow = n, ncol = 1)
	k2all <- k1all
	fb1 <- matrix(0, nrow = n, ncol = 1)
	fb2 <- matrix(0, nrow = n, ncol = 1)
	l2 <- fb2
	l1 <- fb1 
	fb1[zd1] <- s1
	k1all[zd1] <- k1
	l1[zd1] <- fun.cumsum(dla1)
	fb2[zd2] <- s2
	k2all[zd2] <- k2
	l2[zd2] <- la2

for  (i in n2:1){
        if (i==n2){
            if (zd2[i]==n){
                k1all[zd2[i]] <- 0
		}
            else{
			k1all[zd2[i]] <- k1all[zd2[i]+1]
		}            
        }
	  else{
		k1all[zd2[i]] <- k1all[zd2[i]+1]
	  }
}

for  (i in n1:1){
        if (i==n1){
            if (zd1[i]==n){
                k2all[zd1[i]] <- 0
		}
            else{
                k2all[zd1[i]] <- k2all[zd1[i]+1]
		}
	  }
        else{
            k2all[zd1[i]] <- k2all[zd1[i]+1]
	  }
}

for  (i in 1:n2) {
        if (i==1){
            if (zd2[i]==1){
                fb1[zd2[i]] <- 1
                l1[zd2[i]] <- 0
		}
            else {
                fb1[zd2[i]] <- fb1[zd2[i]-1]
                l1[zd2[i]] <- l1[zd2[i]-1]
            }
	  }
        else {
            fb1[zd2[i]] <- fb1[zd2[i]-1]
            l1[zd2[i]]=l1[zd2[i]-1]
        }
}

for  (i in 1:n1){
        if (i==1){
            if (zd1[i]==1){
                fb2[zd1[i]] <- 1
                l2[zd1[i]] <- 0
		}
            else {
                fb2[zd1[i]] <- fb2[zd1[i]-1]
                l2[zd1[i]] <- l2[zd1[i]-1]
		}
        }
        else {
            fb2[zd1[i]] <- fb2[zd1[i]-1]
            l2[zd1[i]] <- l2[zd1[i]-1]
        }
}

dff1 <- fun.Rjudge(fb1,0,'>')
dfb1 <- fb1*dff1+(1-dff1)
u0 <- (1 / dfb1-1) 
ta1 <- sum(dff1)

u0[ta1:n] <- matrix(1, nrow = n-ta1+1, ncol = 1)*u0[ta1]
r0 <- max(u0)
lamc <- log(1+u0)

temp1 <- matrix(c(1:n), nrow = n, ncol = 1)
k1n <- temp1[k1all>0]
ka1 <- length(k1n)

temp2 <- matrix(c(1:n), nrow = n, ncol = 1)
k2n <- temp2[k2all>0]
ka2 <- length(k2n)

kk <- n

if (k2n[ka2]<n){
    k2all[(k2n[ka2]+1):n]=k2all[k2n[ka2]]*matrix(1, nrow = n-k2n[ka2], ncol = 1)
}

if (k1n[ka1]<n){
	k1all[(k1n[ka1]+1):n]=k1all[k1n[ka1]]*matrix(1, nrow = n-k1n[ka1], ncol = 1)
}


length <- list(Num1=n1,Num2=n2)
GroupData <- list(length=length)

fb <- list(Num1=fb1,Num2=fb2)

kall <- list(Num1=k1all,Num2=k2all)

l <- list(Num1=l1,Num2=l2)

#u0 <- list(u0=u0)
#kk <- list(kk=kk)
#fb1 <- list(dfb1=dfb1)

output <- list(GroupData=GroupData,fb=fb,kall=kall,u0=u0,kk=kk,l=l,dfb1=dfb1)
return(output)

}
