"hf.inv2" <-
function (data,binsize=1) 
{

N<-binsize

if (is.list(data)){

nsc<-length(data[[2]])
data<-c(data[[1]],data[[2]][nsc])
	}

#^converts $d, $s list into classical dwt vector form, finest first
#if necessary

    a <- 2
    n <- length(data)
    nhalf <- n/2
    J <- logb(n, 2)

res<-NULL
for (i in 1:J){
    res <- c(data[1:nhalf],res)
	data<-data[(nhalf+1):n]
	n<-n/2
	nhalf<-nhalf/2    
}

#factor <- 2^J
#    for (j in 1:(J - 1)) {
#        factor <- c(factor,rep(2^(J - j), times = 2^j))
#    }
#print(factor)

res<-c(data[1],res)
#print(res)

    nhalf <- 1
    n <- 2

    sm <- rep(0, nhalf)
    det <- sm
den<-sm1<-NULL
    for (i in 1:J) {
        sm[1:nhalf] <- res[1:nhalf]
	sm1<-c(sm[1:nhalf],sm1)
        det[1:nhalf] <- res[(nhalf + 1):n]
#	print(factor[nhalf:(2*nhalf-1)])
	den<-sm[1:nhalf]*(N-sm[1:nhalf])/(N)
#	cat("den............",den,"\n")
	  den[den<=0]<-0
        res[2 * (1:nhalf) - 1] <- a/2 * (sm[1:nhalf] + det[1:nhalf]*sqrt(den))
        res[2 * (1:nhalf)] <- a/2 * (sm[1:nhalf] - det[1:nhalf]*sqrt(den))
        
	n <- 2 * n
        nhalf <- 2 * nhalf
    }
return(res)
	
#    return(list(res=res,sm1=sm1))
}

