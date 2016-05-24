`sturges.freq` <-
function (x,k=0) 
{
    n <- length(x)
    if (k==0) k <- round(1+log(n,2),0)
    p<- floor(log(abs(median(x,na.rm=TRUE)),10))
	x<-x/10^(p-1)
	maximo <- max(x,na.rm=TRUE)
	minimo <- min(x,na.rm=TRUE)
	min1<-floor(minimo)
	max1<-ceiling(maximo)
	amplitud <- max1 - min1
	tic <- round(amplitud/k,1)
	clases <- seq(min1, max1, tic)
    if (maximo > clases[length(clases)]) {
	clases <- c(clases, clases[length(clases)] + tic)
    }
	k <- length(clases)-1
	maximo<-maximo*10^(p-1);minimo<-minimo*10^(p-1);tic=tic*10^(p-1)
	clases<-clases*10^(p-1); amplitude=amplitud*10^(p-1)
    lista <- list(maximum = maximo, minimum = minimo, amplitude = amplitud, 
            classes = k, interval = tic, breaks = clases)
    return(lista)
}

