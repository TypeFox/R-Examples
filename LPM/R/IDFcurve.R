IDFcurve <- function(rain ,g, s, t, stvalue1=1, stvalue2=fre, fre, Tr=200, MP)
{rain <- as.matrix(rain)
Ev1 <- eventi(rain,g,s)
Ev2 <- eventi2(Ev1)
nr <- nrow(Ev2)
vb1 <- vb(nr-1)
ris <- nlm(f3,c(vb1,stvalue1),x3=Ev2,fre=fre,iterlim=1000)
if (MP) ris=nlm(f2,c(ris$estimate,stvalue2),x3=Ev2,fre=fre,iterlim=1000)
mu <- mean(ris$estimate[1:nr])-0.45006*sd(ris$estimate[1:nr])
sigma <- sd(ris$estimate[1:nr])/1.2825
a <- qGumbel(1-1/Tr,sigma=sigma, mu= mu)

if (MP)
{b <- ris$estimate[nr+2]
m <- ris$estimate[nr+1]
min <- ris$minimum}
else
{b <- 0
m <- ris$estimate[nr+1]
min <- ris$minimum}
aH <- Hreg(Ev2,ris$estimate[1:nr],m,b,fre)
h <- a*t/(b+t)^m
i <- a/(b+t)^m
out <- list(par=c(a,m,h,i,min),Curve=aH)
options(warn = -1)
ts.plot(ts(t(Ev2[1:10,]),start=fre,end=fre*g,frequency=1/fre),type="p")
for (w in 1:10)
lines(x=seq(fre,fre*g,fre),y=aH[w,],col="red")
cat("a(Tr) = ", a, "\n")
if (MP)
cat("m = ", m, "\n")
else
cat("n = ", 1-m, "\n")
cat("b = ", b, "\n")
cat("h(t) = ", h, "\n")
cat("i(t) = ", i, "\n")
cat("Offset ", min, "\n")
out
}


eventi <- function (x, g, s) 
{
    n <- nrow(as.matrix(x))
    out <- matrix(0, n, g)
    k = 1
    while (k < (n - g)) {
        if (x[k] > s) {
            out[k, 1] <- x[k]
            for (j in 1:(g - 1)) out[k, j + 1] <- out[k, j] + 
                x[k + j]
            k = k + g
        }
        else {
            out[k, 1] <- 0
            k = k + 1
        }
    }
    out
}


eventi2 <- function (x) 
{
    cont = 0
    x=as.matrix(x)
    n = nrow(x)
    out2 = matrix(0, n, ncol(x))
    for (i in 1:n) {
        if (x[i, 1] > 0) {
            cont = cont + 1
            out2[cont, ] = x[i, ]
        }
    }
    out3 <- out2[1:cont, ]
    out4 <- out3
    for (j in 1:ncol(out3))
    out4[,j]=sort(out3[,j],decreasing = T )
    out4
}

f1 <- function (x1, x2) 
{
    x1 <- as.matrix(x1)
    x2 <- as.matrix(x2)
    v <- 0
    r <- nrow(x1)
    c <- ncol(x1)
    {
        for (j in 1:c) {
            for (i in 1:r) v <- v + (x1[i, j] - x2[i, j])^2
        }
    }
    v
}

f2 <- function (e, x3,fre) 
{
    e <- as.matrix(e)
    r <- nrow(e)
    a <- e[1:(r - 2)]
    m <- e[r-1]
	b <- e[r]

    x2 <- Hreg(x3, a, m, b,fre)
    out <- f1(x3, x2)
    out
}

f3 <- function (e, x3,fre) 
{
    e = as.matrix(e)
    r = nrow(e)
    a = e[1:(r - 1)]
    b = 0
	m = e[r]

    x2 <- Hreg(x3, a, m, b, fre)
    out = f1(x3, x2)
    out
}

Hreg <- function (ev, a, m, b,fre) 
{
    out = matrix(0, nrow(ev), ncol(ev))
    for (i in 1:ncol(ev)) out[, i] = a * (fre * i)/(b + fre * i) ^m
    out
}

vb <- function(s)
{out=c(1)
for (j in 1:s)
out = c(out,1)  
out
}

