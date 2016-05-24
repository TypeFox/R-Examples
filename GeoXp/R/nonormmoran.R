`nonormmoran` <- function(y,x,W)
{
# initialisation

x<-as.matrix(x)
y<-as.vector(y)

res <- dim(x)
n <- dim(x)[1]
k <- dim(x)[2]

# la matrice W est-elle normée
ifelse(apply(W,1,sum)==rep(1,n),is.norm<-TRUE,is.norm<-FALSE)
ifelse(is.norm,nb.term<-n,nb.term<-sum(W))
# calcul des éléments nécessaires au calcul du I de Moran

x1 <- qr(x)
#b <- qr.coef(x1,y)
e <- qr.resid(x1,y)


epe <- e %*% e
mi <- (n/nb.term)*(e %*% W %*% e)/epe
#s<-t(rep(1,n))%*%W%*%rep(1,n)

M <- diag(n) - (x %*% solve(t(x) %*% x) %*% t(x))

mat <-  M %*% W 

tmv<-sum(diag(mat))

meani <- (n/nb.term)*tmv/(n-k)
mat1 <- (M %*% W) %*% (M %*% t(W))
mat2 <- (M %*% W) %*% (M %*% W)

tmw1<-sum(diag(mat1))
tmw2<-sum(diag(mat2))

vari <- tmw1 + tmw2 + tmv*tmv
vari <- (n/nb.term)^2*vari/((n-k)*(n-k+2))
vari <- vari - meani*meani
mis <- (mi-meani)/sqrt(vari)
prob <- 1-pnorm(mis)

return(list(nobs=n,nvar=k,morani=round(mi,4),imean=round(meani,4),istat=mis,ivar=round(vari,4),prob=round(prob,4)))

  }

