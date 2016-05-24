"ifa.init.pca" <-
function(xx,L,scaling=TRUE)
{
if (scaling) xx<-scale(xx)
S<-var(xx)  
R<-cor(xx)
inv.R<-solve(R)
hii <- 1 - (1/diag(inv.R))
hii<-hii*diag(S)
diag(S)<-hii
ALA<-S
if (L>1) {
eig.L<-(diag(eigen(ALA)$values[1:L])^(0.5))
eig.L<-ifelse(is.na(eig.L),runif(1,0,1),eig.L)
H<-(eigen(ALA)$vectors[,1:L])%*%eig.L
psi <- var(xx) - (H %*% t(H))
}

if (L==1) {
eig.L<-((eigen(ALA)$values[1:L])^(0.5))
eig.L<-ifelse(is.na(eig.L),runif(1,0,1),eig.L)
H<-as.matrix(eigen(ALA)$vectors[,1:L])*eig.L
psi <- var(xx) - (H %*% t(H))
}

output<-list(psi=psi,H=H)
}

