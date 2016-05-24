BDM.test <- function (Y, X) 
{
    if(any(is.na(Y))==TRUE) stop("Please remove missing values")
    n <- length(Y)
    n.s <- tapply(Y, X, length)
    Q <- matrix((1/n) * (tapply(rank(Y), X, function(x) {
        mean(x) - (1/2)
    })))
    r <- nlevels(X)
    S.sq <- matrix(ncol = 1, nrow = r)
    ranks <- rank(Y)
    for (i in 1:r) {
        Xn <- as.numeric(X)
        S.sq[i] <- 1/(n^2 * (length(Y[Xn == i]) - 1)) * sum((ranks[Xn == 
            i] - mean(ranks[Xn == i]))^2)
    }
    V <- n * diag(as.vector(S.sq) * (1/n.s))
    M <- diag(1, r, r) - ((1/r) * matrix(1, r, r))
    F.star <- n/sum(diag(M[1] * V)) * (t(Q) %*% M %*% Q)
    nu1 <- M[1]^2 * sum(diag(V))^2/sum(diag(M %*% V %*% M %*% 
        V))
    nu2 <- sum(diag(V))^2/sum(diag(V %*% V %*% diag(1/(n.s - 
        1))))
    res <- list()
    res$head <- "One way Brunner-Dette-Munk test" 	
    res$Q <- data.frame(Levels = levels(X), Rel.effects = Q)
    res$BDM.Table <- data.frame(df1 = nu1, df2 = nu2, F = F.star, 
        p.val = pf(F.star, nu1, nu2, lower.tail = FALSE), row.names = "X")
    colnames(res$BDM.Table) <- c("df1", "df2", "F*", "P(F > F*)")
    if(nrow(res$BDM.Table) == 1) rownames(res$BDM.Table) <- ""	
    class(res) <- "BDM"
    res
}

print.BDM<-function(x,digits= max(3, getOption("digits")), ...){
cat("\n")
cat(x$head,"\n\n")
print(x$BDM.Table,digits=digits)
cat("\n")
invisible(x)
}


BDM.2way<-function(Y,X1,X2){
n<-length(Y)
n.s<-tapply(Y,interaction(X1,X2,lex.order=TRUE),length)
Q<-matrix((1/n)*(tapply(rank(Y),interaction(X1,X2,lex.order=TRUE),function(x){mean(x)-(1/2)})))
ranks<-rank(Y)
J<-nlevels(X1)
K<-nlevels(X2)
L<-J*K
S.sq<-matrix(ncol=1,nrow=L)

for(i in 1:L){
Xn<-as.numeric(interaction(X1,X2))
S.sq[i]<-1/(n^2*(length(Y[Xn==i])-1))*sum((ranks[Xn==i]-mean(ranks[Xn==i]))^2)
}

V<-n*diag(as.vector(S.sq)*(1/n.s))

P.J<-diag(1,J)-(matrix(1,J,J)/J)
P.K<-diag(1,K)-(matrix(1,K,K)/K)
M.A<-kronecker(P.J,matrix(1,K,K)/K)
M.B<-kronecker(matrix(1,J,J)/J,P.K)
M.AB<-kronecker(P.J,P.K)
F.A<-n/sum(diag(M.A[1]*V))*(t(Q)%*%M.A%*%Q)
F.B<-n/sum(diag(M.B[1]*V))*(t(Q)%*%M.B%*%Q)
F.AB<-n/sum(diag(M.AB[1]*V))*(t(Q)%*%M.AB%*%Q)

nu1.A<-M.A[1]^2*sum(diag(V))^2/sum(diag(M.A%*%V%*%M.A%*%V))
nu1.B<-M.B[1]^2*sum(diag(V))^2/sum(diag(M.B%*%V%*%M.B%*%V))
nu1.AB<-M.AB[1]^2*sum(diag(V))^2/sum(diag(M.AB%*%V%*%M.AB%*%V))

nu2<-sum(diag(V))^2/sum(diag(V%*%V%*%diag(1/(n.s-1))))
res<-list()
res$head <- "Two way Brunner-Dette-Munk test" 
res$Q<-data.frame(Levels=levels(interaction(X1,X2,lex.order=TRUE)),Rel.effects=Q)

res$BDM.Table<-data.frame(nu1=c(nu1.A,nu1.B, nu1.AB),nu2=c(nu2,nu2,nu2),F.star=c(F.A,F.B,F.AB),P.val=c(pf(F.A,nu1.A,nu2,lower.tail=FALSE),pf(F.B,nu1.B,nu2,lower.tail=FALSE),
pf(F.AB,nu1.AB,nu2,lower.tail=FALSE)),row.names=c("X1","X2","X1:X2"))
colnames(res$BDM.Table)<-c("df1","df2","F*","P(F > F*)")
class(res) <- "BDM"
res
}