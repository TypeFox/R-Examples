emfa <-
function(data,nbf,x=1,test=x[1],pvalues=NULL,min.err=0.001) {
if (class(data)[1]!="FAMTdata") stop("Class of data should be FAMTdata")
n = ncol(data$expression)
m = nrow(data$expression)
if (nbf == 0) {
   B=NULL
   Psi=rep(1,m)
   Factors=NULL
   commonvar=0
   SelectH0=NULL
   corfactors=NULL
} 
if (nbf > 0) {
   print(paste("Fitting Factor Analysis Model with",nbf,"factors"))
   rdata = residualsFAMT(data,x=x,test=test,pvalues=pvalues)
   SelectH0 = rdata$SelectH0
   P = rdata$P
   rdata = rdata$residuals
   cdata = scale(rdata)
   eig = svd((1/sqrt((n - 1))) * t(cdata))
   evectors = eig$u[, 1:nbf]
   evalues = eig$d^2
   if (nbf > 1) B = evectors[, 1:nbf] %*% diag(sqrt(evalues[1:nbf]))
   if (nbf == 1) B = matrix(evectors, nrow = m, ncol = 1)*sqrt(evalues[1])
   Psi = 1 - apply(B^2, 1, sum)
   crit = 1
   while (crit > min.err) {
        sB0 = scale(t(B[SelectH0, ]), center = FALSE, scale = sqrt(Psi[SelectH0]))
        G = solve(diag(nbf) + sB0 %*% t(sB0))
        sB0 = scale(t(B[SelectH0, ]), center = FALSE, scale = Psi[SelectH0])
        expf = G %*% sB0 %*% t(cdata[, SelectH0])
        expff = data.frame(lapply(1:n, function(i, M, G) as.vector(G + 
            M[, i] %*% t(M[, i])), M = expf, G = G))
        expff = matrix(apply(expff, 1, sum), nrow = nbf,ncol = nbf)
        Bnew = t(expf %*% cdata) %*% solve(expff)
        expfcdata = expf %*% cdata
        Psinew = 1 - (1/n) * unlist(lapply(1:m, function(i, 
            M1, M2) (M1[i, ] %*% M2[, i])[1, 1], M1 = Bnew, M2 = expfcdata))
        crit = mean((Psi - Psinew)^2)
        B = Bnew
        Psi = Psinew
    }
    if (nbf == 1) {
        hat = solve(expff) %*% expf
    }
    if (nbf > 1) {
        rotation = varimax(B)
        hat = t(rotation$rotmat) %*% solve(expff) %*% expf
        B = rotation$loadings
    }
    sB0 = scale(t(B[SelectH0, ]), center = FALSE, scale = sqrt(Psi[SelectH0]))
    G = solve(diag(nbf) + sB0 %*% t(sB0))
    sB0 = scale(t(B[SelectH0, ]), center = FALSE, scale = Psi[SelectH0])
    Factors = G %*% sB0 %*% t(cdata[, SelectH0])
    S1 = t(Factors) %*% hat %*% P
    P1 = diag(ncol(S1)) - S1
    cz = rep(sum(diag((P1 %*% t(P1)))), m)
    Psi = Psi * (n/cz)
    b2 = apply(B^2,1,sum)
    commonvar = sum(b2)/sum(b2 + Psi)
    Factors = t(Factors)
}
res = list(B=B,Psi=Psi,Factors=Factors,commonvar=commonvar,SelectH0=SelectH0)
return(res)    
}
