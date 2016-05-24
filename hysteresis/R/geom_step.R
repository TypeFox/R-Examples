geom_step <- function (pars,x_resid,y_resid,n) {
  #This code has been translated from the matlab work of Gander, Golub, and Strebel.
  S <- sin(pars[1:n])
  C <- cos(pars[1:n])
  rS <- sin(pars[n+1])
  rC <- cos(pars[n+1])
  JA <- cbind(-pars[n+3]*S,C,0,rC,rS)
  JB <- cbind(pars[n+2]*C,0,S,-rS,rC)
  DA <- -pars[n+2]*diag(S)
  DB <- pars[n+3]*diag(C)
  H <- matrix(0,nrow=n+5,ncol=n+5)
  for (i in 1:n) {
    H[i,i] <- cbind(pars[n+2]*C[i],pars[n+3]*S[i])%*%t(cbind(x_resid[i],y_resid[i]))
    H[i,n+1] <- cbind(pars[n+3]*C[i],pars[n+2]*S[i])%*%t(cbind(x_resid[i],y_resid[i]))
  }  
  H[1:n,n+2] <- S*x_resid
  H[1:n,n+3] <- -C*y_resid
  H[n+1,n+1] <- t(c(pars[n+2]*C,pars[n+3]*S))%*%c(x_resid,y_resid)
  H[n+1,n+2] <- t(C)%*%y_resid
  H[n+1,n+3] <- -t(S)%*%x_resid
  H <- H + upper.tri(H)
  DD <- pars[n+2]^2*S^2+pars[n+3]^2*C^2
  J1 <- DA%*%JA + DB%*%JB
  J2 <- rbind(cbind(diag(DD),J1),cbind(t(J1),crossprod(JA)+crossprod(JB)))
  Y2 <- rbind(DA%*%x_resid+DB%*%y_resid,crossprod(JA,x_resid)+crossprod(JB,y_resid))
  change <- solve(J2 + H,Y2)
  pars <- pars + change
  #list(pars,"J"=J2+H)
  list(pars,"J"=J2)
}

