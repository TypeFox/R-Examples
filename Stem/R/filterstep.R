`filterstep` <-
function(z,Fmat,Gmat,Vt,Wt,mx,Cx,XXXcov,betacov,flag)
  {
    a <- Gmat %*% t(mx)
    R <- Gmat %*% Cx %*% t(Gmat) + Wt
    f <- t(Fmat) %*% a
    Q <- t(Fmat) %*% R %*% Fmat + Vt

    if(flag==FALSE)        e <- z - f
    else                 e <- z - (XXXcov %*% betacov) - f

    A <- R %*% Fmat %*% solve(Q)
    m <- a + A%*%e
    C <- R - A%*%Q%*%t(A)

    if (length(z)>1  &&  flag==TRUE)  loglikterm <-log(dmvnorm(as.numeric(z),as.numeric(f+(XXXcov %*% betacov)),Q))
    if (length(z)>1  &&  flag==FALSE)  loglikterm <-log(dmvnorm(as.numeric(z),as.numeric(f),Q))
    if (length(z)==1 &&  flag==TRUE)  loglikterm <- -0.5*( log(2*pi) + log(Q) + (z - (f-XXXcov %*% betacov)) ^2/Q)
    if (length(z)==1 &&  flag==FALSE)  loglikterm <- -0.5*( log(2*pi) + log(Q) + (z-f)^2/Q)

  list(m=m,C=C,loglikterm=loglikterm)
  }

