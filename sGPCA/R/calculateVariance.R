calculateVariance <- function(X,Q,R,U,V,K){
  n = dim(X)[1]
  p = dim(X)[2]
  QXR = Q %*% X %*% R
  denomenator = sum(diag(t(X) %*% QXR))
  var = rep(0,K)


if(K >1){
  for(i in 1:K){
   if(sum(abs(U[,i])) != 0 && sum(abs(V[,i])) != 0){
      Pu = U[,1:i] %*% solve(t(U[,1:i]) %*% Q %*% U[,1:i]) %*% t(U[,1:i])
      Pv = V[,1:i] %*% solve(t(V[,1:i]) %*% R %*% V[,1:i]) %*% t(V[,1:i])
      Xk = Pu %*% QXR %*% Pv
      numerator = sum(diag(t(Xk) %*% Q %*% Xk %*% R))
      var[i] = numerator/denomenator  
     }else if(i > 1){
	var[i]= var[i-1]
	}
   }
}

if(K ==1){
    if(sum(abs(U)) != 0 && sum(abs(V)) != 0){
    Pu = U %*% solve(t(U) %*% Q %*% U) %*% t(U)
    Pv = V %*% solve(t(V) %*% R %*% V) %*% t(V)
    Xk = Pu %*% QXR %*% Pv
    numerator = sum(diag(t(Xk) %*% Q %*% Xk %*% R))
    var[1] = numerator/denomenator
}else if(i >1){
var[i] = var[i-1]
}
}
  return(var)
}
