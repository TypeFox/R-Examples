`Jacob` <- function(X, y, h){
      n = dim(X)[1]
      p = dim(X)[2]
      S = t(X)%*%X
      s = t(X)%*%y
      ai = s
      atildei = matrix(0,p,1)
      Hi = matrix(0,p,p)
      for(i in 1:h){
          if(i == 1){
             dahdy = t(X)
             atildei = s
             datildehdy = dahdy
             qi = t(atildei) %*% S %*% (atildei)
             Hi = (atildei %*% t(atildei)) / as.numeric(qi)
             dbhdy = Hi %*% t(X) + (kronecker(t(s),atildei) %*% datildehdy + 
                     as.numeric(t(s) %*% atildei) * datildehdy) / as.numeric(qi)
                     -2 * as.numeric(t(s) %*% atildei) * Hi %*% S %*% datildehdy / as.numeric(qi)
          }
          if(i > 1){
             imat = diag(1:p)
             dahdy = dahdy - S %*% Hi %*% t(X) - S %*% ((kronecker(t(s),atildei) %*% datildehdy + as.numeric(t(s) %*% atildei)
                     * (imat %*% datildehdy)) - 2 * as.numeric(t(s) %*% atildei) * Hi %*% S %*% datildehdy) / as.numeric(qi)
             ai = ai - S %*% Hi %*% s
             datildehdy = (imat - Hi %*% S) %*% dahdy - ((kronecker(t(ai) %*% S, atildei) %*% datildehdy + 
                          as.numeric(t(ai) %*% S %*% atildei) * imat %*% datildehdy) - 2 * 
                          as.numeric(t(ai) %*% S %*% atildei) * Hi %*% S %*% datildehdy) / as.numeric(qi)
             atildei = (imat - Hi %*% S) %*% ai
             qi = (t(atildei) %*% S %*% atildei)
             Hi = atildei %*% t(atildei) / as.numeric(qi)
             dbhdy = dbhdy + Hi %*% t(X) + ((kronecker(t(s), atildei) %*% datildehdy + as.numeric(t(s) %*% atildei) * 
                     datildehdy) -2 * as.numeric(t(s) %*% atildei) * Hi %*% S %*% datildehdy) / as.numeric(qi)
          }
      }
      structure(list(Jacomatrix = dbhdy))
}
