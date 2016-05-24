##### computes the Q, R and S vectors from Sam's note
#####   input:
#####     - y = Observed response
#####     - A, b, X.E = Object containing all the relevant selection information. Important items: X.E.list (list of reduced predictor matrices in each cluster), A, b (define the selection polyhedron)
compute.QRS.vectors <-
function(y, X.E, groups, test.group, A, b, mu=NULL, type=c("univariate", "multivariate")){
  n = length(y)
  type = type[1] # defaults to univariate test
  
  # construct the projection matrices, depending on the value of mu
  if (is.null(mu)){
    if (type == "univariate"){
      P1 = matrix (1/n, nrow=n, ncol=n)
    }else{ # multivariate
      X1 = cbind (rep(1,n), X.E[,groups!=test.group, drop=FALSE])  
      P1 = X1%*%ginv(X1)
    }
    X.tilde = cbind(rep(1,n), X.E)
    P2 = X.tilde%*%ginv(X.tilde)
  }else{ # multivariate
    if (type == "univariate"){
      P1 = matrix (0, nrow=n, ncol=n)
    }else{
      if (sum(groups != test.group) == 0){ # none of the other group columns selected
        P1 = matrix (0, nrow=n, ncol=n)
      }else{
        X1 = X.E[,groups!=test.group,drop=FALSE]
        P1 = X1%*%ginv(X1)
      }
    }
    P2 = X.E%*%ginv(X.E)
  }
  

  P1.perp = diag(nrow(P1)) - P1
  P2.perp = diag(nrow(P2)) - P2
  df1 = sum(diag(P1.perp)) - sum(diag(P2.perp))
  df2 = sum(diag(P2.perp))
  c = df1/df2
  
  ### all the truncation region paraphernalia
  R1 = as.numeric(P1.perp%*%y)
  R2 = as.numeric(P2.perp%*%y)
  delta = as.numeric(P1%*%y)
  l = sqrt(sum(R1^2))
  
  V.N = as.numeric((P2 - P1)%*%y)
  V.N = V.N/sqrt(sum(V.N^2))
  
  V.D = as.numeric(P2.perp%*%y)
  V.D = V.D/sqrt(sum(V.D^2))
  
  ### construct the necessary vectors
  Q = l*as.numeric(A%*%V.N)
  R = as.numeric(A%*%delta) - b
  S = l*as.numeric(A%*%V.D)
  
  ### F statistic
  F.stat = (sum (R1^2) - sum(R2^2))/sum(R2^2)/c
  
  list(Q=Q, R=R, S=S, R1=R1, R2=R2, c=c, A=A, b=b, F.stat=F.stat, df1=df1, df2=df2)
}
