new.lambda <- function(penden.env) {
  lambda <- get("lambda",penden.env)
  if(get("lambda",penden.env)!=0) {
    eps <- 0.04*get("lambda",penden.env)
  }
  else eps <- 1
  epsdf <- 0
  p <- get("p",penden.env)
  calc <- TRUE
  u <- t(get("ck.val.temp",penden.env))%*%get("DDD.sum",penden.env)%*%get("ck.val.temp",penden.env)
  pen.mat <- get("DDD.sum",penden.env)
  help2 <- get("eigen.pen.mat",penden.env)
  index <- get("index.eigen.pen.mat",penden.env)
  #Utilde <- Matrix(get("Utilde.eigen.pen.mat",penden.env))
  Utilde <- get("Utilde.eigen.pen.mat",penden.env)
  t.Utilde <- t(Utilde)
  diag.help2 <- diag(help2$values[index])
  hh <-1
  while(calc) {
    if(hh==31) {
      assign("df.val",df.val,penden.env)
      break
    }
    Derv2(penden.env,temp=TRUE,lambda=lambda[hh])
    #browser()
    df.val <- sum(diag(x=solve(a=t.Utilde%*%(-get("Derv2.cal.temp",penden.env))%*%Utilde+lambda[hh]*diag.help2,tol=1e-50)%*%(t.Utilde%*%(-get("Derv2.cal.temp",penden.env))%*%Utilde)))
    if(df.val < epsdf) {
      print("df kleiner 0")
      assign("df.val",df.val,penden.env)
      return(lambda[hh-1])
    }
    help <- abs(df.val/u - lambda[hh])
    help2 <- abs(df.val/u - lambda[1])
    if((df.val/u)<0) {
      assign("df.val",df.val,penden.env)
      return(lambda[hh-1])
    }
    if(help<eps) {
      calc <- FALSE
      assign("df.val",df.val,penden.env)
    }
    if((help2<eps)&hh>5) {
      assign("df.val",df.val,penden.env)
      assign("lambda.out",TRUE,penden.env)
      lambda[hh+1]<-df.val/u
      #print("out")
      #print(lambda)
      return(lambda[1])
    }
    else {
      lambda[hh+1] <- df.val/u
      hh <- hh+1
    } 
  }
  return(lambda[hh])
}
