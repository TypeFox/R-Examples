.toleranceSign <- function(dir){
  # Ax<=b+(1-beta)t -> Ax-t(1-beta)<=b
  # Ax>=b-(1-beta)t -> Ax+t(1-beta)>=b
  
  dir2<-(dir=="<=" | dir=="<") # Transform <= or < in 1
  dir3<-(dir==">=" | dir==">") # Transform >= or > in 1
  dir2-dir3 # < and <= become 1, > and >= become -1 and other thing become 0
}

#' @import FuzzyNumbers
.Yager_1 <- function(tfn){
  tfn<-as.TrapezoidalFuzzyNumber(tfn)
  r <- supp(tfn)[1]
  R <- supp(tfn)[2]
  u <- core(tfn)[1]
  U <- core(tfn)[2]

  (r+(u+U)/2+R)/3
}


.Yager_3 <- function(tfn){
  tfn<-as.TrapezoidalFuzzyNumber(tfn)
  r <- supp(tfn)[1]
  R <- supp(tfn)[2]
  u <- core(tfn)[1]
  U <- core(tfn)[2]
  (r+R+u+U)/4
}


.Adamo <- function(tfn, alpha){
  tfn<-as.TrapezoidalFuzzyNumber(tfn)
  R <- supp(tfn)[2]
  U <- core(tfn)[2]
  R-alpha*(R-U)
}


.Average <- function(tfn, lambda, t){
  tfn<-as.TrapezoidalFuzzyNumber(tfn)
  r <- supp(tfn)[1]
  R <- supp(tfn)[2]
  u <- core(tfn)[1]
  U <- core(tfn)[2]
  
  u-(u-r)/(t+1)+lambda*(U-u+((R-r)-(U-u))/(t+1))
}


.Chang <- function(tfn){
  r <- supp(tfn)[1]
  R <- supp(tfn)[2]
  u <- core(tfn)[1]
  U <- core(tfn)[2]
  
  ((U^2-u^2) + (R^2-r^2) + (U*R-u*r)) / 6
}


.Yager_2 <- function(tfn){
  R <- supp(tfn)[2]
  U <- core(tfn)[2]
  
  R/(R-U+1)
}


.Posibi <- function(tfn,w1,w2){
  tfn<-as.TrapezoidalFuzzyNumber(tfn)
  drop(supp(tfn)%*%c(w1,w2))
}


.Interv <- function(tfn,alpha,w1,w2){
  tfn<-as.TrapezoidalFuzzyNumber(tfn)
  r <- supp(tfn)[1]
  R <- supp(tfn)[2]
  u <- core(tfn)[1]
  U <- core(tfn)[2]
  
  zi=u-alpha*(u-r)
  zc=(1/2)*((u+U)-alpha*(u-r)+alpha*(R-U))
  (w1*zi)+(w2*zc)
}


.StratI <- function(tfn, alpha){
  tfn<-as.TrapezoidalFuzzyNumber(tfn)
  r <- supp(tfn)[1]
  u <- core(tfn)[1]
  
  r+alpha*(u-r)
}


.StratS <- function(tfn, alpha){
  tfn<-as.TrapezoidalFuzzyNumber(tfn)
  R <- supp(tfn)[2]
  U <- core(tfn)[2]
  
  R-alpha*(R-U)
}


.MultiobI <- function(tfn, alpha){
  tfn<-as.TrapezoidalFuzzyNumber(tfn)
  r <- supp(tfn)[1]
  u <- core(tfn)[1]
  
  u-alpha*(u-r)
}


.MultiobS <- function(tfn, alpha){
  tfn<-as.TrapezoidalFuzzyNumber(tfn)
  R <- supp(tfn)[2]
  U <- core(tfn)[2]
  
  U+alpha*(R-U)
}

.evalObjective <- function(obj, sol){
  sum=as.TrapezoidalFuzzyNumber(0)
  for (i in 1:length(obj)){
    sum=sum+as.TrapezoidalFuzzyNumber(obj[[i]]*sol[i])
  }
  sum
}
