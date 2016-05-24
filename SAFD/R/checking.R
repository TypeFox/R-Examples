checking <-
function(X,com=1){
 #X...2-dim data.frame with columns "x" and "alpha" 
 #containing the vertexes of the polygonal fuzzy number
 #
 #function checks if X defines a polygonal fuzzy number (steps 1 to 6)
 #and if it is of the desired form - same alpha-levels twice (step 7)
 #function returns 0 in case of error, 1 otherwise
 #
 temp<-rep(1,7)
 #
 #check1: if names are ok
  if(names(X)[1]!=c("x")|names(X)[2]!=c("alpha")){
   temp[1]<-0
      print(paste("input has to be a data frame with column names x and alpha"))
   invisible(c(0))
   }
 #
 #check2: if max and min are ok
 if(temp[1]==1){
  if(min(X$alpha)!=0 | max(X$alpha)!=1) {
  temp[2]<-0
  if(com==1){
   print(paste("input defines no polygonal fuzzy number:"))
    print(paste("max alpha-value has to be 1, min 0"))}
   invisible(c(0))
  }
  }
 #
 #check3: if x is a missing value
  if(min(temp[1:2])==1){
  if(all(complete.cases(X$x))==FALSE) {
    temp[3]<-0
     if(com==1){
    print(paste("input defines missing value"))
     }
    invisible(c(0))
     }
    }
 #
 #check4: if x is increasing
 if(min(temp[1:3])==1){
  dat<-c(X$x[1],X$x)
  shifted<-c(X$x,X$x[nrow(X)])
  if(min(shifted-dat)< -.Machine$double.eps^0.5) {
    temp[4]<-0
    if(com==1){
     print(paste("input defines no polygonal fuzzy number:"))
      print(paste("x values must be non-decreasing"))}
    invisible(c(0))
   }
  }
 #
 #check5: if alpha is increasing from 0 to 1
 if(min(temp[1:4])==1){
  A<-subset(X,X$alpha==1)
  cut1<-min(as.numeric(row.names(A)))
  cut2<-max(as.numeric(row.names(A)))
  fromto<-seq(cut1,cut2,by=1)
  difference<-setdiff(fromto,as.numeric(row.names(A)))
  Left<-X[1:cut1,]
  Right<-X[cut2:nrow(X),]
  dat<-c(Left$alpha[1],Left$alpha)
  shifted<-c(Left$alpha,Left$alpha[nrow(Left)])
  if(min(shifted-dat)<0|length(difference)>0) {
   temp[5]<-0
   if(com==1){
    print(paste("input defines no polygonal fuzzy number:"))
     print(paste("alpha-levels must increase from 0 to 1 and decrease from 1 to 0"))}
    invisible(c(0))
   }
  } 
 #check6: if alpha is decreasing from 1 to 0
 if(min(temp[1:5])==1){ 
  dat<-c(Right$alpha[1],Right$alpha)
  shifted<-c(Right$alpha,Right$alpha[nrow(Right)])
  if(min(dat-shifted)<0) {
   temp[6]<-0
   if(com==1){
   print(paste("input defines no polygonal fuzzy number:"))
     print(paste("alpha-levels must increase from 0 to 1 and decrease from 1 to 0"))}
   invisible(c(0))
  }
 } 
 if(min(temp[1:6])==1){ 
   nl<-nrow(X)/2
   if(nl!=trunc(nl,0)|max(abs(X$alpha[1:nl]-X$alpha[(2*nl):(nl+1)]))>0){
    temp[7]<-0
    if(com==1){
      print(paste("input fuzzy number is not of the desired form:"))
      print(paste("even number nr of rows and X$alpha[1:(nr/2)]=X$alpha[nr:(nr/2+1)]."))
      print(paste("Use translator function to convert in the correct form"))
    }
    invisible(c(0))
   }
 }
 invisible(min(temp))
}
