"fun.theo.bi.mv.gld"<-
function(L1,L2,L3,L4,param1,M1,M2,M3,M4,param2,p1,normalise="N"){

if(length(L1)>1){
p1<-L1[9]
M4<-L1[8]
M3<-L1[7]
M2<-L1[6]
M1<-L1[5]
L4<-L1[4]
L3<-L1[3]
L2<-L1[2]
L1<-L1[1]
}

raw.moment1<-fun.rawmoments(L1,L2,L3,L4,param1)
raw.moment2<-fun.rawmoments(M1,M2,M3,M4,param2)
p2<-1-p1

ex1<-p1*raw.moment1[1]+p2*raw.moment2[1]
ex2<-p1*raw.moment1[2]+p2*raw.moment2[2]
ex3<-p1*raw.moment1[3]+p2*raw.moment2[3]
ex4<-p1*raw.moment1[4]+p2*raw.moment2[4]

result <- rep(NA, 4)

result[1]<-ex1
result[2]<-ex2-ex1^2
result[3]<-(ex3-3*ex2*ex1+2*ex1^3)/(result[2]^1.5)
result[4]<-(ex4-4*ex3*ex1+6*ex2*ex1^2-3*ex1^4)/(result[2]^2)

if(param1=="rs"){
check1<-!as.logical((L3 > -1/1:4) * (L4 > -1/1:4))
}
else if(param1=="fmkl"){
check1<-as.logical(min(L3, L4) < (-1/1:4))
}

if(param2=="rs"){
check2<-!as.logical((M3 > -1/1:4) * (M4 > -1/1:4))
}
else if(param2=="fmkl"){
check2<-as.logical(min(M3, M4) < (-1/1:4))
}

result[check1+check2!=0]<-NA

if(normalise=="Y"){
result[4] <- result[4] - 3
}

names(result) <- c("mean", "variance", "skewness", "kurtosis")
return(result)
}

