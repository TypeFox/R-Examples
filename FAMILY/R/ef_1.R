fm2v<- function(mat){
  diag(mat)[-1]<- NaN
  as.vector(na.omit(as.vector(mat)))
}

v2fm<- function(vec,p){

  fc<- vec[1:(p+1)]
  ind<- (1:p)*p+2
  fr<- c(fc[1], vec[ind])
  
  vec2<- vec[-ind]
  vec2<- vec2[-(1:(p+1))]
  
  mat<- matrix(0, ncol = p, nrow = p)
  diag(mat)<- NaN
  
  for(i in 1:p){
    if(i==1){
      mat[-1,i]<- vec2[1:(p-1)]
    }else if(i==p){
      mat[-p,i] <- tail(vec2,p-1)
    }else{
      mat[1:(i-1), i  ] <- vec2[(i-1)*p - (i-2) + 0:(i-2) ]
      mat[(i+1):p, i ]<- vec2[(i-1)*p +1 + 0:(p-i-1) ]
    } 
    
  }
  
  fin<- matrix(0, ncol = p+1, nrow = p+1)
  fin[1,]<- fr
  fin[,1]<- fc
  fin[-1,-1]<- mat
  return(fin)
}



soft.thres<- function(vec,scal){
  #print(vec)
  sign(vec)*pmax((abs(vec) - scal),0)
}

max_l<- function(l,grp,y_vec,alpha,p){
  n<- length(y_vec)
  vec<- as.vector(t(grp)%*%y_vec)
  for(i in 1:length(l)){
    l[i] <- sum((soft.thres(vec/n,l[i]*alpha))^2) - p*(1-alpha)^2*l[i]^2  
  }
  return(l)
}
