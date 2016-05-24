#' @export 
#' @method print SO3
print.SO3<-function(x,...){
  Rs<-x
  len<-length(Rs)
  
  #if(len%%9!=0)
  #  stop("Input is not of the correct length.")
  
  if(len==9){
    tr<-matrix(Rs,3,3)
    #class(tr)<-"SO3"
    print.default(tr,...)
  }else{

    tRs<-matrix(Rs,dim(Rs))
    
    if(ncol(tRs)==9){
      cnames<-c("R11","R21","R31","R12","R22","R32","R13","R23","R33")
      colnames(tRs)<-cnames
    }
    print.default(tRs,...)
  }
}

#' @export 
#' @method head SO3
head.SO3<-function(x, n = 6L,...){
  
  #The following two lines are from 'head.matrix'
  stopifnot(length(n) == 1L)
  n <- if(n < 0L)  max(nrow(x) + n, 0L)  else min(n, nrow(x))
  
  x[seq_len(n) ,]
  
}

#' @export 
#' @method tail SO3
tail.SO3<-function(x, n = 6L, addrownums = TRUE,...){
  
  x<-matrix(x,dim(x))
  
  if(ncol(x)==9){
    cnames<-c("R11","R21","R31","R12","R22","R32","R13","R23","R33")
    colnames(x)<-cnames
  } 
  
  tail.matrix(x, n, addrownums, ...)

}

#' @export 
#' @method str SO3
str.SO3<-function(object,...){
  object<-matrix(object,length(object)/9,9)
  str(object)
}

#' @export 
#' @method print Q4
print.Q4<-function(x,...){
  digs<-.Options$digits
  x<-round(x,digs)
  Qs<-x
  len<-length(Qs)
  
  #if(len%%4!=0)
  #  stop("Input is not of the correct length.")
  
  if(len==4){
    Qs<-matrix(Qs,1,4)
    negs<-length(which(Qs[2:4]<0))
    
    if(negs==0){ 
      
      print.default(bquote(.(Qs[1])+.(Qs[2])*i+.(Qs[3])*j+.(Qs[4])*k),...)
      
    }else if(negs==1){
      
      if(Qs[2]<0){
        Qs[2]<--Qs[2]
        print.default(bquote(.(Qs[1])-.(Qs[2])*i+.(Qs[3])*j+.(Qs[4])*k),...)
      }else if(Qs[3]<0){
        Qs[3]<--Qs[3]
        print.default(bquote(.(Qs[1])+.(Qs[2])*i-.(Qs[3])*j+.(Qs[4])*k),...)
      }else{
        Qs[4]<--Qs[4]
        print.default(bquote(.(Qs[1])+.(Qs[2])*i+.(Qs[3])*j-.(Qs[4])*k),...)
      }
      
    }else if(negs==2){
      
      if(all(Qs[2:3]<0)){
        
        Qs[2:3]<-abs(Qs[2:3])
        print.default(bquote(.(Qs[1])-.(Qs[2])*i-.(Qs[3])*j+.(Qs[4])*k),...)
        
      }else if(all(Qs[3:4]<0)){
        
        Qs[3:4]<-abs(Qs[3:4])
        print.default(bquote(.(Qs[1])+.(Qs[2])*i-.(Qs[3])*j-.(Qs[4])*k),...)
        
      }else{
        Qs[2:4]<-abs(Qs[2:4])
        print.default(bquote(.(Qs[1])-.(Qs[2])*i+.(Qs[3])*j-.(Qs[4])*k),...)
      }
    }else{ 
      Qs[2:4]<-abs(Qs[2:4])
      print.default(bquote(.(Qs[1])-.(Qs[2])*i-.(Qs[3])*j-.(Qs[4])*k),...)
    }
    
  }else{
    n<-nrow(Qs)
    p<-ncol(Qs)
    tQs<-matrix(Qs,n,p)
    if(is.null(p)) {
      
      print.default(tQs,...)
      
    }else if(p==4){
      
      colnames(tQs)<-c("Real","i","j","k")
      print.default(tQs,...)
      
    }else{
      
      print.default(tQs,...)
      
    }
  }
}

#print.Q4<-function(Qs,...){

#  len<-length(Qs)

#  if(len%%4!=0)
#    stop("Input is not of the correct length.")

#  if(len==4){
#    print.default(sprintf("%f + %f*i+ %f*j+%f*k ",Qs[1],Qs[2],Qs[3],Qs[4]),...)
#  }else{
#    colnames(Qs)<-c("Real","i","j","k")
#    print.default(Qs,...)
#  }
#}

#' @export 
#' @method head Q4
head.Q4<-function(x,n=6L,...){
  
  #The following two lines are from 'head.matrix'
  stopifnot(length(n) == 1L)
  n <- if (n < 0L)  max(nrow(x) + n, 0L)  else min(n, nrow(x))
  
  x[seq_len(n), ]
}

#' @export 
#' @method tail Q4
tail.Q4<-function(x, n = 6L, addrownums = TRUE,...){
  
  x <- matrix(x, dim(x))
  
  if(ncol(x)==4) colnames(x) <- c("Real","i","j","k")
  
  tail.matrix(x, n, addrownums,...)
}

#' @export 
#' @method str Q4
str.Q4<-function(object,...){
  
  object<-matrix(object,length(object)/4,4)
  str(object)
}

#' @export 
#' @method [ SO3
'[.SO3'<-function(x,i,...){
  x<-matrix(x,dim(x))
  x<-x[i,...]
  class(x)<-"SO3"
  return(x)
}

#' @export 
#' @method [ Q4
'[.Q4'<-function(x,i,...){
  x<-matrix(x,length(x)/4,4)
  x<-x[i,...]
  class(x)<-"Q4"
  return(x)
}

#' @export 
#' @method == Q4
'==.Q4'<-function(e1,e2){

  e1<-matrix(e1,dim(e1))
  e2<-matrix(e2,dim(e2))
  if(all(e1==e2) || all(e1==-e2))  return(TRUE) else return(FALSE)

}

#' Arithmetic operators on SO(3)
#' 
#' These binary operators perform arithmetic on rotations in quaternion or rotation matrix form
#' (or objects which can be coerced into them).
#' 
#' The rotation group SO(3) is a multiplicative group so ``adding" rotations \eqn{R_1}{R1} and \eqn{R_2}{R2}
#' results in \eqn{R_1+R_2=R_2R_1}{R1+R2=R2R1}.  Similarly, the difference between rotations \eqn{R_1}{R1} and \eqn{R_2}{R2} is
#' \eqn{R_1-R_2=R_2^\top R_1}{R1-R2=R2'R1}.  With this definition it is clear that 
#' \eqn{R_1+R_2-R_2=R_2^\top R_2R_1=R_1}{R1+R2-R2=R2'R2R1=R1}.  
#' If only one rotation is provided to subtraction then the inverse (transpose) it returned, 
#' e.g. \eqn{-R_2=R_2^\top}{-R2=R2'}.
#' 
#' @name Arithmetic
#' @aliases "+.SO3" "-.SO3" "+.Q4" "-.Q4"
#' @param x first argument
#' @param y second argument (optional for subtraction)
#' @return  \item{+}{the result of rotating the identity frame through x then y}
#'          \item{-}{the difference of the rotations, or the inverse rotation of only one argument is provided}
#' @examples
#' U <- c(1, 0, 0)          #Rotate about the x-axis
#' R1 <- as.SO3(U, pi/8)    #Rotate pi/8 radians about the x-axis
#' R2 <- R1 + R1            #Rotate pi/8 radians about the x-axis twice
#' mis.axis(R2)             #x-axis: (1,0,0)
#' mis.angle(R2)            #pi/8 + pi/8 = pi/4
#' 
#' R3 <- R1 - R1            #Rotate pi/8 radians about x-axis then back again
#' R3                       #Identity matrix
#' 
#' R4 <- -R1                #Rotate in the opposite direction through pi/8
#' R5 <- as.SO3(U, -pi/8)   #Equivalent to R4
#' 
#' M1 <- matrix(R1, 3, 3)   #If element-wise addition is requred,
#' M2 <- matrix(R2, 3, 3)   #translate them to matrices then treat as usual
#' M3 <- M1 + M2
#' 
#' M1 %*% M1                #Equivalent to R2
#' t(M1) %*% M1             #Equivalent to R3
#' t(M1)                    #Equivalent to R4 and R5
#' 
#' #The same can be done with quaternions: the identity rotation is (1, 0, 0, 0)
#' #and the inverse rotation of Q=(a, b, c, d) is -Q=(a, -b, -c, -d)
#' 
#' Q1 <- as.Q4(R1)
#' Q2 <- Q1 + Q1
#' mis.axis(Q2)
#' mis.angle(Q2)
#' 
#' Q1 - Q1                  #id.Q4 = (1, 0, 0, 0)


NULL

#' @rdname Arithmetic
#' @aliases "-.SO3" "+.Q4" "-.Q4"
#' @export 
#' @method + SO3

'+.SO3'<-function(x,y){
  
  if(length(y)==9){
    
    y<-t(matrix(y,3,3))
    return(center(x,y))
    
  }else if(all(dim(x)==dim(y))){
    y<--y
    return(center(x,y))
    
  }else{
    stop("y can be a sinlge rotation or x and y must be of the same dimension")
  }

}

#' @rdname Arithmetic
#' @aliases "+.SO3" "+.Q4" "-.Q4"
#' @export 
#' @method - SO3

'-.SO3'<-function(x,y=NULL){
  
  #If y is left null, return the transpose of
  #each rotation in x
  if(is.null(y)){
    
      n<-length(x)/9
      
      if(n%%1!=0){
        stop("x must be of dimension n-by-9")
      }
      
      xt<-matrix(x,n,9)
      xt<-xt[,c(1,4,7,2,5,8,3,6,9)]
  
    class(xt)<-"SO3"
    return(xt)  
  }else{
  
    return(center.SO3(x,y))
  
  }
}

#' @rdname Arithmetic
#' @aliases "+.SO3" "-.SO3" "-.Q4"
#' @export 
#' @method + Q4

'+.Q4'<-function(x,y){

  return(as.Q4(as.SO3(x)+as.SO3(y)))
}

#' @rdname Arithmetic
#' @aliases "+.SO3" "-.SO3" "+.Q4"
#' @export 
#' @method - Q4

'-.Q4'<-function(x,y=NULL){
  
  x<-matrix(x,length(x)/4,4)
  if(is.null(y)){ 
    x[,2:4]<- -1*x[,2:4]
    class(x)<-'Q4'
    return(x)
  }

  return(as.Q4(as.SO3(x)-as.SO3(y)))
}
