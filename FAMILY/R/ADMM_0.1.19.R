
#This is the main function for updating B in the ADMM algorithm
#We evaluate this using the SVD of \tilde{W} as defined in Appendix I of paper.
update_B.svd<- function(y = numeric(),svd.w,B = matrix(),
                        matD=matrix(),matE=matrix(),matF=matrix(),
                        g.list = list(),rho = 0, quad = TRUE){
  
  #NOTE: The svd.w object is modified in the parent function to be a list which
  #contains the transpose of U and V instead of evaluating it each time.
  p<- ncol(B)-1
  num<- length(y)
  
  #The update of B in Eqn (30) is a sum of two expressions.
  
  #We have the first part in the equation
  first.part<- (svd.w$d)*(svd.w$tu%*%y)
  first.part<- first.part - (svd.w$d^2/(3*rho*num+svd.w$d^2))*first.part
  first.part<- (svd.w$v%*%first.part)/(3*rho*num)
  
  #We now have the second part:
  if(quad == FALSE){
    big.exp<- fm2v(rho*(matD+matE+matF) - g.list[[1]] -g.list[[2]]-g.list[[3]])
  }else{
    big.exp <- as.vector(rho*(matD+matE+matF) - g.list[[1]] -g.list[[2]]-g.list[[3]])
  }
  
  second.part<- svd.w$tv%*%big.exp
  second.part<- (svd.w$d^2/(3*rho*num+svd.w$d^2))*second.part
  second.part<- svd.w$v%*%second.part
  second.part<- (big.exp-second.part)/(3*rho)
  
  #Return the update of B in matrix form
  b.new<- as.vector(first.part+second.part)
  
  if(quad==FALSE){
    B.new<- v2fm(b.new,p = p)
  }else{
    B.new<- matrix(b.new, ncol = ncol(B), nrow = nrow(B) )
  } 
  
  return(B.new)
}



#####################################################################################

#The function I use for the l_infinity prox function.

l_inf_prox<- function(y,lam,na.rm = TRUE){
  had.na<- FALSE
  if(na.rm == TRUE){
    ind.na <- which(is.na(y))
    if(length(ind.na) != 0){
      had.na<- TRUE
      y<- y[-ind.na]
    } 
  }
  
  sorted.y<- sort(abs(y),TRUE)
  ans<- (cumsum(sorted.y)-lam)/(1:length(y))
  
  if(all(ans<=0)){
    ret.obj<- rep(0,length(y)) 
    
    if(had.na == TRUE){
      ret.obj<- append(ret.obj,NaN, after = ind.na-1 )
    }
    return(ret.obj)
  }
  statement<- ans > sorted.y
  if(all(statement==FALSE)){
    c<- tail(ans,1)
    ret.obj<- sign(y)*(pmin(c,abs(y)))
    
    if(had.na == TRUE){
      ret.obj<- append(ret.obj,NaN, after = ind.na-1 )
    }
    return(ret.obj)
    
  }else{
    ind<- which(statement==TRUE)[1]   
    c<- ans[ind-1]
    ret.obj<- sign(y)*(pmin(c,abs(y)))
    
    if(had.na == TRUE){
      ret.obj<- append(ret.obj,NaN, after = ind.na-1 )
    }
    return(ret.obj) 
  }
}

#####################################################################################


#The update of D for ADMM.
update_D<- function(y = numeric(),W = matrix(),B = matrix(),
                    g.list = list(),rho = 0,lambda, norm = "l2"){
  
  if(norm == "l2"){
    D.new<- B
    #First we opimize the first row (row 0 in paper) as done in Step 3(c) of Section 5.1
    D.new[1,]<- B[1,]+g.list[[1]][1,]/rho
    
    #Now for the rest of matrix via soft shrinkage
    new.mat<- B[-1,] + g.list[[1]][-1,]/(rho)
    row.norm<- sqrt(apply(new.mat^2,1,sum,na.rm = T))
    coef.term <- pmax(1-(lambda[1]/rho)/(row.norm) , 0)
    
    new.mat2<- scale(t(new.mat),center = FALSE,scale = 1/coef.term)
    D.new[-1,]<- t(new.mat2)
    return(D.new)
  }else if(norm == "l_inf"){
    D.new<- B
    #First we opimize the first row (row 0 in paper) as done in Step 3(c) of Section 5.1
    D.new[1,]<- B[1,]+g.list[[1]][1,]/rho
    
    #Now for the rest of matrix via soft shrinkage
    new.mat<- B[-1,] + g.list[[1]][-1,]/(rho)
    D.new[-1,]<- t(apply(new.mat,1,l_inf_prox,lam = lambda[1]/rho))
    return(D.new)
  }else{
    stop("Unrecognized norm.")
  }
  
}

#####################################################################################

#The update of E for ADMM
update_E<- function(y = numeric(),W = matrix(),B = matrix(),
                    g.list = list(),rho = 0,lambda,norm = "l2"){
  
  if(norm == "l2"){
    #Again, we take care of the first column
    E.new<- B
    E.new[,1]<- B[,1]+g.list[[2]][,1]/rho
    
    #Now we have the main part.
    new.mat<- B[,-1]+g.list[[2]][,-1]/rho
    col.norm<- sqrt(apply(new.mat^2,2,sum,na.rm = T))
    
    coef.term<- pmax(1-(lambda[2]/rho)/(col.norm),0)
    new.mat2<- scale(new.mat, center = FALSE, scale = 1/coef.term)
    E.new[,-1]<- new.mat2
    return(E.new)
  }else if(norm == "l_inf"){
    #Again, we take care of the first column
    E.new<- B
    E.new[,1]<- B[,1]+g.list[[2]][,1]/rho
    
    #Now we have the main part.
    new.mat<- B[,-1]+g.list[[2]][,-1]/rho
    E.new[,-1]<- apply(new.mat,2,l_inf_prox,lam = lambda[2]/rho)
    return(E.new)
    
  }else{
    stop("Unrecognized norm.")
  }
  
}

#####################################################################################

#Update of F in ADMM.
update_F<- function(y = numeric(),W = matrix(),B = matrix(),
                    g.list = list(),rho = 0,lambda){
  
  #this assignment takes care of the first row and first column. 
  F.new<- B  +  g.list[[3]]/(rho)
  
  #Now we take the part F_{-0,-0}
  temp.F<- F.new[-1,-1]
  temp.F2<-  sign(temp.F)*pmax(abs(temp.F)-lambda[3]/rho,0)
  F.new[-1,-1]<- temp.F2
  return(F.new)
  
}

#####################################################################################

#Update the dual variable here
update_g.list<- function(y = numeric(),W = matrix(),B = matrix(),
                    matD=matrix(),matE=matrix(),matF=matrix(),
                    g.list = list(),rho = 0){
  #Maybe removing this assignment might help 
  g.new<- g.list
  g.new[[1]]<- g.list[[1]]+rho*(B-matD)
  g.new[[2]]<- g.list[[2]]+rho*(B-matE)
  g.new[[3]]<- g.list[[3]]+rho*(B-matF)
  return(g.new)
  
}



#####################################################################################

#A small function I use to generate W 
quick.func<- function(xz = c(),xn){
  as.vector(xz[1:xn]%o%xz[-(1:xn)])
}

#Generate the matrix \widetilde{W} as used in Appendix I for use in the function.
generate.w<- function(X=matrix(),Z=matrix(), quad = TRUE){
  
  p1<- ncol(X)
  p2<- ncol(Z)
  
  if(quad == FALSE){
    p<- p1
    if(p1!=p2) stop("To remove quadtratic terms p1 must be equal to p2")
    ind<- (1:p)*p + (2*(1:p)+1)
  }
  
  #Just in case we have only one oberservation? Not sure why I did this
  if(is.vector(X)) p1<- length(X)
  if(is.vector(Z)) p2<- length(Z)
  
  #Add the intercept 
  x<- cbind(1,X)
  z<- cbind(1,Z)
  
  W<- t(apply(cbind(x,z),1,quick.func,xn= p1+1))
  if(quad == FALSE){
    W<- W[,-ind]
  }
  return(W)
  
}



########################################################################################

#Main function I use within this package though this function should not be visible to the user
estimate.SH<- function(x,z,y=numeric(),w,svd.w,
                      lambda = c(1,1,1), rho = 1, 
                      B,matD,matE,matF,g.list,iter = 100,
                      e.abs = 1e-4, e.rel=1e-4,quad = TRUE,norm = "l2"){
  
  XequalZ = FALSE
  if(all(dim(x)==dim(z))){
    if(sum(abs(x-z)) < 1e-5){
      XequalZ = TRUE
    }
  }
  
  num<- length(y)  
  for(i in 1:iter){
    #print(i)
    #Update each variable
    newB<- update_B.svd(y,svd.w,B,matD,matE,matF,g.list,rho,quad = quad)
    newmatD<- update_D(y ,w,newB ,g.list,rho ,lambda,norm = norm)
    newmatE<- update_E(y ,w,newB ,g.list,rho ,lambda,norm = norm) 
    newmatF<- update_F(y ,w,newB ,g.list,rho ,lambda)
    newg.list<- update_g.list(y ,w,newB ,newmatD,newmatE,newmatF,g.list,rho)
  
    #dual residual s. As defined in Section 4.3 of paper 
    d.diff<- sum((-rho*(newmatD-matD))^2,na.rm = TRUE)
    e.diff<- sum((-rho*(newmatE-matE))^2,na.rm = TRUE)
    f.diff<- sum((-rho*(newmatF-matF))^2,na.rm = TRUE)
    s.norm<- sqrt(d.diff+e.diff+f.diff) 

    
    
    
    #primal residual r. As defined in the paper
    d.dif2<- sum((newB-newmatD)^2,na.rm = TRUE)
    e.dif2<- sum((newB-newmatE)^2,na.rm = TRUE)
    f.dif2<- sum((newB-newmatF)^2,na.rm = TRUE)
    r.norm<- sqrt(d.dif2+e.dif2+f.dif2) 
  
    other<- c(newmatD, newmatE, newmatF)
    crit1<- max(sqrt(sum(newB^2)), sqrt(sum(other^2)))*e.rel + sqrt(length(B))*e.abs
    crit2<- sqrt(sum(unlist(newg.list)^2))*e.rel + sqrt(length(B))*e.abs
    
    #print(c(s.norm, r.norm))
    #Stopping criteria as defined in Boyd paper using primal and dual residual.
    if(r.norm< crit1 & s.norm < crit2 ){
      finB<- newB*(newmatD!=0)*(newmatF!=0)*(newmatE!=0)
      if(quad == FALSE){
        diag(finB)[-1]<- 0
      }
      returnlist<- list("finB" = finB,"B"= newB, "D" = newmatD, "E" = newmatE, 
                        "F" = newmatF,"glist" = g.list,
                        "rho" = rho,"conv" = TRUE,iters = i,"equal" = XequalZ)
      class(returnlist)<- "interac"
      return(returnlist)
    }else{
      #If not converged then move on to the next update
      B<- newB
      matD<- newmatD
      matE<- newmatE
      matF<- newmatF
      g.list<- newg.list
      
      #Update rho by criteria descriped in paper APpendix. 
      if( r.norm > 10*s.norm ){
        rho<- 2*rho
      }else if(r.norm*10 < s.norm ){
        rho<- rho/2
      }
      
    }
  }
  
  finB<- newB*(newmatD!=0)*(newmatF!=0)*(newmatE!=0)
  if(quad == FALSE){
    diag(finB)[-1]<- 0
  }
  returnlist<- list("finB" = finB,"B"= B, "D" = matD, "E" = matE, "F" = matF, "glist" = g.list,
                    "rho" = rho, "conv" = FALSE,iters=i,"equal" = XequalZ)
  
  class(returnlist)<- "interac"
  return(returnlist)
  
}

#####################################################################################

# THe main function will we a list of interac objects: the output of the previosu functions
#Thses objects will display the following:
print.interac<- function(x, ...){
  if(x$conv){
    cat(paste("ADMM algorithm converged in ",x$iters, " iterations.\n", sep = ""))
    cat("The fitted B matrix is given by:\n")
    print(x$finB)
    cat("\n To view all objects in the list use names() and use plot() to for a graphic of B\n")
  }else{
    paste("Algorithm DID NOT converge in ",x$iters, " iterations.",sep = "")
    cat("The fitted B matrix is given by:\n")
    print(x$finB)
    cat("\n To view all objects in the list use names() and use plot() to for a graphic of B\n")
  }
  invisible(x)
}
#####################################################################################
draw.heatmap<- function(mat,equal = FALSE){
  
  rownames(mat)<- c("inter" , paste("X",1:(nrow(mat)-1),sep = ""))
  if(equal){
    colnames(mat)<- rownames(mat)
  }else{
    colnames(mat)<- c("inter" , paste("Z",1:(ncol(mat)-1),sep = ""))
  }
  
  
  hm2 <- pheatmap(as.matrix(mat), scale="none", 
                  cluster_rows=F, cluster_cols=F,
                  asp = 1)
}


plot.interac<- function(x, ...){
  
  draw.heatmap(x$finB,equal = x$equal)
  
  
}

#####################################################################################

#Now moving on to logistic regression
#The update of B for logistic regression
update_B.logistic<- function(y = numeric(),w = matrix(),svd.w = list(),B = matrix(),
                             matD=matrix(),matE=matrix(),matF=matrix(),
                             g.list = list(),rho = 0,MAX = 100,B.tol, quad = TRUE){
  p1<- nrow(B)
  p2<- ncol(B)
  
  num<- length(y)
  l.dot<- function(y,Xb){
    1/(1+exp(-Xb))-y 
  }
  
  #Initialize B
  
  if(quad== FALSE){
    current.b<- fm2v(B)
  }else{
    current.b<- as.vector(B)
  }
  
  #Formulate new response vector which we need to solve for
  #We break up response vector into two parts for computational feasibility
  W.current.b<- w%*%current.b
  response.vec1<- (1/(2*sqrt(num)))*(W.current.b - 4*l.dot(y,W.current.b))
  
  if(quad == FALSE){
    response.vec2<- fm2v( (rho*(matD+matE+matF)- g.list[[1]] -g.list[[2]]-g.list[[3]])/(sqrt(3*rho)) )
  }else{
    response.vec2<- as.vector( (rho*(matD+matE+matF)- g.list[[1]] -g.list[[2]]-g.list[[3]])/(sqrt(3*rho)) )
  }
  
  #This is the same as the first and second part for linear regression
  #In this case the second part does not change for every iteration of newton raphson
  second.part<- svd.w$tv%*%response.vec2
  second.part<- (svd.w$d^2/(12*num*rho+svd.w$d^2))*second.part
  second.part<- svd.w$v%*%second.part
  second.part<- (response.vec2 - second.part)/sqrt(3*rho)
  
  #Main loop for the iterative algorithm
  for(i in 1:MAX){
    #Again, we have the update for first part only now this is iterative for B
    first.part<- (svd.w$d/(2*sqrt(num)))*(svd.w$tu%*%response.vec1)
    first.part<- first.part - (svd.w$d^2/(12*num*rho+svd.w$d^2))*first.part
    first.part<- (svd.w$v%*%first.part)/(3*rho)
    
      
    #Define updated B
    b.new<- as.vector(first.part + second.part)
    
    #Check for convergence
    if(all(abs(current.b-b.new)< B.tol)){
      if(quad==FALSE){
        return.mat<- v2fm(b.new, p = p1-1)
      }else{
        return.mat<- matrix(b.new, ncol = p2,nrow = p1)
      }
      
      return(list("ans" = return.mat,"iter" = i))
    }else{
      current.b<- b.new
      W.new.b<- w%*%b.new
      response.vec1<- (1/(2*sqrt(num)))*(W.new.b - 4*l.dot(y,W.new.b))
    }  
  }
  if(quad==FALSE){
    return.mat<- v2fm(current.b, p = p1-1)
  }else{
    return.mat<- matrix(current.b, ncol = p2,nrow = p1)
  }
  
  return(list("ans" = return.mat,"iter" = MAX))
  
}


#####################################################################################

#estimate model for logistic regression
#Note: that this could easily be combined with the previous function to get 
#one main function. 
estimate.logistic<- function(x,z,y=numeric(),w,svd.w,
                             lambda = c(1,1,1), rho = 1, 
                             B,matD,matE,matF,g.list,iter = 100,e.abs = 1e-4,
                                e.rel= 1e-4,maxiter.B = 10,tol.B = 1e-4,quad = TRUE,
                             norm = "l2"){
  
  XequalZ = FALSE
  if(all(dim(x)==dim(z))){
    if(sum(abs(x-z)) < 1e-5){
      XequalZ = TRUE
    }
  }
  
  for(i in 1:iter){
    
    #Using the logistoc regression update for B
    newB<- update_B.logistic(y,w,svd.w,B ,matD,matE,matF,g.list,rho,MAX = maxiter.B,B.tol = tol.B,
                             quad = quad)$ans 
    #print(dim(newB))
    newmatD<- update_D(y ,w,newB ,g.list,rho ,lambda,norm = norm)  
    newmatE<- update_E(y ,w,newB ,g.list,rho ,lambda,norm = norm) 
    newmatF<- update_F(y ,w,newB ,g.list,rho ,lambda)
    newg.list<- update_g.list(y ,w,newB ,newmatD,newmatE,newmatF,g.list,rho)
    
    #dual residual s. As defined in the paper Appendix 
    d.diff<- sum((-rho*(newmatD-matD))^2,na.rm = TRUE)
    e.diff<- sum((-rho*(newmatE-matE))^2,na.rm = TRUE)
    f.diff<- sum((-rho*(newmatF-matF))^2,na.rm = TRUE)
    s.norm<- sqrt(d.diff+e.diff+f.diff) 
    
    
    #primal residual r. As defined in the paper
    d.dif2<- sum((newB-newmatD)^2,na.rm = TRUE)
    e.dif2<- sum((newB-newmatE)^2,na.rm = TRUE)
    f.dif2<- sum((newB-newmatF)^2,na.rm = TRUE)
    r.norm<- sqrt(d.dif2+e.dif2+f.dif2) 
    
    
    #Stopping creiteria
    if(r.norm< e.abs & s.norm < e.rel ){
      finB<- newB*(newmatD!=0)*(newmatF!=0)*(newmatE!=0)
      if(quad == FALSE){
        diag(finB)[-1]<- 0
      }
      
      returnlist<- list("finB" = finB,"B"= newB, "D" = newmatD, "E" = newmatE, 
                        "F" = newmatF,"glist" = g.list,
                        "rho" = rho,"conv" = TRUE,iters = i,"equal" = XequalZ)
      class(returnlist)<- "interac"
      return(returnlist)
    }else{
      #If not converged then move on
      B<- newB
      matD<- newmatD
      matE<- newmatE
      matF<- newmatF
      g.list<- newg.list
      
      #Update rho by criteria descriped in paper APpendix. 
      if( r.norm > 10*s.norm ){
        rho<- 2*rho
      }else if(r.norm*10 < s.norm ){
        rho<- rho/2
      }
      
    }
    
  }
  
  finB<- newB*(newmatD!=0)*(newmatF!=0)*(newmatE!=0)
  if(quad == FALSE){
    diag(finB)[-1]<- 0
  }
  
  returnlist<- list("finB" = finB,"B"= B, "D" = matD, "E" = matE, "F" = matF, "glist" = g.list,
                    "rho" = rho, "conv" = FALSE,iters=i,"equal" = XequalZ)
  class(returnlist)<- "interac"
  return(returnlist)
  
}




########################################################################################

#The MAIN function. This will solve for a path of variables.
#This is main function a use will use by supplying a vecto of values for alpha and a vector
#of values for lambda along which one wishes to solve. 
FAMILY<- function(X,Z,Y =numeric(),lambdas = c(),alphas = c(),
                      family = c("gaussian","binomial"),
                      rho = 1, B = NULL,norm = "l2" , quad = TRUE,iter = 500,e.abs = 1e-3,
                      e.rel= 1e-3,maxiter.B = 50, tol.B = 1e-4,
                      verbose = FALSE){
  
  start.time<- proc.time()[3]
  #Some simple errors to verify!
  
  if(!is.numeric(X) | !is.numeric(Z) | !is.numeric(Y)){
    stop("Data must be a numeric matrix/vector")
  }
  if(!is.vector(Y)){
    stop("Response Y must be a vector")
  }
  if(any(lambdas<0)){
    stop("Lambda values must be positive")
  }
  if(any(alphas<=0) | any(alphas>=1)){
    stop("Alpha values must be between 0 and 1")
  }
  if(nrow(X)!= nrow(Z)| nrow(X)!= length(Y) | nrow(Z)!= length(Y)){
    stop("Mismatched dimensions of data. Make sure nrow(X)=nrow(Z)=length(Y)")
  }
  
  family = family[1]
  if(family == "binomial"){
    if(!all(sort(unique(Y))==c(0,1))){
      stop("Response vector must be a binary vector")
    }
  }
  
  if(!(family== "gaussian" | family=="binomial")){
    stop("Unrecognized family name. Family must be 'gaussian' or 'binomial' ")
  }
  
  if(is.null(B)){
    B<- matrix(0, ncol = ncol(Z)+1,nrow = ncol(X)+1)
    if(quad == FALSE){
      diag(B)[-1] = NaN
    }
  }
  
  if(quad == FALSE) {
    diag(B)[-1] <- NaN
  }
  matD = matE = matF = B
  g1<- rep(0,ncol(B)*nrow(B))
  g2<- rep(0,ncol(B)*nrow(B))
  g3<- rep(0,ncol(B)*nrow(B))
  g.list<- list(matrix(g1,ncol= ncol(B)),matrix(g2,ncol= ncol(B)),
                matrix(g3,ncol= ncol(B)))
  
  
  length.alpha<- length(alphas) 
  length.lambda<- length(lambdas)
  
  fin.result<- vector("list", length(alphas))
  
  cat("Computing w...")
  if(quad){
    w<- generate.w(X=X,Z=Z,quad)
  }else{
    w<- generate.w(X=X,Z=Z,quad)
    w.full<- generate.w(X=X,Z=Z, TRUE)
  }
  
  cat("done.\n")
  num<- length(Y)  
  #SVD decomposition for the main matrix. 
  #Main slow part
  cat("Starting svd...")
  #print(dim(w))
  
  svd.w<- svd(w)
  svd.w$tu<- t(svd.w$u)
  svd.w$tv<- t(svd.w$v)
  
  cat("done.\n")
  
  
  p1<- ncol(X)
  p2<- ncol(Z)
  
  for(i in 1:length.alpha){
    alpha<- alphas[i]
#    #FIND MAX LAMBDA FOR A FIXED ALPHA
#     croot<- 0
#     for(rs in 2:(p1+1)){
#       ind<- (0:p2)*(p1+1)+ rs
#       if(quad){
#         temp<- w[,ind]
#       }else{
#         #print(ind[rs])
#         ind<- ind[-rs]
#         temp<- w.full[,ind]
#       }
#       
#       root<- uniroot.all(max_l, interval = c(1e-5,50), 
#                   grp = temp, y_vec = Y,alpha = alpha, p = p2+1)
#       croot<- max(root,croot)
#     }
#     for(cs in 2:(p2+1)){
#       ind<- ((cs-1)*(p1+1)+1):(cs*(p1+1))
#       if(quad){
#         temp<- w[,ind]
#       }else{
#         #print(ind[rs])
#         ind<- ind[-rs]
#         temp<- w.full[,ind]
#       }
#       root<- uniroot.all(max_l, interval = c(1e-5,50), 
#                          grp = temp, y_vec = Y,alpha = alpha, p = p2+1)
#       croot<- max(root,croot)
#     }
#     
    
    b.array<- vector("list", length.lambda)
    
    if(verbose) cat("Fitting model for alpha =", round(alpha,2) ,
                    "and lambda =", round(lambdas[length.lambda],2), "\n")
    
    if(family=="gaussian"){
      b.array[[length.lambda]]<- estimate.SH(X,Z,Y,w,svd.w = svd.w, 
                                  c(lambdas[length.lambda]*(1-alpha)*sqrt(p2),
                                  lambdas[length.lambda]*(1-alpha)*sqrt(p1),alpha*lambdas[length.lambda]) ,
                                  rho = rho, 
                                  B,matD,matE,matF,g.list,iter = iter,e.abs = e.abs, 
                                  e.rel=e.rel, quad = quad,norm = norm)
      b.array[[length.lambda]]$alpha = alpha
      b.array[[length.lambda]]$lambda = lambdas[length.lambda]
      
    }else{
      b.array[[length.lambda]]<- estimate.logistic(X,Z,Y,w,svd.w = svd.w, 
                                  c(lambdas[length.lambda]*(1-alpha)*sqrt(p2),
                                  lambdas[length.lambda]*(1-alpha)*sqrt(p1),alpha*lambdas[length.lambda]) ,
                                  rho = rho, 
                                  B,matD,matE,matF,g.list,iter = iter,e.abs = e.abs, 
                                  e.rel=e.rel,quad = quad,norm = norm)  
      b.array[[length.lambda]]$alpha = alpha
      b.array[[length.lambda]]$lambda = lambdas[length.lambda]
    }
    
    if(!b.array[[length.lambda]]$conv)warning(paste("The algorithm did not converge for alpha= ",
                                                    alpha," and lambda= " , lambdas[length.lambda],"did not converge"))
    for(lam in (length.lambda-1):1){
      if(verbose) cat("Fitting model for alpha =", round(alpha,2) ,"and lambda =", round(lambdas[lam],2), "\n")
      
      if(family=="gaussian"){
        
        b.array[[lam]]<- estimate.SH(X,Z,Y,w,svd.w = svd.w, 
                                     c(lambdas[lam]*(1-alpha)*sqrt(p2),
                                       lambdas[lam]*(1-alpha)*sqrt(p1),alpha*lambdas[lam]) ,
                                     rho = b.array[[lam+1]]$rho, b.array[[lam+1]]$B,
                                     b.array[[lam+1]]$D,b.array[[lam+1]]$E,b.array[[lam+1]]$F,
                                     b.array[[lam+1]]$glist,iter = iter,e.abs = e.abs, 
                                     e.rel = e.rel,quad = quad,norm = norm)
        b.array[[lam]]$alpha = alpha
        b.array[[lam]]$lambda = lambdas[lam]
        
      }else{
        b.array[[lam]]<- estimate.logistic(X,Z,Y,w,svd.w = svd.w, 
                                     c(lambdas[lam]*(1-alpha)*sqrt(p2),
                                       lambdas[lam]*(1-alpha)*sqrt(p1),alpha*lambdas[lam]) ,
                                     rho = b.array[[lam+1]]$rho, b.array[[lam+1]]$B,
                                     b.array[[lam+1]]$D,b.array[[lam+1]]$E,b.array[[lam+1]]$F,
                                     b.array[[lam+1]]$glist,iter = iter,e.abs = e.abs, 
                                     e.rel = e.rel,quad = quad,norm = norm)
        b.array[[lam]]$alpha = alpha
        b.array[[lam]]$lambda = lambdas[lam]
        
      }
      
      if(!b.array[[lam]]$conv) warning(paste("The algorithm did not converge for alpha= ",
                                                alpha," and lambda= " , lambdas[lam],"did not converge"))
    }
    
    fin.result[[i]]<- b.array
    
  }
  
  end.time<- proc.time()[3]
  run.time<- end.time-start.time
  output.obj<- list("Estimate" = fin.result,"alpha" = alphas,"lambda" = lambdas,
                    "Y.train" = Y, "X.train" = X,"Z.train" = Z,"family" = family,"quad" = quad,
                    "call" = match.call(),"time" = run.time,"w" = w)
  
  class(output.obj)<- "FAMILY"
  return(output.obj)
  
}
########################################################################################
print.FAMILY<- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\nThe function ran for",x$time , "seconds for", length(x$alpha), "alpha values:\n")
  cat("alpha =", round(x$alpha,3) )
  cat("\nand", length(x$lambda),"lambda values:\n")
  cat("lambda =", round(x$lambda,3) ,"\n")
  
  rand.alpha<- sample.int(length(x$alpha), size = 1)
  rand.lambda<- sample.int(length(x$lambda), size = 1)
  cat("\nUse $Estimate to get the fitted model for given alpha and lambda value. For instance
to get the estimated model for alpha = ", round(x$alpha[rand.alpha],3) ," and lambda = ", 
round(x$lambda[rand.lambda],3), " we use fit$Estimate[[",rand.alpha,"]][[",rand.lambda,"]].",
sep = "")
  invisible(x)
}


########################################################################################

#This function extracts the coefficients
#I am not sure how to allow extra parameters in the function using '...'

coef.FAMILY<- function(object,XequalZ = FALSE,Bias.corr = FALSE,...){
  #This function also returns a 2 dimensional list just like the interact.SH function
  coef.result<- object$Estimate
  
  #I wish to take care of the case when $X=Z$ and when $X\not=Z$
  if(!XequalZ){
    for(a in 1:length(object$alpha)){
      for(b in 1:length(object$lambda)){
        #print(c(a,b))
        obj<- object$Estimate[[a]][[b]]
        mainsZ<- which(obj$finB[1,-1]!=0)
        mainsX<- which(obj$finB[-1,1]!=0)
        interacts<- numeric()
        interacANS<- obj$finB[-1,-1]
        coef.interacts<- c()
        for(r in 1:nrow(interacANS)){
          for(c in 1:ncol(interacANS)){
            if(interacANS[r,c]!=0 ){
              interacts<- rbind(interacts,c(r,c))
              coef.interacts<- c(coef.interacts, obj$finB[r+1,c+1])
            }
          }
        }
        
        if(length(interacts!=0)){
          interacts<- as.matrix(cbind(interacts, coef.interacts))
          colnames(interacts)<- c("X","Z","Coef. est")
        }
        
        mainsX<- as.matrix(cbind(mainsX,(obj$finB[-1,1])[mainsX]))
        colnames(mainsX)<- c("X","Coef. est")
        
        mainsZ<- as.matrix(cbind(mainsZ,(obj$finB[1,-1])[mainsZ]))
        colnames(mainsZ)<- c("Z","Coef. est")
        
        intercept = obj$finB[1,1]
        if(Bias.corr){
          myX<- c()
          if(length(mainsX) !=0 ){
            myX<- object$X.train[,mainsX[,1]]
          }
          if(length(mainsZ) !=0 ){
            myX<- cbind(myX, object$Z.train[,mainsZ[,1]])
          }
          
          #If we have interactions seletced
          if(length(interacts)!=0){
            if(length(interacts)==3){
              myX<- cbind(myX, 
                          object$X.train[,interacts[1]]*object$Z.train[,interacts[2]])
            }else{
              myX<- cbind(myX, 
                          (object$X.train[,interacts[,1]])*(object$Z.train[,interacts[,2]]))
            }
            
          }
          
          myY<- object$Y.train
          if(length(myX)==0){
            if(object$family=="binomial"){
              fit.lm <- glm(myY ~ 1, family = "binomial")
            }else{
              fit.lm <- lm(myY ~ 1)
            }
            
          }else{
            if(object$family=="binomial"){
              fit.lm <- glm(myY ~ myX, family = "binomial",maxit = 500)
            }else{
              fit.lm <- lm(myY ~ myX)
            }
            
            coeffs.lm<- fit.lm$coefficients[-1]  
            num.X<- length(mainsX)/2
            num.Z<- length(mainsZ)/2
            num.In<- length(interacts)/3
            if(num.X!=0)indX<- 1:num.X
            if(num.Z!=0)indZ<- (num.X+1):(num.X+num.Z)
            if(num.In!=0) indIn<- (num.X+num.Z+1):(num.X+num.Z+num.In) 
            
            if(num.X!=0) mainsX[,2]<- coeffs.lm[indX]
            if(num.Z!=0) mainsZ[,2]<- coeffs.lm[indZ]
            if(num.In!=0) interacts[,3]<- coeffs.lm[indIn]
            
          }
          
          intercept = fit.lm$coefficients[1]  
        }
        
        coef.result[[a]][[b]]<- list("intercept" = intercept ,"mainsX" = mainsX,
                                     "mainsZ" = mainsZ,"interacts" = interacts,
                                     "alpha" = object$alpha[a], "lambda" = object$lambda[b])
      }
    }
    len.a<- length(coef.result)
    len.l<- length(object$lambda)
    
    names(coef.result)<- paste("alpha = ",object$alpha)
    for(i in 1:len.a){
      names(coef.result[[i]])<- paste("lambda = ",object$lambda)
    }
    
    class(coef.result)<- "coefMISH"
    return(coef.result)
  }else{
    
    #This is the case of $X=Z$.
    for(a in 1:length(object$alpha)){
      for(b in 1:length(object$lambda)){
        obj<- object$Estimate[[a]][[b]]
        mains<- which(obj$finB[-1,1]!=0)
        interacts<- numeric()
        interacANS<- obj$finB[-1,-1]
        coef.interacts<- c()
        for(r in 1:nrow(interacANS)){
          for(c in r:ncol(interacANS)){
            if(interacANS[r,c]!=0 ){
              interacts<- rbind(interacts,c(r,c))
              if(r==c){
                coef.interacts<- c(coef.interacts, obj$finB[r+1,c+1])
              }else{
                coef.interacts<- c(coef.interacts, 2*obj$finB[r+1,c+1])
              }
              
            }
          }
        }
        
        if(length(interacts!=0)){
          interacts<- as.matrix(cbind(interacts, coef.interacts))
          colnames(interacts)<- c("X","Z","Coef. est")
        }
        
        mains<- as.matrix(cbind(mains,(2*obj$finB[-1,1])[mains]))
        colnames(mains)<- c("X","Coef. est")
        
        intercept = obj$finB[1,1]
        if(Bias.corr){
          myX<- c()
          if(length(mains) !=0 ){
            myX<- object$X.train[,mains[,1]]
          }
          
          #If we have interactions seletced
          if(length(interacts)!=0){
            if(length(interacts)==3){
              myX<- cbind(myX, 
                          object$X.train[,interacts[1]]*object$Z.train[,interacts[2]])
            }else{
              myX<- cbind(myX, 
                          (object$X.train[,interacts[,1]])*(object$Z.train[,interacts[,2]]))
            }
            
          }
          
          
          myY<- object$Y.train
          if(length(myX)==0){
            fit.lm <- lm(myY ~ 1)
          }else{
            fit.lm <- lm(myY ~ myX)
            
            coeffs.lm<- fit.lm$coefficients[-1]  
            num.X<- length(mains)/2
            num.In<- length(interacts)/3
            if(num.X!=0)indX<- 1:num.X
            if(num.In!=0) indIn<- (num.X+1):(num.X+num.In) 
            
            if(num.X!=0) mains[,2]<- coeffs.lm[indX]
            if(num.In!=0) interacts[,3]<- coeffs.lm[indIn]
            
          }
          
          intercept = fit.lm$coefficients[1]  
        }
        
        coef.result[[a]][[b]]<- list("intercept" = intercept,"mains" = mains, 
                                     "interacts" = interacts,
                                     "alpha" = object$alpha[a], "lambda" = object$lambda[b])
      }
    }
    
    len.a<- length(coef.result)
    len.l<- length(object$lambda)
    
    names(coef.result)<- paste("alpha  = ",object$alpha)
    for(i in 1:len.a){
      names(coef.result[[i]])<- paste("lambda  = ",object$lambda)
    }
    class(coef.result)<- "coefMISH"
    return(coef.result)
  }
  
}
#########################################################################################
print.coefMISH<- function(x, ...){
  
  cat("\nThis function outputs a list with the same dimensions as the fitted model with ", length(x),"
      different alpha values and ", length(x[[1]])," lambda values.\n")
  
  cat("\nEach element in the 2 dimensional list contains matrices for main effects and interaction terms.\n")
  
  rand.alpha<- sample.int(length(x), size = 1)
  rand.lambda<- sample.int(length(x[[1]]), size = 1)
  cat("\nAs an example, to get the estimated model for alpha = ",x[[rand.alpha]][[rand.lambda]]$alpha ," and lambda = ", 
      x[[rand.alpha]][[rand.lambda]]$lambda, " we use object[[",rand.alpha,"]][[",rand.lambda,"]].",
      sep = "")
  invisible(x)
}

#########################################################################################


predict.FAMILY<- function(object,new.X,new.Z, Bias.corr = FALSE,XequalZ = FALSE,...){
  
  #Generate w for prediction
  w.test <- generate.w(X=new.X,Z=new.Z)
  
  
  n<- nrow(new.X)
  predict.result<- array(0, c(n,length(object$lambda),length(object$alpha)))
  
  #List of dimension names
  if(object$family == "gaussian") yhat.name<- rep("yhat",n)
  if(object$family == "binomial") yhat.name<- rep("phat",n)
  alpha.name<- c()
  lambda.name<- c()
  
  alpha.name<- paste("alpha = ", round(object$alpha,2))
  lambda.name<- paste("lam = ", round(object$lambda,2))
  
  #Obtain coefficient object
  coef.obj<- coef(object, XequalZ)
  
  if(!Bias.corr){
    #This is the simple case when we do not have a bias correction step
    for(a in 1:length(object$alpha)){
      for(b in 1:length(object$lambda)){
        obj<- object$Estimate[[a]][[b]]
        
        predict.result[,b,a] <-  w.test%*%as.vector(obj$finB) 
        if(object$family == "binomial"){
          predict.result[,b,a] <-  1/(1+exp(-predict.result[,b,a]))
        }
      }
    } 
  }else{
    #Here we will have two cases of $X=Z$ and $Xnot= Z$
    #Note that the coef.obj will take care of most of the work for this
    #dichotomy
    
    for(a in 1:length(object$alpha)){
      for(b in 1:length(object$lambda)){
        #print(c(a,b))
        current.obj<- coef.obj[[a]][[b]]
        myX<- rep(1,nrow(object$X.train))
        mytestX<- rep(1,nrow(new.X))

        
        if(XequalZ){
          mains<- current.obj$mains[,1]
          mainsZ<- numeric()
        }else{
          mains<- current.obj$mainsX[,1]
          mainsZ<- current.obj$mainsZ[,1]
        }

        
        #If we have some main effects selected
        if(length(mains)!=0){
          myX<- object$X.train[,mains]  
          mytestX<- new.X[,mains]
        }
        if(length(mainsZ)!=0){
          myX<- cbind(myX,new.Z[,mainsZ])  
          mytestX<- cbind(mytestX,new.Z[,mainsZ])
        }

        #If we have interactions seletced
        if(length(current.obj$interacts)!=0){
          interactions<- (current.obj$interacts[,1:2])
          if(length(interactions)==2){
            interactions<- matrix(interactions,ncol = 2)
          }else{
            myX<- cbind(myX, 
                        (object$X.train[,interactions[,1]])*(object$Z.train[,interactions[,2]]))
            mytestX<- cbind( mytestX, (new.X[,interactions[,1]])*(new.Z[,interactions[,2]]))
          }
          
        }
        
        #Re-define myX, mytestX and myvalX as dataframes
        myX<- as.data.frame(myX)
        colnames(myX)<- 1:(ncol(myX))
        mytestX<- as.data.frame(mytestX)
        colnames(mytestX)<- 1:(ncol(mytestX))
        if(object$family == "gaussian"){
          myY<- object$Y.train
          fit<- lm(myY ~ . , data = cbind(myY,myX))
          p.te<- predict(fit,newdata = mytestX)
          
        }else if(object$family == "binomial"){
          myY<- object$Y.train
          fit<- glm(myY ~ . , data = cbind(myY,myX),family  ="binomial",maxit = 1000)
          p.te<- predict.glm(fit,newdata = mytestX,type = "response")
          
        }
        
        predict.result[,b,a] <-  p.te

      }
    }
    
  }
  
  dimnames(predict.result)<- list(yhat.name,lambda.name,alpha.name)
    
  return(predict.result)
}

