# R programme for the book:
#
# Marginal models
# Models for Dependent, Clustered, and Longitudinal Categorical Data
# www.cmm.st
# 
# Wicher P. Bergsm
# Marcel Croon
# Jacques A. Hagenaars
#
# Springer, Statistics for the Social Sciences
#
# This R-programme is written by Wicher Bergsma and Andries van der Ark, 2009
#
#
##############################################################################
##############################################################################
##############################################################################
#
# List of main functions
#
#   Fitting:
#
#       MarginalModelFit( n, model ) (TO DO: Method="MDI")
#       SampleStatistics( n, coeff ) 
#       ModelStatistics( n, m, model, coeff )
#
#   Matrix functions:
#
#       MarginalMatrix( var, marg, dim )
#       DesignMatrix( var, marg, dim )
#       ConstraintMatrix( var, marg, dim )
#       DirectSum( mat1, mat2 )
# 
#
#   Specification of coefficients and models
#
#       SpecifyCoefficient( name, parameters, rep )  (TO DO: "KendallTau" etc)
#       JoinModels( model1, model2 )
#
#
#   Data formatting:
# 
#       RecordsToFrequencies( rec, dim )
# 
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
# matrix functions: desmat,MarginalMatrix


##############################################################################
# DirectSum (modified from web), takes variable no of args

DirectSum <- function (...){
     p.tr = 0;p.ll = 0;   ## p.tr and p.ll are padding values: originally given as optional args to function
     matlist = list(...);
     nmat = length(matlist);
     m1 = matlist[[1]];
     matlist = if(nmat==1 && is.list(m1)) m1 else matlist # check if list of matrices is given and amend accordingly
     nmat = length(matlist);                              # ,,
     m1 = matlist[[1]];                                   # ,,
     if(nmat==1) return(m1);
     for(i in 2:nmat){ 
        m2 = matlist[[i]];
        topleft <- m1
        topright <- matrix(p.tr, nrow(m1), ncol(m2))
        colnames(topright) <- colnames(m2)
        lowleft <- matrix(p.ll, nrow(m2), ncol(m1))
        lowright <- m2
        m1 = rbind(cbind(topleft, topright), cbind(lowleft, lowright))
     }
     m1
} 

#a=matrix(1:12,3,4);b=matrix(rep(2,6),3,2);DirectSum(a);



##############################################################################
#subsets

#a recursive function returning a list of all subsets of size n of list elements.
#used for function allsubsets
allnsubsets1<-function(x,n,index=1){
   rtn<-list()
   if(length(x)<n){return(NULL)}
   if(length(x)==n){return(list(x))}
   if(length(x)>n){
      #trigger recursion
      for(i in 1:length(x)){
         if(i>=index){
            rtn<-c(rtn,allnsubsets1(x[-i],n,i))
         }
      }   
   }
   return(rtn)
} 
allnsubsets=function(x,n,index=1){rev(allnsubsets1(x,n,index))}

#allnsubsets(c(1,2,3,4),2)


#finds all subsets of a set or a list of sets
allsubsets1 <- function(x){
   if ( data.class(x) == "NULL" ) { subsets <- c() }
   else
   if (data.class(x)=="numeric" || data.class(x)=="character"){
      n=length(x)
      subsets=allnsubsets(x,1)
      if (n>1) for(i in 2:n) subsets=c(subsets,allnsubsets(x,i))
   }
   else{
      subsets=allsubsets1(x[[1]])
      k=length(x)
      if (k>1) for(i in 2:k) subsets=union(subsets,allsubsets1(x[[i]]))
   }
   return(subsets)
} 

allsubsets = function(x){c("NULL",allsubsets1(x))}

#allsubsets(list(c(1,2),c(2,3)))

#allsubsets(c(1,2,3))
#allsubsets(list(c(1,2,3)))
#allsubsets(list(c(1,2),c(2,3),c(1,4)))




##############################################################################
#MarginalMatrix and DesignMatrix


desmat1 = function(var, subvar, dim, coding ){

   #(*lmat2 produces kxk-matrix*)
   lmat2 = function(type,k){
      x = c(1:k) - mean(c(1:k));
      switch( type,
         "Identity"  = diag(k),
         "Nominal"   = rbind(diag(k)[1:k-1,]),
         "Linear"    = rbind(x),
         "Quadratic" = rbind(x,c(1:k)^2),
         "Cubic"     = rbind(x,c(1:k)^2,c(1:k)^3),
         "Quartic"   = rbind(x,c(1:k)^2,c(1:k)^3,c(1:k)^4),
         "Quintic"   = rbind(x,c(1:k)^2,c(1:k)^3,c(1:k)^4,c(1:k)^5)
         )
   }
   
   nvar=length(var);
   coding = if( is.character(coding)&&length(coding)==1 ) rep(coding,length(subvar)) else coding;
   coding = lapply( var, function(x){if(is.element(x,subvar)) coding[subvar==x][[1]] else rbind(rep(1,dim[var==x]))} );

   matlist = list();
   for(i in 1:nvar){
      matlist[[i]] = 
         if      (is.vector(coding[[i]])&&!is.character(coding[[i]])) rbind(coding[[i]])
         else if (is.matrix(coding[[i]])&&!is.character(coding[[i]])) coding[[i]]
         else if (is.element(var[[i]],subvar))                        lmat2(coding[[i]],dim[[i]])
         else                                                         rbind(rep(1,dim[[i]]))
   }
   mat=matlist[[1]];
   if(nvar>1) for (i in 2:nvar){mat = kronecker(mat,matlist[[i]])}

   return(t(mat))
}

#desmat1(c(1,2),c(1,2),c(3,3),c("Nominal","Quintic"))
#desmat1(c(1,2),c(1,2),c(3,3),list("Nominal",rbind(c(1,2,3))))
#desmat1(c(1,2),c(1,2),c(3,3),c("Nominal","Quintic"))
#desmat1(c("A","B"),c("A","B"),c(3,3),c("Nominal","Quintic"))
#desmat1(c("A","B"),c(),c(2,2),c("Nominal"))
#desmat1(c("A"),c("A"),c(2),c("Nominal"))


DesignMatrix = function(var, suffconfigs, dim, SubsetCoding="Automatic", MakeSubsets=TRUE){

   suffconfigs = if(is.list(suffconfigs)) suffconfigs else list(suffconfigs)
   marglist = if(MakeSubsets) allsubsets(suffconfigs) else suffconfigs;
   nmarg=length(marglist);
   
   #put in standard form
   q = SubsetCoding[[1]];
   SubsetCoding = if( is.character(q) || is.numeric(q) ) list(SubsetCoding) else SubsetCoding

   nsubs = length(SubsetCoding);

   #coding gives specification for each marginal
   coding = 
      if( isTRUE(SubsetCoding=="Automatic") ) rep("Nominal",nmarg)   
      else if( is.character(SubsetCoding[[1]]) && length(SubsetCoding)==1 && length(SubsetCoding[[1]])==1) rep(SubsetCoding,nmarg) 
      else lapply(marglist, function(x){ 
         for(i in 1:nsubs) {
            ss1 = SubsetCoding[[i]][[1]];
            if(length(x)==length(ss1)) if(isTRUE(prod(x==ss1)==1)) return(SubsetCoding[[i]][[2]])
            }
         return("Nominal")
         }  )
   deslist = list();
   for(i in 1:nmarg) deslist[[i]] = desmat1(var,marglist[[i]],dim,coding[[i]])
   des = deslist[[1]];
   if(nmarg>1) for(i in 2:nmarg) des = cbind(des,deslist[[i]])
   return(des)
}
#DesignMatrix(c("P", "R"), list(c("P","R")), c(3,3), SubsetCoding = list(c("P", "R"),list("Linear", rbind(c(1,1,0),c(1,1,0),c(0,0,1)))))
#DesignMatrix(c("P", "R"), list(c("P","R")), c(3,3), SubsetCoding = list(c("P", "R"),list("Linear", rbind(c(1,1,0),c(1,1,0),c(0,0,1)))))
#DesignMatrix(c("A","B"),list(c("A"),c("B")),c(2,2))
#DesignMatrix(c(1,2,3),list(c(1,2),c(2,3)),c(3,3,3),SubsetCoding=list(list(c(1,2),c("Linear","Linear")),list(c(1,2),c("Linear","Linear"))))
#t(DesignMatrix(c(1,2),c(1,2),c(3,3),SubsetCoding="Identity",MakeSubsets=FALSE))



EMarginalMatrix <- function(X, marg=NULL, var=var.default, dim=NULL){
  var.default <- 1:ncol(X)
  X <- X[,var]
  if(class(marg)!="list")marg <- list(marg)
  match.c <- function(strng,positions,criteria) all(unlist(strsplit(strng,NULL))[positions] == criteria)
  # dim.default <- apply(cbind(X,NA),2, unique)[-(ncol(X)+1)]
  # dim <- lapply(dim, sort)
  
  n <- table(apply(X,1,paste, collapse=""))
  lab.n <- names(n)
  D <- NULL; i <- 1; j <- 1
  for (i in 1:length(marg)){
    if(is.null(marg[[i]])) D <- rbind(D, rep(1,length(lab.n)))
    else{
      c.values <- expand.grid(dim[marg[[i]]])
      if(ncol(c.values) > 1) c.values <- c.values[,ncol(c.values):1]
      c.positions <- marg[[i]]
      for (j in 1:nrow(c.values)) D <- rbind(D,sapply(lab.n, match.c, c.positions, c.values[j,])*1)
    }
  }
  return(D)
}

MarginalMatrix = function(var,marg,dim,SubsetCoding="Identity",SelectCells="All"){
   if( isTRUE(SelectCells=="All") )
      t(DesignMatrix(var,marg,dim,SubsetCoding=SubsetCoding,MakeSubsets=FALSE))
   else
      EMarginalMatrix(SelectCells,marg,var,dim)
}
#MarginalMatrix(c(1,2),list(c(1),c(2)),c(3,3))

new.null <- function (M) {
#   library(MASS)
   tmp <- qr(M)
   set <- if (tmp$rank == 0) 1:ncol(M) else -(1:tmp$rank)
   qr.Q(tmp, complete = TRUE)[, set, drop = FALSE]
}


#returns hierarchical model constraint matrix, ie nullspace of design matrix
ConstraintMatrix = function(var,suffconfigs,dim, SubsetCoding="Automatic"){
   return(t(new.null(DesignMatrix(var,suffconfigs,dim,SubsetCoding))))
   }

#ConstraintMatrix(c(1,2),list(c(1),c(2)),c(2,2))

#compute the matrix for obtaining the beta parameters under different coding schemes

betamat = function(var, subvar, dim, coding){

   #lmat1 produces kx1-vector, lmat2 produces kxk-matrix

   polmat = function(x,k){
      powmat = rbind(rep(1,k));
      if(k>1) for(i in 1:(k-1)) powmat = rbind(powmat, x^i - mean(x^i))
      invmat = solve(t(powmat))
      return(rbind(invmat[2:k,],rep(0,k)))     }
   
   lmat1 = function(type2, k){
      x = if(length(type2)>1) type2[[2]] else NULL
      type = if(length(type2)>1) type2[[1]] else type2 #second arg is constant for dummy coding, vector for polynomial
      switch( type,
         "Effect"     = rbind(rep(1/k,k)),
         "Simple"     = rbind(rep(1,k)),
         "Dummy"      = {if(length(x)==0) x=1; rbind(diag(k)[x,])}, 
         "Polynomial" = rbind(rep(1,k))
      ) }  

   lmat2 = function(type2,k){
      x = if(length(type2)>1) type2[[2]] else NULL
      type = if(length(type2)>1) type2[[1]] else type2 #second arg is constant for dummy coding, vector for polynomial
      switch( type,
         "Effect"     = diag(k)-1/k,
         "Simple"     = diag(k),
         "Dummy"      = {if(length(x)==0) x=1; diag(k) - t(matrix(rep(diag(k)[x,],k),k,k))},
         "Polynomial" = {if(length(x)==0) x=c(1:k)-(k+1)/2; polmat(x,k)}
         ) }
         
   nvar = length(var);
   coding = if( is.character(coding)&&length(coding)==1 ) as.list(rep(coding,nvar)) else coding;

   matlist = list();
   for( i in 1:nvar ){matlist[[i]] = if(is.element(var[[i]],subvar)) lmat2(coding[[i]],dim[[i]]) else lmat1(coding[[i]],dim[[i]]) }

   mat = matlist[[1]] 
   if(nvar==1) mat else for(i in 2:nvar) mat = kronecker(mat,matlist[[i]])

   return(mat)

   }

#betamat(c(1,2),c(1,2),c(2,2),list(list("Dummy",2),"Effect"))
#betamat(c(1,2),c(1,2),c(3,3),"Effect")
#betamat(c(1,2),c(1,2),c(2,2),list(list("Polynomial",c(-2,2)),"Polynomial"))


#   
#   lmat2["Polynomial", k_] := 
#   lmat2[{"Polynomial", Range[k] - (k + 1)/2}, k];

























###############################################################################################
###############################################################################################
###############################################################################################
# specification of coefficients


#probabilities in honogeneous form
spec.prob <- function(dim){
  k = prod(dim);
  id <- diag(rep(1,k))
  one <- t(rep(1,k))
  m1 <- rbind(id,one)
  m2 <- t(rbind(id,-one))
  matlist <- list(id,m2,m1)
  funlist <- list("exp","log","identity")
  coeff <- list(matlist, funlist)
  return(coeff)
}

#log probabilities in honogeneous form
spec.logprob <- function(dim){
  k = prod(dim);
  id <- diag(rep(1,k))
  one <- t(rep(1,k))
  m1 <- rbind(id,one)
  m2 <- t(rbind(id,-one))
  matlist <- list(m2,m1)
  funlist <- list("log","identity")
  coeff <- list(matlist, funlist)
  return(coeff)
}


spec.condprob <- function(z){  # z=(var,condvar,dim)
  var = z[[1]];
  condvar = z[[2]];
  dim = z[[3]];
  ident = diag(prod(dim));
  at1 = MarginalMatrix(var,condvar,dim);
  at1 = t(at1) %*% at1;
  at1 = rbind(ident,at1);
  at2 = cbind(ident,-ident);
  at3 = ident;
  matlist = list(at3, at2, at1);
  funlist = list("exp", "log", "identity");
  list(matlist, funlist)
}
  
#spec.condprob(list(c(1,2),c(1),c(2,2)))



spec.mean <- function(scores){
   at <- t(matrix(scores))
   prob <- spec.prob(length(scores))
   matlist <- c(list(at), prob[[1]])
   funlist <- c(list("identity"), prob[[2]])
   coeff <- list(matlist, funlist)
   return(coeff)
}


spec.ginimeandifference <- function(scores){
   k <- length(scores)
   q=diag(k)
   m1=matrix(nrow=k^2,ncol=k)
   for(i in 1:k) for(j in 1:k) m1[(i-1)*k + j,] = q[i,] + q[j,]
   m2=matrix(nrow=1,ncol=k^2)
   for(i in 1:k) for(j in 1:k) m2[,(i-1)*k + j] = abs(scores[i] - scores[j])
   prob <- spec.prob(k)
   matlist <- c( list( m2, m1 ), prob[[1]] )
   funlist <- c(list( "exp", "log" ), prob[[2]])
   coeff <- list(matlist, funlist)
   return(coeff)
}

spec.entropy <- function(k){
   at <- t(matrix(rep(1,k)))
   prob <- spec.prob(k)
   matlist <- c(list(at), prob[[1]])
   funlist <- c(list("xlogx"), prob[[2]])
   coeff <- list(matlist, funlist)
   return(coeff)
}


spec.diversityindex <- function(k){
   at <- t(matrix(rep(1,k)))
   prob <- spec.prob(k)
   matlist <- c(list(at), prob[[1]])
   funlist <- c(list("xbarx"), prob[[2]])
   coeff <- list(matlist, funlist)
   return(coeff)
}


spec.variance <- function(scores){
   newscores <- scores - min(scores) + 1
   m1 <- rbind(newscores^2, newscores, rep(1,length(scores)))
   m2 <- rbind(c(1, 0, -1),c(0, 2, -2))
   m3 <- t(matrix(c(1,-1)))
   matlist <- list(m3, m2, m1)
   funlist <- list("exp", "log", "identity")
   coeff <- list(matlist, funlist)
   return(coeff)
}

spec.standarddeviation <- function(scores){
   varspec <- spec.variance(scores)
   m1 <- t(matrix(c(1)))
   matlist <- c(list(m1), varspec[[1]])
   funlist <- c(list("sqrt"), varspec[[2]])
   coeff <- list(matlist, funlist)
   return(coeff)
}

#coeff=spec.standarddeviation(c(1,2))
#gfunction(c(2,2),coeff)

spec.covariance <- function(xyscores){
   xscores <- xyscores[[1]]
   xscores <- xscores - min(xscores) + 1
   yscores <- xyscores[[2]]
   yscores <- yscores - min(yscores) + 1
   k <- length( xscores )
   l <- length( yscores )
   m1a <- MarginalMatrix( c(1,2), list(c(1,2),c(1),c(2)), c(k,l) )
   m1b1 <- as.vector( as.vector(xscores) %o% as.vector(yscores) )
   m1b1 <- rbind( rep( 1, k*l ), m1b1 )
   m1b2 <- DirectSum( rbind(xscores), rbind(yscores) )
   m1b <- DirectSum( m1b1, m1b2 )
   m1 <- m1b %*% m1a
   m2 <- rbind(c(1, 1, 0, 0),c(0, 0, 1, 1))
   m3 <- t(matrix(c(1,-1)))
   matlist <- list(m3, m2, m1)
   funlist <- list("exp", "log", "identity")
   list(matlist, funlist)
}

   
spec.correlation <- function(xyscores){
   xscores <- xyscores[[1]]
   xscores <- xscores - min(xscores) + 1
   yscores <- xyscores[[2]]
   yscores <- yscores - min(yscores) + 1
   k <- length( xscores )
   l <- length( yscores )
   m1a <- MarginalMatrix( c(1,2), list(c(1,2),c(1),c(2)), c(k,l) )
   m1b1 <- as.vector( as.vector(xscores) %o% as.vector(yscores) )
   m1b1 <- rbind( rep( 1, k*l ), m1b1 )
   m1b2 <- DirectSum( rbind(xscores,xscores^2), rbind(yscores,yscores^2) )
   m1b <- DirectSum( m1b1, m1b2 )
   m1 <- m1b %*% m1a
   m2 <- rbind(c(1, 1, 0, 0, 0, 0), c(0, 0, 1, 0, 1, 0), c(0, 0, 0, 1, 0, 0), c(0, 0, 2, 0, 0, 0), c(0, 0, 0, 0, 0, 1), c(0, 0, 0, 0, 2, 0) )
   m3 <- rbind( c(1, 0, 0, 0, 0, 0), c(0, 1, 0, 0, 0, 0), c(0, 0, 1, 0, -1, 0), c(0, 0, 0, 1, 0, -1) )
   m4 <- rbind( c(1, 0, -.5, -.5), c(0, 1, -.5, -.5) )
   m5 <- rbind( c(1, -1) )
   matlist <- list( m5, m4, m3, m2, m1 )
   funlist <- list( "exp", "log", "exp", "log", "identity" )
   list(matlist, funlist)
}

#SpecifyCoefficient("Correlation",list(c(1,2),c(1,2)))

spec.diffprop <- function(c,Response="Columns",DropLast=FALSE){
   id1 <- diag(c)
   id2 <- diag(2*c)
   
   m1b = if(Response == "Columns") MarginalMatrix(c(1, 2), c(1), c(2,c)) else MarginalMatrix(c(1, 2), c(2), c(c,2))
   m1b = t(m1b) %*% m1b;
   m1 = rbind(id2, m1b);

   m2 = cbind(id2, -id2);
   m3 = cbind(id1, -id1);
   m3 = if(DropLast) m3[1:dim(m3)[1]-1,] else m3
   list( list(m3, m2, m1), list("exp", "log", "identity"))
}

#spec.diffprop(2,Response="Rows",DropLast=TRUE)

spec.goodmankruskaltau <- function(rc,Response="Columns"){
   r=rc[[1]];
   c=rc[[2]];

   id = diag(r*c);
   idc = diag(c);
  
   m1a = MarginalMatrix(c(1,2),c(1),c(r,c))
   m1b = MarginalMatrix(c(1,2),c(2),c(r,c))
   m1 = if(Response=="Columns") rbind(id,t(m1a) %*% m1a,m1b) else rbind(id,t(m1b) %*% m1b,m1a)

   z1=matrix(0,r*c,c);
   z2=matrix(0,r*c,r*c);
   m2 = rbind(cbind(2*id, -id, z1), cbind(id, z2, z1), cbind(t(z1), t(z1), 2*idc));

   m3a = matrix(1,1,r*c);
   m3b = matrix(1,1,c);
   zero= matrix(0,1,r*c);
   m3 = rbind(cbind(m3a, zero, -m3b), cbind(zero, m3a, -m3b));

   m4 = matrix(c(1,-1),nrow=1,ncol=2);
   m5 = matrix(c(1),nrow=1,ncol=1);

   funlist = list("exp", "log", "exp", "log", "identity");
   matlist = list(m5, m4, m3, m2, m1);

   prob <- spec.prob(r*c)
   matlist <- c(matlist, prob[[1]])
   funlist <- c(funlist, prob[[2]])
   coeff <- list(matlist, funlist)
   return(coeff)
}

#coeff=spec.goodmankruskaltau(3,2)
#gfunction(c(4,2,6,4,6,7),coeff)


spec.uncertaintycoefficient <- function(rc){
   r=rc[[1]];
   c=rc[[2]];

   m1 = MarginalMatrix(c(1,2),list(c(1, 2),c(1),c(2)), c(r,c));
   m2rc = matrix(1,1,r*c);
   m2r = matrix(1,1,r);
   m2c = matrix(1,1,c);
   zrc= matrix(0,1,r*c);
   zr = matrix(0,1,r);
   m2 = rbind(cbind(m2rc, -m2r, -m2c), cbind(zrc, zr, -m2c));
   m3 = matrix(c(1,-1),nrow=1,ncol=2);
   m4 = matrix(c(1),nrow=1,ncol=1);

   funlist = list("exp","log","xlogx","identity");
   matlist = list(m4, m3, m2, m1);

   prob <- spec.prob(r*c)
   matlist <- c(matlist, prob[[1]])
   funlist <- c(funlist, prob[[2]])
   coeff <- list(matlist, funlist)
   return(coeff)
}


spec.cohenkappa <- function(k){
   m1a = c(diag(k));
   m1b = MarginalMatrix(c(1,2), list(c(),c(1),c(2)), c(k,k));
   m1one = matrix(1,1,k*k);
   m1 = rbind(m1a, m1b, m1one);
   m2a = DirectSum(DirectSum(matrix(1,1,1),matrix(1,1,1)),MarginalMatrix(c(1,2),c(2),c(2,k)));
   m2b = apply(m2a,1,sum);
   m2 = cbind(m2a,-m2b);
   m3 = DirectSum(matrix(1,1,1), rbind(c(0,rep(1, k)), c(1,rep(-1, k))));
   m4 = rbind(c(1, 0, -1), c(0, 1, -1));
   m5 = matrix(c(1,-1),nrow=1,ncol=2);
   list(list(m5, m4, m3, m2, m1),list("exp","log","exp","log","identity"))
}

#coeff=spec.cohenkappa(2)
#gfunction(c(4,2,6,4),coeff)

all.patterns <- function(J,m){
  grid <- list()
  j <- 0;
  p <- m^J
  for (j in 1:J){
    grid <- c(grid, j)
    grid[[j]] <- 0:(m-1)
  }
  X <- t(expand.grid(grid))
  dimnames(X) <- NULL
  return(X[J:1,])
}


spec.cronbachalpha <- function(arg, data){
# arg = list(
#         list(items alpha 1, items alpha 2, ..., items alpha k),
#         list(    g alpha 1,     g alpha 2, ...,     g alpha k),
#         list(group alpha 1, group alpha 2, ..., group alpha k)
#       )
# 
  eps <- 1e-5
# Functions for matrices filled with ones (U) and zeroes (O)
  U <- function(a,b) matrix(1,a,b)
  O <- function(a,b) matrix(0,a,b)
# Function converts string to integer
  string2integer <- function(s) as.numeric(unlist(strsplit(s,NULL)))
# Function modified DirectSum
  direct.sum <- function(L){
    if (length(L)==1) return(L[[1]])
    C <- L[[1]]
    for (i in 2:length(L)) C <- DirectSum(C,L[[i]])
    return(C)
  }

  if(is.null(arg)) stop("'arg' not specified")
  if(!is.list(arg)) stop("'arg' is not a list")
  if(!is.list(arg[[1]])) stop("First element of 'arg' is not a list")
  method <- ifelse(is.null(data),1,2)  
  rep <- length(arg[[1]])
  mats <- list()
  for (i in 1:5)  mats[[i]] <- list()
  for (i in 1:rep){
    items <- arg[[1]][[i]]
    J <- length(items)
    g <- arg[[2]][[i]]
    group.var <- NULL
    if(length(arg)==3) group.var <- 1
    if (method==1){
      R <- all.patterns(J,g)
    } else {
      if(is.null(group.var)) Y <- data[,items] else Y <- data[data[,group.var]==i,items]
      if(any(apply(Y,2,var) < eps))stop("One or more variables have zero variance")  
      lab.n <- matrix(names(table(apply(Y,1,paste, collapse=""))))
      R <- apply(lab.n,1,string2integer)
    }
    S <- U(1,J) %*% R
    C5 <- rbind(R,S,R^2,S^2,U(1,ncol(R)))
    C4 <- rbind(
            cbind(O(J+1,J+1),  diag(J+1), U(J+1,1)),
            cbind(diag(2,J+1), O(J+1,J+1),O(J+1,1))
          )
    C3 <- rbind(
           cbind(U(1,J), 0,      -U(1,J), 0),
           cbind(O(2,J), U(2,1),  O(2,J), -U(2,1))
          )
    C2 <- matrix(c(0,1,1,-1,-1,0),nrow=2)
    C1 <- matrix(c(J/(J-1),-J/(J-1)),nrow=1)
    mats[[1]][[i]] <- C1
    mats[[2]][[i]] <- C2
    mats[[3]][[i]] <- C3
    mats[[4]][[i]] <- C4
    mats[[5]][[i]] <- C5
  }
  C1 <- direct.sum(mats[[1]])   
  C2 <- direct.sum(mats[[2]])   
  C3 <- direct.sum(mats[[3]])   
  C4 <- direct.sum(mats[[4]])   
  C5 <- direct.sum(mats[[5]])   
  matlist = list(C1,C2,C3,C4,C5)
  funlist = list("exp","log","exp","log","identity")
  return(list(matlist,funlist))
}

spec.mokkenh <- function(arg, data){
# arg = list(
#         list(items H 1, items H 2, ..., items H k),
#         list(    g H 1,     g H 2, ...,     g H k),
#         list(group H 1, group H 2, ..., group H k)
#       )
# 
  eps <- 1e-5
# Functions for matrices filled with ones (U) and zeroes (O)
  U <- function(a,b) matrix(1,a,b)
  O <- function(a,b) matrix(0,a,b)
# Function converts string to integer
  string2integer <- function(s) as.numeric(unlist(strsplit(s,NULL)))
# Function modified DirectSum
  direct.sum <- function(L){
    if (length(L)==1) return(L[[1]])
    C <- L[[1]]
    for (i in 2:length(L)) C <- DirectSum(C,L[[i]])
    return(C)
  }
# Function prp
  prp <- function(A,B){
    if(nrow(A)!=nrow(B) | ncol(A)!=ncol(B)) stop('A and B have unequal rows or columns')
    n <- nrow(A)
    C <- matrix(nrow=n*(n-1)/2,ncol=ncol(A))
    i <- 0; j <- 0
    for(i in (1:(n-1))) for(j in ((i+1):n)){
      k <- (i-1)*n - sum(0:i) + j
      C[k,] <- A[i,] * B[j,]
    }
    return(C)
  }

  if(is.null(arg)) stop("'arg' not specified")
  if(!is.list(arg)) stop("'arg' is not a list")
  if(!is.list(arg[[1]])) stop("First element of 'arg' is not a list")
  method <- ifelse(is.null(data),1,2)  
  rep <- length(arg[[1]])
  mats <- list()
  for (i in 1:5)  mats[[i]] <- list()
  for (i in 1:rep){
    items <- arg[[1]][[i]]
    J <- length(items)
    g <- arg[[2]][[i]]
    group.var <- NULL
    if(length(arg)==3) group.var <- 1
    if (method==1){
      R <- all.patterns(J,g)
    } else {
      if(is.null(group.var)) Y <- data[,items] else Y <- data[data[,group.var]==i,items]
      if(any(apply(Y,2,var) < eps))stop("One or more variables have zero variance")  
      lab.n <- matrix(names(table(apply(Y,1,paste, collapse=""))))
      R <- apply(lab.n,1,string2integer)
    }
    p <- J*(J-1)/2
  
    C5 <- rbind(U(1,ncol(R)),1-R,R,prp(1-R,R)) 
    Q1 <- O(p,2*J)
    for (ii in (1:(J-1))) for (jj in (ii+1):J) Q1[J*(ii-1) - ii*(ii+1)/2 + jj,c(ii,J+jj)] <- 1
    C4 <- DirectSum(U(1,1),Q1,diag(p))
    C3 <- DirectSum(U(2,1),U(1,p),U(1,p))
    C2 <- matrix(c(1,-1,0,0, 0, 1,-1,1),2,4,byrow=TRUE)
    C1 <- matrix(c(1,-1),1,2)
    mats[[1]][[i]] <- C1
    mats[[2]][[i]] <- C2
    mats[[3]][[i]] <- C3
    mats[[4]][[i]] <- C4
    mats[[5]][[i]] <- C5
  }
  C1 <- direct.sum(mats[[1]])   
  C2 <- direct.sum(mats[[2]])   
  C3 <- direct.sum(mats[[3]])   
  C4 <- direct.sum(mats[[4]])   
  C5 <- direct.sum(mats[[5]])   
  matlist = list(C1,C2,C3,C4,C5)
  funlist = list("exp","log","exp","log","identity")
  return(list(matlist,funlist))
}

spec.allmokkenhj <- function(arg, data){
# arg = list(
#         list(items Hj 1, items Hj 2, ..., items Hj k),
#         list(    g Hj 1,     g Hj 2, ...,     g Hj k), # SHOULD BE 2
#         list(group Hj 1, group Hj 2, ..., group Hj k)
#       )
# 
  eps <- 1e-5
# Functions for matrices filled with ones (U) and zeroes (O)
  U <- function(a,b) matrix(1,a,b)
  O <- function(a,b) matrix(0,a,b)
# Function converts string to integer
  string2integer <- function(s) as.numeric(unlist(strsplit(s,NULL)))
# Function modified DirectSum
  direct.sum <- function(L){
    if (length(L)==1) return(L[[1]])
    C <- L[[1]]
    for (i in 2:length(L)) C <- DirectSum(C,L[[i]])
    return(C)
  }
# Function prp
  prp <- function(A,B){
    if(nrow(A)!=nrow(B) | ncol(A)!=ncol(B)) stop('A and B have unequal rows or columns')
    n <- nrow(A)
    C <- matrix(nrow=n*(n-1)/2,ncol=ncol(A))
    i <- 0; j <- 0
    for(i in (1:(n-1))) for(j in ((i+1):n)){
      k <- (i-1)*n - sum(0:i) + j
      C[k,] <- A[i,] * B[j,]
    }
    return(C)
  }

  if(is.null(arg)) stop("'arg' not specified")
  if(!is.list(arg)) stop("'arg' is not a list")
  if(!is.list(arg[[1]])) stop("First element of 'arg' is not a list")
  method <- ifelse(is.null(data),1,2)  
  rep <- length(arg[[1]])
  mats <- list()
  for (i in 1:5)  mats[[i]] <- list()
  for (i in 1:rep){
    items <- arg[[1]][[i]]
    J <- length(items)
    g <- arg[[2]][[i]]
    group.var <- NULL
    if(length(arg)==3) group.var <- 1
    if (method==1){
      R <- all.patterns(J,g)
    } else {
      if(is.null(group.var)) Y <- data[,items] else Y <- data[data[,group.var]==i,items]
      if(any(apply(Y,2,var) < eps))stop("One or more variables have zero variance")  
      lab.n <- matrix(names(table(apply(Y,1,paste, collapse=""))))
      R <- apply(lab.n,1,string2integer)
    }
    p <- J*(J-1)/2
  
    C5 <- rbind(U(1,ncol(R)),1-R,R,prp(1-R,R)) 
    Q1 <- O(p,2*J)
    for (ii in (1:(J-1))) for (jj in (ii+1):J) Q1[J*(ii-1) - ii*(ii+1)/2 + jj,c(ii,J+jj)] <- 1
    C4 <- DirectSum(U(1,1),Q1,diag(p))
    Q2 <- t(Q1[,1:J] + Q1[,(J+1):(2*J)])

    C3 <- DirectSum(U(2,1),Q2,Q2)
    C2 <-  rbind(
        cbind(matrix(c(1,-1),1,2),O(1,2*J)),
        cbind(U(J,1),O(J,1),-diag(J),diag(J))
      )
    C1 <- cbind(U(J,1),diag(-1,J))
    mats[[1]][[i]] <- C1
    mats[[2]][[i]] <- C2
    mats[[3]][[i]] <- C3
    mats[[4]][[i]] <- C4
    mats[[5]][[i]] <- C5
  }
  C1 <- direct.sum(mats[[1]])   
  C2 <- direct.sum(mats[[2]])   
  C3 <- direct.sum(mats[[3]])   
  C4 <- direct.sum(mats[[4]])   
  C5 <- direct.sum(mats[[5]])   
  matlist = list(C1,C2,C3,C4,C5)
  funlist = list("exp","log","exp","log","identity")
  return(list(matlist,funlist))
}


spec.hij <- function(){
  cat("For computing Loevinger's scalability coefficients it is assumed that the items are ordered by there p-value;",fill=TRUE)
  cat("item 1 is the easiest/most popular item",fill=TRUE)
  R <- matrix(c(0,0,0,1,1,0,1,1),2,4)
  C3 <- rbind(matrix(1,2,4),1-R,R,matrix(c(0,1,0,0),1,4))
  # C3 . n = (n+++, n+++, n0++, n+0+, n++0, n1++, n+1+, n++1, n01+, n0+1, n+01)'
  C2 <- matrix(c(1,1,-1,0,0,-1,0,0,0,0,0,-1,0,1),2,7)
  # exp(C2. log(C1 .n)) = 1, 1-H12, 1-H13, 1-H14
  C1 <- matrix(c(1, -1),1,2)
  # C1 %*% exp(C2 %*% log(C3 %*% n))
  matlist = list(C1,C2,C3)
  funlist = list("exp","log","identity")
  return(list(matlist,funlist))
}


spec.concdiscprobs <- function(rc){
   r=rc[[1]];
   c=rc[[2]];

   m1a=matrix(0,(r-1)*(c-1),r*c);
   m1b=matrix(0,(r-1)*(c-1),r*c);
   m1c=matrix(0,(r-1)*(c-1),r*c);
   m1d=matrix(0,(r-1)*(c-1),r*c);
   for(i in 1:(r-1)) for(j in 1:(c-1)) for(k in (i+1):r) for(l in (j+1):c) {m1a[[(i-1)*(c-1)+j,(k-1)*c+l]]=1}
   for(i in 1:(r-1)) for(j in 1:(c-1)) for(k in (i+1):r) for(l in (j+1):c) {m1b[[(i-1)*(c-1)+j,(i-1)*c+j]]=1}
   for(i in 1:(r-1)) for(j in 2:c) for(k in (i+1):r) for(l in 1:(j-1)) {m1c[[(i-1)*(c-1)+j-1,(k-1)*c+l]]=1}
   for(i in 1:(r-1)) for(j in 2:c) for(k in (i+1):r) for(l in 1:(j-1)) {m1d[[(i-1)*(c-1)+j-1,(i-1)*c+j]]=1}
   m1 = rbind(m1a, m1b, m1c, m1d);
   m2 = diag((r-1)*(c-1));
   m2 = cbind(m2, m2);
   m2 = DirectSum(m2, m2);
   m3 = rbind(4 * rep(1,(r-1)*(c-1)));
   m3 = DirectSum(m3, m3);
   funlist = list("exp","log","identity");
   matlist = list(m3, m2, m1);

   prob <- spec.prob(r*c)
   matlist <- c(matlist, prob[[1]])
   funlist <- c(funlist, prob[[2]])
   coeff <- list(matlist, funlist)
   return(coeff)
}

#coeff=spec.concdiscprobs(2,3)
#gfunction(c(4,2,6,4,6,7),coeff)


spec.kendalltau <- function(rc){
   concdiscspec <- spec.concdiscprobs(rc)
   m1 = matrix(c(1,-1),nrow=1,ncol=2);
   matlist <- c(list(m1), concdiscspec[[1]])
   funlist <- c(list("identity"), concdiscspec[[2]])
   coeff <- list(matlist, funlist)
   return(coeff)
}

spec.gamma <- function(rc){

   m1 = rbind(c(1, 0), c(0, 1), c(1, 1));
   m2 = rbind(c(1, 0, -1), c(0, 1, -1));
   m3 = rbind(c(1, -1));

   funlist = list("exp","log","identity");
   matlist = list(m3, m2, m1);

   concdiscspec <- spec.concdiscprobs(rc)
   matlist <- c(matlist, concdiscspec[[1]])
   funlist <- c(funlist, concdiscspec[[2]])
   coeff <- list(matlist, funlist)
   return(coeff)
}

#coeff=spec.gamma(2,3)
#gfunction(c(4,2,6,4,6,7),coeff)

spec.loglinearpars = function(z){
   if(length(z)==2&&!is.numeric(z)) 
      {dim = z[[1]]; coding = z[[2]]}
   else 
      {dim = z; coding = "Effect"}
   var = c(1:length(dim));
   mat = betamat(var, var, dim, coding)
   list(list(mat),list("log"))
}

#spec.loglinearpars(list(c(2,2),"Dummy"))
#spec.loglinearpars(c(2,2))


spec.logOR <- function(z){
   list(list(rbind(c(1,-1,-1,1))),list(c("log")))
}


#unitvec(5,2) yields c(0,1,0,0,0)
unitvec <- function(n,k){ a=rep(0,n); a[k]=1; a }

spec.ordinal.loc.L <- function(z){

   nvar = z[[1]];
   ncat = z[[2]];
   
   spec0 = spec.condprob(list(c(1,2),c(1),c(nvar,ncat)))
   matlist0 = spec0[[1]];
   funlist0 = spec0[[2]];

   m1a = list(); k=0;
   for(a in 1:nvar) for(b in 1:nvar) for(i in 1:(ncat-1)) for(j in (i+1):(ncat)) 
      {k=k+1; m1a[[k]]=unitvec(nvar*ncat,(a-1)*ncat+i)+unitvec(nvar*ncat,(b-1)*ncat+j)}
   m1a = do.call(rbind,m1a);
   m1b = list(); k=0;
   for(a in 1:nvar) for(b in 1:nvar) for(j in 1:(ncat-1)) for(i in (j+1):(ncat)) 
      {k=k+1; m1b[[k]]=unitvec(nvar*ncat,(a-1)*ncat+i)+unitvec(nvar*ncat,(b-1)*ncat+j)}
   m1b = do.call(rbind,m1b);
   m1 = rbind(m1a,m1b);
   one = t(as.matrix(rep(1,ncat*(ncat-1)/2)));
   m2 = one; for(i in 1:(nvar^2-1)) m2 = DirectSum(m2,one)
   m2 = DirectSum(m2,m2)
   m3 =diag(dim(m2)[[1]]/2)
   m3 = cbind(m3,-m3)
   
   matlist1 = list(m3,m2,m1);
   funlist1 = list("log","exp","log");

   list(append(matlist1,matlist0),append(funlist1,funlist0))
}

spec.ordinal.loc.A <- function(z){ #z=c(nvar,ncat)
   spec0 = spec.ordinal.loc.L(z);
   spec0[[2]][[1]] = "identity";  
   spec0
}
#mat=spec.ordinal.loc.A(c(2,3))


spec.ordinal.disp.L <- function(z){

   nvar = z[[1]];
   ncat = z[[2]];
   if(ncat<3){print("Number of categories must be at least 3");stop()};
   
   spec0 = spec.condprob(list(c(1,2),c(1),c(nvar,ncat)))
   matlist0 = spec0[[1]];
   funlist0 = spec0[[2]];

   m1a = list(); z=0;
   for(a in 1:nvar) for(b in 1:nvar) for(i in 1:(ncat-2)) for(j in (i+1):(ncat-1))  for(k in (j):(ncat-1))  for(l in (k+1):(ncat)) 
      {z=z+1; m1a[[z]]=unitvec(nvar*ncat,(a-1)*ncat+i)+unitvec(nvar*ncat,(b-1)*ncat+j)+unitvec(nvar*ncat,(b-1)*ncat+k)+unitvec(nvar*ncat,(a-1)*ncat+l)}
   m1a = do.call(rbind,m1a);
   m1b = list(); z=0;
   for(a in 1:nvar) for(b in 1:nvar) for(i in 1:(ncat-2)) for(j in (i+1):(ncat-1))  for(k in (j):(ncat-1))  for(l in (k+1):(ncat)) 
      {z=z+1; m1b[[z]]=unitvec(nvar*ncat,(b-1)*ncat+i)+unitvec(nvar*ncat,(a-1)*ncat+j)+unitvec(nvar*ncat,(a-1)*ncat+k)+unitvec(nvar*ncat,(b-1)*ncat+l)}
   m1b = do.call(rbind,m1b);
   m1 = rbind(m1a,m1b);
   one = t(as.matrix(rep(1,(2*ncat-ncat^2-2*ncat^3+ncat^4)/24))); #length found with MMA
   m2 = one; for(i in 1:(nvar^2-1)) m2 = DirectSum(m2,one)
   m2 = DirectSum(m2,m2)
   m3 = diag(dim(m2)[[1]]/2)
   m3 = cbind(m3,-m3)
   
   matlist1 = list(m3,m2,m1);
   funlist1 = list("log","exp","log");

   list(append(matlist1,matlist0),append(funlist1,funlist0))

}
#spec.ordinal.disp.L(c(2,3))

spec.ordinal.disp.A <- function(z){ #z=c(nvar,ncat)
   spec0 = spec.ordinal.disp.L(z);
   spec0[[1]][[1]] = 4*spec0[[1]][[1]];  
   spec0[[2]][[1]] = "identity";  
   spec0
}
#spec.ordinal.disp.A(c(2,3))


MultiCoefficient <- function(coeff,k){
    if( k == 1 ) newcoeff <- coeff
    else{
        matlist <- coeff[[1]] 
        funlist <- coeff[[2]]
        q <- length(matlist)
        newmatlist <- list()
        for (i in 1:q) {
            newmatlist[[i]] <- matlist[[i]]
            for (j in 2:k) newmatlist[[i]] <- DirectSum(newmatlist[[i]],matlist[[i]])   } 
        newcoeff <- list( newmatlist, funlist) }
    return(newcoeff)
}


SpecifyCoefficient <- function(
      name,        # Name of the coefficient 
      arg = NULL,  # Coefficient specific arguments
      rep = 1,
      data = NULL  # Data set. If provided, matrices designed for MEL are provided.
      ){
    coeff <- switch( name,
        "Mean" = spec.mean(arg),
        "Variance" = spec.variance(arg),
        "StandardDeviation" = spec.variance(arg),
        "Entropy" = spec.entropy(arg),
        "DiversityIndex" = spec.diversityindex(arg),
        "GiniMeanDifference" = spec.ginimeandifference(arg),
        "Covariance" = spec.covariance(arg),
        "Correlation" = spec.correlation(arg),
        "CohenKappa" = spec.cohenkappa(arg), 
        "CronbachAlpha" = spec.cronbachalpha(arg, data),
        "Hij" = spec.hij(),
        "KendallTau" = spec.kendalltau(arg),
        "GoodmanKruskalGamma" = spec.gamma(arg),
        "DifferenceInProportions" = spec.diffprop(arg), #arg has the form (var,condvar,dim)
        "ConditionalProbabilities" = spec.condprob(arg),
        "Probabilities" = spec.prob(arg),
        "LogProbabilities" = spec.logprob(arg),
        "LoglinearParameters" = spec.loglinearpars(arg), #arg=dim
        "LogOddsRatio" = spec.logOR(arg),
        # NOT YET AVAILABLE "MokkenHij" = spec.mokkenhij(),
        # NOT YET AVAILABLE "AllMokkenHij" = spec.allmokkenhij(arg,data),
        # NOT YET AVAILABLE "MokkenHj" = spec.mokkenhj(arg),
        "AllMokkenHj" = spec.allmokkenhj(arg,data),
        "MokkenH" = spec.mokkenh(arg, data),
        "OrdinalLocation-A" = spec.ordinal.loc.A(arg),
        "OrdinalLocation-L" = spec.ordinal.loc.L(arg),
        "OrdinalDispersion-A" = spec.ordinal.disp.A(arg),
        "OrdinalDispersion-L" = spec.ordinal.disp.L(arg)
        )
    MultiCoefficient(coeff,rep)
}



#OLD
#JoinModels <- function(model1,model2){
#
#  model1 = tomodelform(model1);
#  model2 = tomodelform(model2);#
#
#  bt1 = model1[[1]][[1]]
#  matlist1 = model1[[1]][[2]][[1]]
#  funlist1 = model1[[1]][[2]][[2]]
#  at1 = model1[[1]][[3]]
#  x1 = model1[[2]]

#  bt2 = model2[[1]][[1]]
#  matlist2 = model2[[1]][[2]][[1]]
#  funlist2 = model2[[1]][[2]][[2]]
#  at2 = model2[[1]][[3]]
#  x2 = model2[[2]]#

#  bt = DirectSum(bt1,bt2)
#  matlist = list();
#  for(i in 1:length(matlist1)) matlist[[i]] = DirectSum(matlist1[[i]],matlist2[[i]])
#  funlist = if( length(funlist1)==length(funlist2) && funlist1[[1]]==funlist2[[1]] ) funlist1 else break;
#  at = rbind(at1,at2)
#  model = list(bt,list(matlist,funlist),at);
#  
#  return(model)
#
#}

#design matrix not yet implemented
JoinModels <- function(...){
   modellist = list(...);
   modellist = lapply(modellist,tomodelform);
   nmod = length(modellist);
   if(nmod==1) return(modellist);
   
   btlist = lapply(modellist,function(x){x[[1]]$bt});
   matlistlist = lapply(modellist,function(x){x[[1]]$coeff$matlist});
   funlistlist = lapply(modellist,function(x){x[[1]]$coeff$funlist});
   atlist = lapply(modellist,function(x){x[[1]]$at});
   
   bt = do.call("DirectSum",btlist);
   at = do.call("rbind",atlist);
   
   rr=list();
   for(i in 1:length(matlistlist[[1]])) rr[[i]] = lapply(matlistlist,function(x){x[[i]]})
   matlist = lapply(rr,function(x){do.call("DirectSum",x)})
   
   funlist = funlistlist[[1]];
   
   joinmod = list()
   joinmod$bt = bt
   joinmod$coeff$matlist = matlist
   joinmod$coeff$funlist = funlist
   joinmod$at = at
   
   return(joinmod)
}
#a=matrix(1:12,4,3);b=matrix(rep(2,6),3,2);m1=list(a,"log",b);JoinModels(m1,m1,m1)



###############################################################################################
###############################################################################################
###############################################################################################
# Misc functions (Andries)


#simplification using ftable
RecordsToFrequencies <- function(data){ c(t(ftable(data))) }

#dat=rbind(c(1,1,2),c(2,1,3))
#RecordsToFrequencies(dat,c(2,2),SelectColumns=c(1,2))
#RecordsToFrequencies(dat,c(2,2,3))


# 4. Required functions
Gfunction <- function( m, coeff ){
# requires functions (1) phi, (2) dphi,
  matlist <- coeff$matlist     # read design matrices from the model
  funlist <- coeff$funlist     # read algebraic operations from the model
  q <- length(funlist)
  g <- list()
  g[["g.0"]] <- m
  for (i in 1:q) g[[i+1]] <- phi(matlist[[q+1-i]], g[[i]], funlist[[q+1-i]])
  Gfun <- list()
  Gfun[["Gfun.0"]] <- diag(length(m))
  for (i in 1:q) Gfun[[i+1]] <- dphi(matlist[[q+1-i]], g[[i]], Gfun[[i]], funlist[[q+1-i]])
  return( Gfun[[q+1]] )
}

gfunction <- function(m,coeff){
# requires functions (1) phi
    matlist <- coeff$matlist     # read design matrices from the model
    funlist <- coeff$funlist     # read algebraic operations from the model
    q <- length(funlist)
    g <- list()
    g[["g.0"]] <- m
    for (i in 1:q) {g[[i+1]] <- phi(matlist[[q+1-i]], g[[i]], funlist[[q+1-i]])}
    return( g[[q+1]] )
}

phi <- function(A,f, action){
# Numerical values are translations h(A %*% f) = A %*% f -
    eps=1E-80;
    switch(action,
    "identity" = A %*% f,
    "exp"      = A %*% exp(f),
    "log"      = A %*% log(abs(f)+eps),
    "sqrt"     = A %*% sqrt(f),
    "xlogx"    = A %*% (-f*log(f+eps)),
    "xbarx"    = A %*% (f*(1-f))  # x(1-x)
    )
}

dphi <- function(A,f,df, action){
  eps=1E-80;
  switch(action,
  "identity" = A %*% df,
  "exp"      = A %*% (as.numeric(exp(f)) * df),
  "log"      = A %*% (as.numeric(1/(f+eps)) * df),
  "sqrt"     = A %*% (as.numeric(1/(2*sqrt(f))) * df),
  "xlogx"    = A %*% (as.numeric(-1-log(f+eps)) * df),
  "xbarx"    = A %*% (as.numeric(1-2*f) * df),  # x(1-x)
  )
}

pds <- function(n,m,lambda=1){
# Power divergence statistics (Cressie + Read)
   if(length(m)!=length(n))stop("m and n have unequal lengths")
   n[n < 1e-100] <- 1e-100
   m[m < 1e-100] <- 1e-100
   if(lambda==0) return(2*sum(n*log(n/m)))
   if(lambda==-1) return(2*sum(m*log(m/n)))
   return(2/(lambda*(lambda+1))*sum(n*(((n/m)^lambda) - 1)))
}



###############################################################################################
###############################################################################################
###############################################################################################
# Printing statistics


printalgorithmdetails = function(fit){
   estmethod=switch(fit$Method,
      "ML"="Maximum Likelihood",
      "MDI"="Minimum Discrimination Information",
      "GSK"="GSK");
   cat(fit$Title,"\n")
   cat("\n")
   cat("Estimation method     = ",estmethod,"\n")
   cat("Time taken            = ",fit$TimeTaken," seconds","\n")
   cat("Number of iterations  = ",fit$NumberOfIterations,"\n") 
   cat("Convergence criterion = ",fit$ConvergenceCriterion,"\n")
}

printbasicstatistics = function(stats,showeig,satmod){
   if(satmod){
       cat("Sample size = ",stats$SampleSize,"\n")
       cat("Eigenvalues sample covariance matrix","\n")
       cat("\t",stats$Eigenvalues,"\n")
      } else
   if (stats$Method=="ML") { 
          cat("Loglikelihood ratio (G^2) = ",stats$LogLikelihoodRatio,"\n");
          cat("Chi-square (X^2)          = ",stats$ChiSquare,"\n")} else if 
      (stats$Method=="MDI"){ 
          cat("Discrimination information = ",stats$DiscriminationInformation,"\n")
          cat("Loglikelihood ratio (G^2)  = ",stats$LogLikelihoodRatio,"\n")} else if
      (stats$Method=="GSK"){ 
          cat("Wald statistic             = ",stats$WaldStatistic,"\n")}
   if(!satmod){
      cat("Degrees of freedom         = ",stats$DegreesOfFreedom,"\n")
      cat("p-value     = ",stats$PValue,"\n")
      cat("Sample size = ",stats$SampleSize,"\n")
      cat("BIC         = ",stats$BIC,"\n")
      if(showeig) {cat("\n"); cat("Eigenvalues of inverse covariance matrix of Lagrange multipliers: ", "\n", stats$Eigenvalues, "\n")}
   }
}


printadvancedstatistics = function( stats, satmod, ShowCorrelations ){

    obsval   = stats$ObservedCoefficients 
    fitval   = stats$FittedCoefficients 
    se       = stats$CoefficientStandardErrors
    zscores  = stats$CoefficientZScores
    adjres   = stats$CoefficientAdjustedResiduals
    cormat   = stats$CoefficientCorrelationMatrix
    coeffdim = stats$CoefficientDimensions
    varlabels = stats$CoefficientTableVariableLabels
    catlabels = stats$CoefficientTableCategoryLabels

    coeffdim  = rev(coeffdim);
    catlabels = rev(catlabels);

    nvar=length(coeffdim)
    if(nvar==1) {  
       if(!satmod) {
          cat("   Observed values             = ",obsval,sep="\t","\n")
          cat("   Fitted values               = ",fitval,sep="\t","\n");
          cat("   Ese of fitted values        = ",se,sep="\t","\n");
          cat("   Fitted values / std errors  = ",zscores,sep="\t","\n");
          cat("   Adjusted residuals          = ",adjres,sep="\t","\n") }
       else{
          cat("   Observed values              = ",obsval,sep="\t","\n")
          cat("   Ese of observed values       = ",se,sep="\t","\n");
          cat("   Observed values / std errors = ",zscores,sep="\t","\n")} }
    else {    
       if(!satmod) {
          cat("   Observed values: ","\n");  print(array(obsval,coeffdim,dimnames=catlabels))
          cat("   Fitted values: ","\n");  print(array(fitval,coeffdim,dimnames=catlabels))
          cat("   Ese of fitted values: ","\n");  print(array(se,coeffdim,dimnames=catlabels))
          cat("   Fitted values / std errors: ","\n");  print(array(zscores,coeffdim,dimnames=catlabels))
          cat("   Adjusted residuals: ","\n");  print(array(adjres,coeffdim,dimnames=catlabels))  }
       else{
          cat("   Observed values: ","\n");  print(array(obsval,coeffdim,dimnames=catlabels))
          cat("   Ese of observed values: ","\n");  print(array(se,coeffdim,dimnames=catlabels))
          cat("   Observed values / std errors: ","\n");  print(array(zscores,coeffdim,dimnames=catlabels)) }
    }

    if(ShowCorrelations && length(se)>1) {
       cat("Correlation matrix:","\n")
       cat(cormat,"\n","\n") }
}


tocatlabels = function(lab,dim){
   nvar = length(dim);
   catlabels=vector("list",nvar)  #labels for categories
   for(i in 1:nvar){ 
       for(j in 1:dim[[i]]) {catlabels[[i]][[j]]=
          if(is.list(lab[[i]])) paste(lab[[i]][[1]],"=",lab[[i]][[2]][[j]],sep="")
          else paste(lab[[i]],"=",j,sep="")}}
   catlabels
   }

#tocatlabels(c(1,2),c(2,2))
#tocatlabels(list(list("G",c("men","women")),"T"),c(2,3))

printThetaAndBeta = function(stats, showcoeffs, showpars, ShowCorrelations, ParameterCoding, satmod, modeltype="Manifest"){

   obsval = stats$ObservedCoefficients 
   fitval = stats$FittedCoefficients 
   coeffdim = stats$CoefficientDimensions
   varlabels = stats$CoefficientTableVariableLabels
   catlabels = stats$CoefficientTableCategoryLabels
   nvar=length(coeffdim)

   if(showcoeffs){
      cat("Statistics for the coefficients: ", "\n")
      cat("Variables = ",varlabels," (dim = ",coeffdim,")","\n",sep="")
      printadvancedstatistics(stats, satmod, ShowCorrelations)
   }

   printBeta = function(beta){
      efflab = beta$CoefficientTableVariableLabels
      effdim = beta$CoefficientDimensions
      cat("Effect = ");cat(as.character(efflab,sep=","),sep=",");cat(" (dim = ");cat(effdim,sep=",");cat(")","\n",sep="")
      printadvancedstatistics(beta, satmod, ShowCorrelations)  
   }

   if(showpars){
      cat("\n")
      cat("Statistics for the parameters: ", "\n")
      subvar = allsubsets(c(1:nvar));
      for(i in 1:length(subvar)) printBeta(stats$Parameters[[i]])
   }

}




###############################################################################################
###############################################################################################
###############################################################################################
# Functions relating to model specification


#puts marginal model into standard form (bt,theta,at,d)
tomargmodform = function(margmodel){

   #case 0: margmodel=None
   if( isTRUE(margmodel == "None") ) return(margmodel)

   #find out if d is given in model, if so, drop it from margmodel
   if({lm =length(margmodel);is.vector(margmodel[[lm]])&&is.numeric(margmodel[[lm]])&&!is.matrix(margmodel)}){
         d = margmodel[[lm]];
         margmodel = margmodel[-lm]
         }
   else d = "None";

   mm = list()  #marginal model in standard form

   #case 1: single matrix is given
   if( data.class(margmodel)=="matrix" ){
      mm$at = margmodel;
      id = diag(nrow(mm$at));
      mm$bt = diag(nrow(mm$at));
      mm$coeff$matlist = list(id)
      mm$coeff$funlist = list("identity") }
   else if(length(margmodel)>4) {print("Error in model specification: length of model should be at most four, in the form list(bt,coeff,at,d) representing the equation bt.coeff(at.pi)=d");break}

#   if(length(margmodel)==4) {print("Error in model specification: fourth element should be a numeric constant vector d from the formula bt.coeff(at.pi)=d ");break}
   
   #case 2: form is {bt,Log,at}
   else if(length(margmodel) == 3 && data.class(margmodel[[2]])=="character"){
      mm$bt = margmodel[[1]];
      mm$at = margmodel[[3]];
      if(ncol(mm$bt)!=nrow(mm$at)) {cat("Incompatible matrix dimensions in model specification",dim(mm$bt),"and",dim(mm$at),"\n");break;};
      mm$coeff$matlist = list(diag(ncol(mm$bt)))
      mm$coeff$funlist = list(margmodel[[2]])
   }
   
   #case 3: form is standard: {bt,coeff,at}, with coeff={matlist,funlist}
   else if(length(margmodel) == 3 && !data.class(margmodel[[2]])=="character" ) {
      mm$bt = margmodel[[1]];
      mm$coeff$matlist = margmodel[[2]][[1]];
      mm$coeff$funlist = margmodel[[2]][[2]];
      mm$at = margmodel[[3]];
   }
   
   #hereafter: length[margmodel]=2 ################################
   
   #form is eg {Log,at}
   else if(data.class(margmodel[[1]])=="character" && data.class(margmodel[[2]])=="matrix"){
      id=diag(nrow(margmodel[[2]]));
      mm$bt = id;
      mm$coeff$matlist = list(id);
      mm$coeff$funlist = list(margmodel[[1]]);
      mm$at = margmodel[[2]];
   }

   #form is {bt,Log}
   else if(data.class(margmodel[[2]])=="character" && data.class(margmodel[[1]])=="matrix"){
      id=diag(ncol(margmodel[[1]]));
      mm$bt = margmodel[[1]];
      mm$coeff$matlist = list(id);
      mm$coeff$funlist = list(margmodel[[2]]);
      mm$at = id;
   }
    
   
   #form is {bt,coeff}
   else if( data.class(margmodel[[1]])=="matrix" && length(margmodel[[2]]) == 2 ) list( margmodel[[1]], margmodel[[2]], "None" )
   
   #form is {coeff,None}
   else if( isTRUE(margmodel[[2]] == "None") && length(margmodel[[1]]) == 2 ){
      mm$coeff$matlist = margmodel[[1]][[1]];
      mm$coeff$funlist = margmodel[[1]][[2]];
      df = nrow(mm$coeff$matlist[[1]]);
      mm$bt = diag(df);
      mm$at = "None"
   } 
   
   #form is {coeff,at}
   else if( data.class(margmodel[[2]])=="matrix" && length(margmodel[[1]]) == 2 ){
      mm$coeff$matlist = margmodel[[1]][[1]];
      mm$coeff$funlist = margmodel[[1]][[2]];
      mm$at = margmodel[[2]];
      df = nrow((mm$coeff$matlist[[1]]));
      mm$bt = diag(df);
   }
   
   #form is coeff with coeff={matlist,funlist}
   else if( length(margmodel) == 2 && data.class(margmodel[[1]][[1]])=="matrix" && data.class(margmodel[[2]])=="list" ){
      mm$coeff$matlist = margmodel[[1]];
      mm$coeff$funlist = margmodel[[2]];
      k=length(mm$coeff$matlist);k=length(mm$coeff$matlist[[k]][1,])
      mm$at = diag(k)
      df = nrow((mm$coeff$matlist[[1]]));
      mm$bt = diag(df);
   }
   
   else{ cat("Cannot recognize model specification", margmodel); break}

   #now homogenize margmodel specification
   k = ncol(mm$coeff$matlist[[length(mm$coeff$matlist)]]);
   prob <- spec.prob(k);
   mm$coeff$matlist <- c(mm$coeff$matlist, prob[[1]]);
   mm$coeff$funlist <- c(mm$coeff$funlist, prob[[2]]);

   #add a constant vector d
   d = if( isTRUE(d=="None") ) rep(0,nrow(mm$bt)) else d;
   mm$d=d

   #check if matrix dimensions are correct
   dim0 = list(c(0,length(d)));
   dimfirst = list(dim(mm$bt));
   dimlast  = list(dim(mm$at));
   dimlist  = lapply(mm$coeff$matlist,dim);
   dimlist  = append(dimlist,dimfirst,after=0);
   dimlist  = append(dimlist,dim0,after=0);
   dimlist  = if(length(dimlast[[1]])>0) append(dimlist,dimlast) else dimlist
   for(i in 2:length(dimlist)){if(dimlist[[i-1]][[2]]!=dimlist[[i]][[1]]){print("Incompatible matrix dimensions. Dimensions are: ");print(dimlist);stop()}}
   
   return(mm)
}

#mod=list(matrix(1:6,2,3),"log",matrix(1:12,3,4),c(1,2))
#tomargmodform2(mod)

#margmodel="None"
#tomargmodform2(margmodel)

#mod=rbind(c(1,2))
#tomargmodform2(mod)

#margmodel=list( rbind(c(1,2)), "log", rbind(c(1,2)) )
#tomodelform2(margmodel)



#puts model into standard form: (margmodel,x) where x denotes design matrix for loglinear model for full table
tomodelform = function(model){

   desmatQ = function(mat){ if( !(data.class(mat)=="matrix") ) FALSE else if(nrow(mat)>ncol(mat)) TRUE else FALSE }
  
   standardformmodel =
  
      if(desmatQ(model)) list("None",model)                                                    #form is model=x
  
      else if( isTRUE(model[[1]]=="None") && desmatQ(model[[2]]) ) model                       #form is model={None,x}
     
      else if( length(model)==2 && desmatQ(model[[2]]) && data.class(model[[1]])=="character" ){   #form is eg {Log,x}: only used for coefficients, no df for margmod
         x = model[[2]]
         id = diag(nrow(x))
         list( tomargmodform(list(id,list(list(id),list(model[[1]])), "None")), x )}
        
      else if( length(model) == 2 && (desmatQ(model[[2]]) || isTRUE(model[[2]] == "None") ) ){ #form is {margmodel,x} or {margmodel,None}
         margmod = tomargmodform(model[[1]])
         x = model[[2]]
         list(margmod, x)}
        
      else list(tomargmodform(model), "None")                                                  #form is "margmod", no design matrix
      
   return(standardformmodel)
}

#x=rbind(c(1,2))
#tomodelform(x)

modeldf = function(model,n){
   tblsize = tablesize(model);
   latdim = tblsize / length(n)
   margmod = model[[1]];
   x = model[[2]];
   df1 = if (isTRUE(margmod == "None")) 0 else nrow(margmod$bt)
   df2 = if (isTRUE(x == "None")) 0 else max(0, tblsize/latdim - ncol(x));
   c(df1,df2)
   }

tablesize = function(model) {
   if(length(model)==2){
      margmod = model[[1]];   x = model[[2]];
      if      (!isTRUE(x == "None")) nrow(x)
      else if (!isTRUE(margmod$at == "None")) ncol(margmod$at)
      else    {k=length(margmod$coeff$matlist); ncol(margmod$coeff$matlist[[k]])}  }
   else  #model is of the form margmod
      if(isTRUE(model$at=="None")){k=length(model$coeff$matlist); ncol(model$coeff$matlist[[k]])}
      else ncol(model$at)
}

#mod=list(matrix(1:2,2,2),"log",matrix(1:6,2,3))
#mod=list(matrix(1:2,2,2),"log","None")
#mod2=tomargmodform(mod)
#tablesize(mod2)



###############################################################################################
###############################################################################################
###############################################################################################
# Fitting procedures

margmodfit = function(n, model, MaxSteps=1000, MaxError=1e-25, StartingPoint="Automatic", ShowProgress=TRUE, MaxStepSize=1,Method="ML"){

    eps= 1.e-80;
    starttime = proc.time()[3]

    margmod = model[[1]]
    x     <- model[[2]]
    bt    <- margmod$bt
    coeff <- margmod$coeff
    at    <- margmod$at
    d     <- margmod$d
        
    step <- 1             # Initial step size (no modification yet)
    k <- 0                # index for iterations k = 1, 2, .... . k = 1 is initial iterration
    w <- list()
    H <- matrix()

    m  <- if(isTRUE(StartingPoint=="Automatic") && !is.matrix(x) ){ n + 1e-3 } 
          else if(is.matrix(x)) rep(sum(n)/length(n),length(n))
          else StartingPoint        # initial estimates for expected frequencies
    m  <- as.numeric(m)
    logm <- log(m)                      # logarithm of initial expected frequencies
    atm  <- if(data.class(at)!="matrix") m else at %*% m
    g  <- matrix(gfunction(atm,coeff))  # initial values of g in a C x 1 vector
    Gt <- if(data.class(at)!="matrix") Gfunction(atm,coeff) else Gfunction(atm,coeff) %*% at   # initial values of G in an L x C matrix
    h  <- bt %*% g - d
    Ht <- bt %*% Gt
    H  <- t(Ht)
    df <- nrow(bt)                      # number of constraints
    mu <- matrix(0,df,1)                # initial estimates of Lagrange multipliers in a C x 1 vector
    olderror <- 100000                  # initial values of error

    if(!isTRUE(!ShowProgress)) cat("Iteration      ","step size      ","G^2             ","Error","\n")
    olderror <- 1.e100;

    # Iteration = 1, ..., k
    repeat{
        k <- k + 1
        atm <- if(data.class(at)!="matrix") m else at %*% m
        g <- gfunction(atm,coeff)
        Gt <- if(data.class(at)!="matrix") Gfunction(atm,coeff) else Gfunction(atm,coeff) %*% at
        h <- bt %*% g - d
        Ht <- (bt %*% Gt )
        H <- t(Ht)
        HDH <- Ht %*% ( m * H )
        ndivm <- (m+eps)/(n+eps)
        dl = switch(Method,
            "ML"  = n - m,        #derivative of multinomial likelihood wrt log(m)
            "MDI" = -m*log(ndivm))#derivative of MDI criterion wrt log(m)
        dl2 = switch(Method,
            "ML"  = ndivm - 1,   #derivative of multinomial likelihood wrt m
            "MDI" = -log(ndivm)) #derivative of MDI criterion wrt m
        lambda <- - solve( HDH, Ht %*% dl + h )
        loglinincr = if(isTRUE(x=="None")) 0 else x %*% (solve(t(x) %*% (m * x)) %*% (t(x) %*% dl)) - dl2
        incr <- dl/m + H %*% lambda + loglinincr  #cannot replace dl/m by dl2 for some numerical reason
        newerror <-  (t(incr) %*% (m * incr))[1]
        errorratio = if(newerror==0.) 1 else sqrt(olderror/newerror);
        newstep = if(errorratio<MaxStepSize && errorratio>.01*MaxStepSize) errorratio/2 else MaxStepSize; 
        olderror = newerror;
        logm <- logm + newstep*incr
        m <- exp(logm)
        m <- as.numeric(m)
        m[m < 1e-80] <- 1e-80
        if( isTRUE(ShowProgress) || (data.class(ShowProgress)=="numeric" && isTRUE(k%%ShowProgress==0)) ){ cat(k,"             ",newstep,"             ",pds(n,m,0),"         ",newerror,"\n") }
        if ( k >= MaxSteps || abs(newerror) < MaxError ) break
    }

    fit=list();   fit$FittedFrequencies=m;  fit$ConvergenceCriterion=newerror;    fit$NumberOfIterations=k;    fit$TimeTaken <- proc.time()[3] - starttime
    return(fit) 
}


margmodfitEM <- function( n, latdim, model, maxoutersteps=1000, StartingPoint="Automatic", MaxError=1e-20, MaxInnerSteps=2, ShowProgress=TRUE, MaxStepSize=1){

    starttime = proc.time()[3]
    eps <- 1e-80
    k <- 0
    collapse = function(q){ apply( matrix( q, nrow = latdim ), 2, sum ) }
    expand = function(q){ rep(q, each=latdim ) }
    v <- c(1:(length(n)*latdim))
    mhat = expand(n/latdim)
    mhat = mhat + sum(n)/100/length(mhat)  #startingpoint

    cat("Iteration      ","G^2             ","Error","\n")

    repeat{ k <- k + 1
        nhat <- mhat * expand( (n+eps) / (collapse(mhat)+eps) ) # E-step
        newmhat <- margmodfit( nhat, model, MaxSteps=MaxInnerSteps, StartingPoint=mhat, ShowProgress=FALSE, MaxStepSize=MaxStepSize )$FittedFrequencies  # M-step
        error <- sum( abs( mhat - newmhat ) )
        mhat <- newmhat
        if( isTRUE(ShowProgress) || (data.class(ShowProgress)=="numeric" && isTRUE(k%%ShowProgress==0)) ){ 
             cat("  ",k,"             ",pds(n,collapse(newmhat),0),"         ",error,"\n") }
        if ( k >= maxoutersteps || abs(error) < MaxError ) break
    }

    fit=list();   fit$FittedFrequencies=mhat;  fit$ConvergenceCriterion=error;    fit$NumberOfIterations=k;    fit$TimeTaken <- proc.time()[3] - starttime
    return(fit) 
}


gsk <- function( n, model, CoefficientDimensions="Automatic", Labels="Automatic" ){

    model = tomodelform(model)
    margmod = model[[1]]
    x     <- model[[2]]
    bt    <- margmod$bt
    coeff <- margmod$coeff
    at    <- margmod$at
    d     <- margmod$d

    atn <- if(data.class(at)!="matrix") n else at %*% n

    h  <- bt %*% gfunction(atn,coeff) - d
    Gt <- if(data.class(at)!="matrix") Gfunction(atn,coeff) else Gfunction(atn,coeff) %*% at
    Ht <- (bt %*% Gt )
    H <- t(Ht)

    GDG <- Gt %*% (n * t(Gt)) 
    GDH <- Gt %*% (n * H)
    HDH <- Ht %*% (n * H)
    
    covresid <- GDH %*% solve(HDH) %*% t(GDH) 
    covtheta  <- GDG - covresid
    obsval <- gfunction(atn,coeff)
    fitval <- covtheta %*% solve( GDG, obsval )
    adjres = (obsval - fitval)/sqrt(diag(covtheta))

    stats = list()
    stats$Method = "GSK"
    stats$LogLikelihoodRatio = 0
    stats$DegreesOfFreedom <- length(h)
    stats$WaldStatistic <- h %*% solve( HDH, h )
    stats$SampleSize <- sum(n)
    stats$BIC  <- stats$WaldStatistic - stats$DegreesOfFreedom * log( sum(n) )
    stats$PValue <- signif(1-pchisq(stats$WaldStatistic,stats$DegreesOfFreedom),5)
    stats$Eigenvalues = eigen(HDH, only.values = TRUE)$values
    stats = c(stats,coefficientstats(obsval,fitval,covtheta,covresid,CoefficientDimensions,Labels,FALSE))

    return(stats)
}
#  model <- list( bt, coeff2, at, c(0.1) )
#  m <- gsk(n,model)



#parameters
getBetas = function(efflab,effcatlab,effdim,obsval1,fitval1,covtheta1,covresid1,coeffdim,varlabels,satmod,noresid,ParameterCoding="Effect"){
       effmat = betamat(varlabels,efflab,coeffdim,ParameterCoding);
       obsval   = effmat %*% obsval1;
       fitval   = effmat %*% fitval1
       covtheta = effmat %*% covtheta1 %*% t(effmat)
       covresid = if(satmod) 0 else effmat %*% covresid1 %*% t(effmat)

       var <- diag(covtheta)
       fitval[abs(fitval) < 1e-20] <- 0
       var[var <= 0] <- 1e-80
       se <- sqrt(var)
       zscores <- fitval / se
       if(noresid==1||isTRUE(covresid==list(0))) adjres=NULL else {
          vr <- diag(covresid);
          vr[vr <= 0] <- 1e-80
          adjres <- if(is.null(obsval)||is.null(fitval)||isTRUE(vr==0)) 0 else (obsval - fitval)/sqrt(vr)
          adjres[abs(adjres) < 1e-39] <- 0  }
       zscores[abs(zscores) < 1e-39] <- 0
       cormat = if(length(se)==1) 1 else diag(1/se) %*% covtheta %*% diag(1/se)

       beta = list()
       beta$ObservedCoefficients = obsval
       beta$FittedCoefficients = fitval
       beta$CoefficientStandardErrors = se
       beta$CoefficientZScores = zscores
       beta$CoefficientAdjustedResiduals = adjres
       beta$CoefficientCovarianceMatrix = covtheta
       beta$CoefficientCorrelationMatrix = cormat
       beta$CoefficientAdjustedResidualsCovarianceMatrix = covresid
       beta$CoefficientDimensions = effdim
       beta$CoefficientTableVariableLabels = efflab
       beta$CoefficientTableCategoryLabels = effcatlab

       return(beta)
}


coefficientstats = function(obsval, fitval, covtheta, covresid, coeffdim, varlabels,satmod,eig="Automatic",ParameterCoding="Effect"){
    noresid = prod(obsval==fitval)

    var <- diag(covtheta)
    fitval[abs(fitval) < 1e-20] <- 0
    var[var <= 0] <- 1e-80
    se <- sqrt(var)
    zscores <- fitval / se
    if(noresid==1) adjres=NULL else {
       vr <- diag(covresid);
       vr[vr <= 0] <- 1e-80
       adjres <- if(is.null(obsval)||is.null(fitval)||isTRUE(vr==0)) 0 else (obsval - fitval)/sqrt(vr)
       adjres[abs(adjres) < 1e-39] <- 0  }
    zscores[abs(zscores) < 1e-39] <- 0
    cormat = if(length(se)==1) 1 else diag(1/se) %*% covtheta %*% diag(1/se)

    coeffdim = if(isTRUE(coeffdim=="Automatic")) c(length(obsval)) else coeffdim;
    nvar=length(coeffdim)
    varlabels = if(isTRUE(varlabels=="Automatic")) apply(rbind(c(1:nvar)),1,function(x){paste("Var",x,sep="")}) else varlabels;
    catlabels=tocatlabels(varlabels,coeffdim);
    if(prod(coeffdim)!=length(obsval)){print("Error: CoefficientDimensions incorrect, the product of the dimensions should equal the number of coefficients"); stop()}
    if(length(coeffdim)!=length(varlabels)){print("Error: CoefficientDimensions has different length than Labels"); stop()}

    varlabels1 = list();   #varlabels1 contains only variable labels, not category labels
    for(i in 1:nvar){ varlabels1[[i]] = if(is.list(varlabels[[i]])) varlabels[[i]][[1]] else varlabels[[i]] }
    varlabels1 = as.character(varlabels1)

    nvar = length(varlabels)
    subvar = allsubsets(c(1:nvar));
    beta = list()
    for(i in 1:length(subvar)) beta[[i]] = getBetas(varlabels1[subvar[[i]]],catlabels[subvar[[i]]],coeffdim[subvar[[i]]],obsval,fitval,covtheta,covresid,coeffdim,varlabels1,satmod,noresid,ParameterCoding=ParameterCoding)

    stats = list()
    stats$ObservedCoefficients = obsval
    stats$FittedCoefficients = fitval
    stats$CoefficientStandardErrors = se
    stats$CoefficientZScores = zscores
    stats$CoefficientAdjustedResiduals = adjres
    stats$CoefficientCovarianceMatrix = covtheta
    stats$CoefficientCorrelationMatrix = cormat
    stats$CoefficientAdjustedResidualsCovarianceMatrix = covresid
    stats$CoefficientDimensions = coeffdim
    stats$CoefficientTableVariableLabels = varlabels1
    stats$CoefficientTableCategoryLabels = catlabels
    stats$Parameters = beta

    class(stats) = "CMM"
    return(stats)
}


"summary.CMM" = function(object,...,ShowCoefficients=TRUE,ShowParameters=FALSE,ShowCorrelations=FALSE,ParameterCoding="Effect",Title=""){
    if(!is.null(object$Title)) print(object$Title)
    satmod = if(object$LogLikelihoodRatio==0) TRUE else FALSE
    printbasicstatistics(object,showeig=TRUE,satmod)
    cat("\n")
    printThetaAndBeta( object, ShowCoefficients, ShowParameters, ShowCorrelations, ParameterCoding, satmod )
}


getsamplestats = function(dat,coeff,CoefficientDimensions,ParameterCoding,Labels){
    eps <- 1e-80
    n   = if(data.class(dat)=="numeric") dat else c(t(ftable(dat)));
    coeff = tomargmodform(coeff); #yields list(bt,coeff,at,d)

    atn <- if(data.class(coeff$at)!="matrix") n else coeff$at %*% n
    Gt <- if(data.class(coeff$at)!="matrix") Gfunction( atn, coeff$coeff ) else  Gfunction( atn, coeff$coeff ) %*% coeff$at 
    GDG <- Gt %*% (n * t(Gt))

    stats = list()
    stats$Eigenvalues = eigen(GDG, only.values = TRUE)$values;
    stats$SampleSize <- sum(n)
    stats$LogLikelihoodRatio = 0

    covresid <- NULL
    covtheta <- GDG
    obsval = gfunction( atn + eps, coeff$coeff )
#   cat("\n")
    stats = c(stats,coefficientstats( obsval, obsval, covtheta, covresid, CoefficientDimensions, Labels, satmod=TRUE, ParameterCoding=ParameterCoding ))
}

SampleStatistics = function(dat, 
                            coeff, 
                            CoefficientDimensions="Automatic", 
                            Labels="Automatic", 
                            ShowCoefficients=TRUE, 
                            ParameterCoding="Effect", 
                            ShowParameters=FALSE, 
                            ShowCorrelations=FALSE,  
                            Title="", 
                            ShowSummary = TRUE
                   ){
    stats = getsamplestats(dat,coeff,CoefficientDimensions,ParameterCoding,Labels)
    class(stats) = "CMM"
    if(ShowSummary) summary(stats, ShowParameters=ShowParameters, ShowCorrelations=ShowCorrelations, ParameterCoding=ParameterCoding)
    return(stats)
}


getbasicstatistics = function(n,manmhat,df,eig){
   stats = list()
   stats$SampleSize = sum(n);
   stats$DegreesOfFreedom = df[1]+df[2];
   stats$LogLikelihoodRatio = pds(n,manmhat,0)
   stats$ChiSquare = pds(n,manmhat,1)
   stats$DiscriminationInformation = pds(n,manmhat,-1)
   stats$BIC = stats$LogLikelihoodRatio - stats$DegreesOfFreedom * log(sum(n))
   stats$PValue = signif(1-pchisq(stats$LogLikelihoodRatio,stats$DegreesOfFreedom),5)
   stats$Eigenvalues = eig
   stats$ManifestFittedFrequencies = manmhat
   return(stats)
}

getmodelstats = function(dat, mhat, model, coeff, CoefficientDimensions, Labels, Method, satmod=FALSE, ParameterCoding="Effect"){
    eps = 1e-80
    stats = list()
    stats$Method = Method
    mm=coeff

    n   = if(data.class(dat)=="numeric") dat else c(t(ftable(dat)));

    margmod = model[[1]]

    latdim = length(mhat) / length(n);
    collapse = function(q){ apply( matrix( q, nrow = latdim ), 2, sum ) }
    expand = function(q){ rep(q, each=latdim ) }
    nhat <- if(latdim==1) n else mhat * expand( (n+eps) / (collapse(mhat)+eps) ) 
    manmhat = if(latdim == 1) mhat else collapse(mhat);

    Gfun = if(data.class(mm$at)!="matrix") Gfunction( mhat, mm$coeff ) else Gfunction( mm$at %*% mhat, mm$coeff ) %*% mm$at
    Gt   = if(is.null(mm$bt)) Gfun else mm$bt %*% Gfun
    Ht <- if(data.class(margmod$at)!="matrix")  margmod$bt  %*% Gfunction( mhat, margmod$coeff )  else margmod$bt  %*% Gfunction( margmod$at  %*% mhat, margmod$coeff )  %*% margmod$at
    obsval <- if(data.class(mm$at)!="matrix") gfunction( nhat, mm$coeff) else gfunction( mm$at %*% nhat, mm$coeff)
    fitval <- if(data.class(mm$at)!="matrix") gfunction( mhat, mm$coeff) else gfunction( mm$at %*% mhat, mm$coeff)

    GDG <- Gt %*% (mhat * t(Gt)) 
    GDH <- Gt %*% (mhat * t(Ht))
    HDH <- Ht %*% (mhat * t(Ht))
    covresid <- GDH %*% solve(HDH) %*% t(GDH) 
    covtheta <- GDG - covresid
    eig = eigen(HDH, only.values = TRUE)$values;
    df = modeldf(model,n)

    stats = c(stats,getbasicstatistics(n,manmhat,df,eig))
    stats = c(stats,coefficientstats( obsval, fitval, covtheta, covresid, CoefficientDimensions, Labels, satmod, ParameterCoding=ParameterCoding ))
    return(stats)
}

ModelStatistics <- function(dat, fitfreq, model, coeff, CoefficientDimensions="Automatic",
    Labels="Automatic",Method="ML",ShowCoefficients=TRUE,ShowParameters=FALSE, ParameterCoding="Effect", ShowCorrelations=FALSE, Title="" ){
    
    stats = if(isTRUE(model=="SaturatedModel")) 
                    getsamplestats(dat,coeff,CoefficientDimensions,ParameterCoding,Labels)
            else {
                    model = tomodelform(model); 
                    coeff = tomargmodform(coeff); #yields list(bt,coeff,at,d)
                    getmodelstats(dat,fitfreq,model,coeff,CoefficientDimensions,Labels,Method)}
    class(stats) = "CMM"
    summary(stats,ShowCoefficients=ShowCoefficients,ShowParameters=ShowParameters, ParameterCoding=ParameterCoding, ShowCorrelations=ShowCorrelations)
    return(stats)
}


MarginalModelFit = function(dat, model, 
    ShowSummary=TRUE,
    MaxSteps=1000, MaxStepSize=1, MaxError=1e-20, StartingPoint="Automatic", MaxInnerSteps=2, ShowProgress=TRUE, CoefficientDimensions="Automatic",
    Labels="Automatic",ShowCoefficients=TRUE,ShowParameters=FALSE, ParameterCoding="Effect", ShowCorrelations=FALSE, Method="ML", Title="Summary of model fit"){

    n  = if(data.class(dat)=="numeric") dat else c(t(ftable(dat)));
    model = tomodelform(model); #put model in standard form: "list(margmod,x)"
    latdim = tablesize(model) / length(n)
    if (latdim != floor(latdim) || latdim < 1) cat("Error: incorrect dimensions of the matrix specification of the model","\n") 

    fit <- switch( Method, 
        "ML" = if (latdim == 1) 
                    margmodfit( n, model, MaxSteps=MaxSteps, MaxStepSize=MaxStepSize, StartingPoint=StartingPoint, MaxError=MaxError,ShowProgress=ShowProgress,Method="ML")
               else margmodfitEM( n, latdim, model, maxoutersteps=MaxSteps, MaxInnerSteps=MaxInnerSteps, MaxStepSize=MaxStepSize, StartingPoint=StartingPoint, 
                                    MaxError=MaxError, ShowProgress=ShowProgress),
        "MDI" = if (latdim != 1)
                    cat("Error: cannot do latent variable models for Method='MDI'","\n") 
            else    margmodfit( n, model, MaxSteps=MaxSteps, MaxStepSize=MaxStepSize, StartingPoint=StartingPoint, MaxError=MaxError,ShowProgress=ShowProgress,Method="MDI"),
        "GSK" = gsk( n, model) 
    )

    if( Method!="GSK"){
        coeff = model[[1]]
        coeff$bt=NULL
        fit = c(fit, getmodelstats(dat,fit$FittedFrequencies,model,coeff,CoefficientDimensions,Labels,Method))
    }

    fit$Method = Method 
    fit$Title = Title
    class(fit) = "CMM"

    if( ShowSummary){
        if(Method!="GSK") printalgorithmdetails(fit)
        summary(fit, ShowCoefficients=ShowCoefficients,ShowParameters=ShowParameters, ShowCorrelations=ShowCorrelations, ParameterCoding=ParameterCoding)
        cat("\n")
    }

    return(fit)
}
