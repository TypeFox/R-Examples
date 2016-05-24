# sparse matrix version of the crossprod function
crossprodspam<-function(X){
  t(X) %*% X
}

# create an empty sparse matrix of class spam
make_sparse<-function(ncol, nrow = ncol) spam(0, nrow = nrow, ncol = ncol)

# create a square matrix of 1s of class spam (useful for block kronecker operations)
make_one<-function(ncol, nrow = ncol) spam(1, nrow = nrow, ncol = ncol)

# coerce to spam with given tolerance level  
make_spam<-function(M){
  if(class(M) == "matrix") as.spam(M)
  if(class(M) == "spam") as.spam(M)
  else as.spam(as.spam.dgCMatrix(as(M, "dgCMatrix")))
}

# take list of lists and collapse to a single longer list
make_flat<-function(List){
  new.list<-vector("list")
  for(i in 1:length(List)){
    if(!class(List[[i]]) == "list") new.list[[length(new.list)+1]]<-List[[i]]
    else
      for(j in 1:length(List[[i]])){
        new.list[[length(new.list)+1]]<-List[[i]][[j]]
      }
  }
  new.list
}