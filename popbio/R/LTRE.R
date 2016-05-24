LTRE <-
function(trts, ref)
{
  ## treatments can be a single matrix or a list of matrices
  if(is.list(trts))
  {
    n<-length(trts)
    cont<-vector("list", n)
    if(!is.null(names(trts))){names(cont)<-names(trts)}
  }
  else if(is.matrix(trts)){ n<-1 ; trts<-list(trts) }
  else{stop("Treatments must be a matrix or list of matrices")}
        
  ## loop through each treatment 
  for (i in 1:n)
  {
     A1<-trts[[i]]
     ## matrix of differences (p 262)
     Dm<-A1-ref
     ## matrix halfway between treatment and reference
     Ac<-(A1+ref)/2
     ## sensitivity of Ac
     SAc<-sensitivity(Ac)
     # matrix of contributions
     Cm<-Dm*SAc
     ## output matrix or list of matrices depending on number of trts     
     if(n==1){cont<-Cm}
     else{ cont[[i]]<-Cm}
   }
  cont
}

