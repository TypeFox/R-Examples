
ANDA.crreg <-
function(formula , data , k = 25 , m = 10 , status , trans , cens.code ){
  
  #function parameters
  prep<-preproc.crreg( data , m = m , trans = trans , status = status , cens.code = cens.code )
  data_int<-prep$data2
  data_fix<-prep$data1
  or<-prep$or
  I<-prep$I
  
  r2 <- as.character(rep( data[ , status ] , m ) )
  r2<-factor(r2,levels=unique(c(cens.code,unique(r2))))
  
  mm<-model.matrix(formula , data)
  nc<-sapply(sapply(strsplit( as.character(formula)[2] , '\\+' ) , function(x) gsub(" " , "", x) ),nchar)
  nc2<-sapply(colnames(mm)[-1],nchar)
  sub1<-substr( colnames(mm)[-1] ,  nc+1 , nc2 )
  sub2<-paste(names(nc),sub1,sep=': ')
  colnames(mm) <- c( colnames( mm )[ 1 ] , sub2 )
  dim_beta<-ncol(mm)-1
  beta<-matrix(0,ncol=dim_beta,nrow=1)
  dn<-dimnames(mm)[[2]][-1]
  sigma<-matrix(rep(1e-3,dim_beta^2),ncol=dim_beta,dimnames=list(dn,dn))  
  beta_AN<-mvrnorm(n=m,mu=beta,Sigma=sigma)
  
  #Step1
  # generate sets for the algorithm initialization
  ci0 <- MI.ci_1( m = m , status = status , trans = trans , data = data , cens.code = cens.code , conf.int = FALSE )$est
  ci0$diff <- c( 0 , diff( ci0$est ) )
  
  #initial linear predictors
  Z<-apply( beta_AN , 1 , function(x)  mm[,-1]%*%as.matrix(x))
  
  #Step 2  
  # make matrix to get results
  cat('\nIterates\n' )
  beta_iter<-matrix(NA,ncol=k,nrow=dim_beta,dimnames=list(dn,1:k))
  sigmac_iter<-matrix(NA,ncol=k,nrow=dim_beta,dimnames=list(dn,1:k))
  W_iter<-matrix(NA,ncol=k,nrow=dim_beta,dimnames=list(dn,1:k))
  th_iter<-matrix(NA,ncol=k,nrow=dim_beta,dimnames=list(dn,1:k))
  
  #progression bar
  pb <- txtProgressBar(style = 2 , char = '.')
  i<-0
  repeat {
  i<-i+1 
  setTxtProgressBar(pb ,  i%%4/150 ) 
   
  if( i > k ){setTxtProgressBar(pb ,  0.02 )
  break }
  
  ss1<-apply( data_int , 1 , function(x ) subset( ci0 , time >=  as.numeric(x['left']) & time <= as.numeric(x['right']) ) )
  tk2<-lapply(seq_len(nrow(data_int)) ,function(X)  ss1[[X]]$time)
  
  samples <- sapply( seq_len( nrow( data_int ) ) , function( X ) {
  pk2 <- sapply( 1:m , function( x ) ss1[[ X ]]$diff^exp( Z[ I ,  ][ X , x ] ) ) 
  pk2 <- matrix(pk2,ncol=m) 
  apply( pk2 , 2  , function( x ){
  if( sum(x) & length(x) > 1  ) sample(  tk2[[ X ]] , size = 1 , prob = ( x ) ) 
  else  data_int[ X , 'right' ]  }  ) } )    
  
  samples <- matrix( unlist( samples ) , ncol = m , byrow = T )  
    
  if(dim(data_fix)[2]){
    samples2 <- rbind( samples , data_fix )[ or , ]
  }else{ samples2 <- samples }  
    
    
  times<-as.vector(samples2)
  ci<-survival::Surv( time = times , event = r2 , type = 'mstate')
  fitCI<-survival::survfit( ci ~ 1 , weights = rep( 1 , length( times ) ) / m , conf.type = 'none')  
  
  w <- which( fitCI$states == trans )
  pr <- fitCI$prev[ , w ]
  t0 <- fitCI$time
  
    
  #lines(t0,pr,type='s')
  #Get estimates
  est_1<-apply(samples2 , 2 , get.est.cr.cox , data = data , status = status , trans = trans , cens = cens.code ,  formula = formula )
    
  #Get the beta from the successive sets
  betas<-matrix(unlist(sapply( est_1  ,  function(x) x$beta))  ,  nrow  =  dim_beta  ,  dimnames  =  list(dn  ,  1:m))
  
  #Get mean of beta over samples
  beta <- rowMeans(matrix(unlist(sapply(est_1 , function(x) x$beta))  , nrow  =  dim_beta , dimnames  = list(dn  ,  1:m)))
  
  #Get the sigma square in an d3 array
  sigma<-array( sapply( est_1 , function(x) x$sigma)  ,  dim  =  c(dim_beta  ,  dim_beta  ,  m) )
  W<-apply(sigma,c(1,2),mean)
  
  #update de variance covariance matrix
  B<- ( 1 + ( 1/m ) ) * ( (betas-beta) %*% t(betas-beta) / (m-1) )
  
  #update de variance matrix
  sigmac <- W + B
  
  #sample beta_j from a possibliy multivariate normal with mean means of the betas and sigma obtained before
  beta_AN<-mvrnorm( n = m , mu = beta , Sigma = sigmac )
  
  #get new linear predictor
  Z<-apply( beta_AN , 1 , function(x)  mm[,-1]%*%as.matrix(x))
  
  ##update the CIF0
  #Estimate the survival function of the censoring distribution g_hat
  keep <- as.vector(sapply(strsplit( as.character(formula)[2] , '\\+' ) , function(x) gsub(' ' , '', x) ))
  adf<-as.data.frame(sapply(keep,function(x) as.data.frame((rep(data[,x],m)))))
  colnames(adf)<-keep
  adf2<-data.frame(times,status=r2,adf)
  ci0<-BBCI(formula = formula , time = 'times' , status = status , trans = trans , cens.code = cens.code , data =  adf2 ,beta =  beta  )
  ci0 <- rbind( c(time = 0, est = 0 ) , ci0 , c( time = max(times) , tail(ci0$est,1)) )
  ci0$diff <- c( 0 , diff( ci0$est ) )

  
  #save estimate of beta and sigma for the prensent iteration
  beta_iter[,i]<-beta
  sigmac_iter[,i]<-diag(sigmac)
    
  }
  close(pb)
  ret<-list( beta = beta_iter , sigmac = sigmac_iter , vcov = sigmac , ci0 = ci0 )
  return(ret)
}

