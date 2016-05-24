PMDA.coxph <-
function(formula , data , k , m ){
  #function parameters
  prep<-preproc.coxph( data , m )
  data_int <- prep$data2
  data_fix <-prep$data1
  or<-prep$or
  I<-prep$I 
  
  
  mm<-model.matrix(formula , data)
  nc<-sapply(sapply(strsplit( as.character(formula)[2] , '\\+' ) , function(x) gsub(" " , "", x) ),nchar)
  nc2<-sapply(colnames(mm)[-1],nchar)
  sub1<-substr( colnames(mm)[-1] ,  nc+1 , nc2 )
  sub2<-paste(names(nc),sub1,sep=': ')
  colnames(mm) <- c( colnames( mm )[ 1 ] , sub2 )
  dim_beta<-ncol(mm)-1
  beta<-matrix(0,ncol=dim_beta,nrow=1)
  dn<-dimnames(mm)[[2]][-1]
  
  #Step1
  s0 <- MI.surv_1( m = m , data =  data , conf.int = F )$est
  s0$diff<-c(0,diff(1-s0$surv))
  
  #initial linear predictors
  Z<-mm[,-1]%*%t(beta)
    
  #Step 2  
  cat('\nIterates\n' )
  beta_iter<-matrix(NA,ncol=k,nrow=dim_beta,dimnames=list(dn,1:k))
  sigmac_iter<-matrix(NA,ncol=k,nrow=dim_beta,dimnames=list(dn,1:k))
    
  #progression bar
  i<-0
  pb <- txtProgressBar(style = 2 , char = '.' )
  
 
  repeat {
    i<-i+1 
    setTxtProgressBar(pb ,  i%%4/150 ) 
    if( i > k ){setTxtProgressBar(pb ,  0.02 )
      break }
  
  #print(i)    
  ss1<-apply(data_int , 1 , function(x ) subset( s0 , time >=  as.numeric(x['left']) & time <= as.numeric(x['right']) ) )
  tk2<-sapply(seq_len(nrow(data_int)) ,function(X)  ss1[[X]]$time)
  
  Z <- matrix( rep( Z[,1] , m ) , ncol = m , byrow = F )
  #Get the samples from drown betas from the normal mixture
  samples <- sapply( seq_len( nrow( data_int ) ) , function( X ) {
  #for(X in 1:nrow( data_int )){
  #print(X)
  pk2 <- sapply( 1:m , function( x ) ss1[[ X ]]$diff^exp( Z[ I ,  ][ X , x ] ) ) 
  pk2 <- matrix(pk2,ncol=m) 
  apply( pk2 , 2  , function( x ){
  if( sum(unlist(x)) & length(unlist(x)) > 1  ) sample(  tk2[[ X ]] , size = 1 , prob = ( x ) ) 
  else  data_int[ X , 'right' ]  }  ) } ) 
  
  samples <- matrix( unlist( samples ) , ncol = m , byrow = T )  
  samples2<-rbind(samples,data_fix)[or,]
  times<-as.vector(samples2)
    
  #surv<-Surv( time = times , event = rep( data$right != Inf , m )  , type = 'right')
  #surv2<-Surv( time = times , event = rep( data$right != Inf , m )  , type = 'mstate')
  #surv2[,2] <- surv[,2]
  #fitCI<-survfit( surv2 ~ 1 , weights = rep( 1 , length( times ) ) / m , conf.type = 'none' )  
  #pr <- fitCI$prev
  #t0<- fitCI$time
  #ne<- fitCI$n.event
  #nt<- rep(t0,ne)
  #nt  
  est_1 <- apply( samples2  ,  2  ,  get.est  ,  data  ,   formula )
      
  #get the betas in a matrix (x * k)
  betas <- matrix( unlist( sapply( est_1  , function( x ) x$beta ) )  ,  nrow  =  dim_beta  ,  dimnames  =  list( dn  ,  1:m ) )
  #get the mean of the betas over augmented datasets
  beta <- rowMeans( matrix( unlist( sapply( est_1  ,  function( x ) x$beta ) )  ,  nrow  =  dim_beta  ,  dimnames  =  list( dn  ,  1:m) ) )
      
  #Get the sigma in an d3 array
  sigma <- array( sapply( est_1 , function( x ) x$sigma )  ,  dim  =  c( dim_beta  ,  dim_beta  ,  m ) )
    
  #within variance: W
  W <- apply( sigma , c( 1 , 2 ) , mean )
    
  #update the CIF
  #compute exp (sum of B'Z )
  keep <- as.vector(sapply(strsplit( as.character(formula)[2] , '\\+' ) , function(x) gsub(' ' , '', x) ))
  adf<-as.data.frame(sapply(keep,function(x) as.data.frame((rep(data[,x],m)))))
  colnames(adf)<-keep
  r2<-as.numeric(rep( data$right != Inf , m ))
  adf2<-data.frame(times,status=r2,adf)
  s0<-BBS(formula = formula , time = 'times' , status = 'status' , data = adf2 , beta = beta) 
  s0 <- rbind( c( time = 0 , surv = 1 ) , s0 , c(max(times),tail(s0$surv,1)))
  s0$surv[is.na(s0$surv)] <- 0
  s0$diff<-c(0,diff(1-s0$surv))
  
  #Betwin variance: B with inflation factor 1/k 
  B <- ( 1 + ( 1 / m ) ) * ( ( betas - beta ) %*% t( betas - beta ) / ( m - 1 ) )
  #update de variance matrix
  sigmac <- W + B
  
  #update linear predictor
  Z <- mm[ , -1 ]%*%as.matrix( beta )
  beta_iter[ , i ] <- beta
  sigmac_iter[ , i ] <- diag( sigmac )
  }

  close( pb )
  ret<-list( beta = beta_iter , sigmac = sigmac_iter , vcov = sigmac , s0 = s0 )
  return( ret )
  }