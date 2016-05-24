CVLwls1D <- function(y, t, kernel, npoly, nder, dataType, kFolds = 5 ){

  # If 'y' and 't' are vectors "cheat" and break them in a list of 10 elements
  if ( is.vector(y) && is.vector(t) && !is.list(t) && !is.list(y) ){
    if (length(t) < 21) {
      stop("You are trying to use a local linear weight smoother in a vector with less than 21 values.\n")
    }
    myPartition =   c(1:10, sample(10, length(t)-10, replace=TRUE));
    y = split(y, myPartition)
    t = split(t, myPartition)
    dataType = 'Sparse';
  } 
  
  # Make everything into vectors
  ncohort = length(t);
  tt  = unlist(t);
  yy  = unlist(y);
  ind = unlist(lapply( 1:ncohort, function(j) rep(j, times=length(t[[j]]))));
  yyn = yy[order(tt)];
  ind = ind[order(tt)];
  ttn = sort(tt);
  
  # Get minimum reasonable bandwidth
  a0=ttn[1];
  b0=ttn[length(ttn)];
  rang = b0-a0;
  dstar = Minb(tt, npoly+2);
  if (dataType != 'Dense'){
   h0 = 2.5*dstar;
  } else {
   h0 = dstar;
  }
  if (h0 > rang/4){
    h0 = h0*.75;
    warning(sprintf("Warning: the min bandwith choice is too big, reduce to %f !", (h0)  ))
  }    
  
  # Get the candidate bandwidths
  nbw = 11;
  bw = rep(0,nbw-1);
  n = length(unique(tt));
  for (i in 1:(nbw-1)){
    bw[i]=2.5*rang/n*(n/2.5)^((i-1)/(nbw-1)); # Straight from MATLAB
  }
  bw = bw-min(bw)+h0;
  
  #ave = rep(0, length(t[[1]]));
  #
  #if (dataType == 'Dense'){
  #  for (i in 1:ncohort){
  #    ave = ave + t[[i]]/ncohort;
  #  }
  #}

  cv = c();
  count = c();
  theFolds =  CreateFolds(unique(ind), k= kFolds)

  for (j in 1:(nbw-1)){
    cv[j]=0;
    count[j]=0;
    #for (i in 1:ncohort){
    for (i in 1:kFolds){
      
      xout= ttn[ ind %in% theFolds[[i]]];
      obs = yyn[ ind %in% theFolds[[i]]];
      xin = ttn[!ind %in% theFolds[[i]]];
      yin = yyn[!ind %in% theFolds[[i]]];
      
      win=rep(1,length(yin));
      #win[ind==i] = NA;        
      #if(dataType=='Dense') {
      #  yyn=(ave*ncohort-t[[i]])/(ncohort-1);
      #  ttn=t[[1]];
      #  win=pracma::ones(1,length(t[[1]]));    
      #  yyn = yyn[order(ttn)]
      #  ttn = sort(ttn)           
      #}  
 
      mu = tryCatch(
        Lwls1D(bw= bw[j], kernel_type = kernel, npoly=npoly, nder= nder, xin = xin, yin= yin, xout=xout, win = win), 
        error=function(err) {
        warning('Invalid bandwidth during CV. Try enlarging the window size.')
        return(Inf)
      })
        
      cv[j] = cv[j]+ sum((obs-mu)^2)
      if(is.na(cv[j])){
        cv[j] = Inf;
      }
      #count[j] = count[j]+1;
    }
  }
  #cv = cv[(count/ncohort>0.90)];
  #bw = bw[(count/ncohort>0.90)];
  if(min(cv) == Inf){
    stop("All bandwidths resulted in infinite CV costs.")
  }
  
  bopt = bw[(cv==min(cv))];
  
  return(bopt)

}

