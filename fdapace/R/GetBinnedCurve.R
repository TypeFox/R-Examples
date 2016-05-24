GetBinnedCurve <- function(x, y, M = 10, isMnumBin = TRUE, 
    nonEmptyOnly = FALSE, limits = c(min(x),max(x)) ){
  
  # x : 1 * n vector of values to be binned
  # y : 1 * n vector of values corresponding to x
  # M : positive interger denotes the number of bins to be use  or 
  #             positive value of the width of each equidistant bin
  # isMnumBin : use number of bin (TRUE) or bandwith h (FALSE)
  # nonEmptyOnly : output only non-empty bins (TRUE)
  # limits : lower and upper limit of domain x (a0, b0)
  
  # Auxilary function 'GetBins'
  GetBins <-  function(x,y, xx, N){
    count = rep(0,N-1);
    newy = count;
    
    for (i in  2:(N-1)){    
      ids = ((x >= xx[i-1]) &  (x < xx[i]));
      if (all(ids == 0)){
        count[i-1] = 0;
        newy[i-1] = NaN;      
      }
      else{
        count[i-1] =   sum(ids);
        newy[i-1] =  mean(y[ids]); 
      }
    }
    
    # print('GetBins used')
    # for the last bin, include the left and right end point
    ids =  ((x >= xx[i]) &(x <= xx[i+1]));
    if (all(ids == 0)){
      count[i] = 0;
      newy[i] = NaN;
    }else {
      count[i] = sum(ids);
      newy[i] = mean(y[ids]);
    }
    
    zList = list(newy = newy, count = count)
    return( zList );     
  }
  
 
  # Function 'GetBinnedCurve' starts here  
  if( M<0 ){  
    stop("GetBinnedCurve is aborted because the argument: M is invalid!\n");  
  }
  
  if(!isMnumBin){
    # Use bandwidth h
    h = M;
    if ( h >= diff(range(x))){
      res = getResMisOne(x,y,h)
      return(res);
    }
    
    xx =  seq(limits[1],limits[2],by=h)
    M = length(xx)-1;
    N = M + 1;
    midpoint = xx[1:(N-1)]+(h/2);
    print(midpoint)
    zList = GetBins(x,y,xx,N);
    newy = zList$newy
    count = zList$count
    
    if (nonEmptyOnly){
      midpoint = midpoint[count > 0];
      newy = newy[count > 0];
      count = count[count > 0];
      M = length(midpoint);
    }
    
    res = list(midpoint, newy, count, M, h)
    return(res);
  }  else {   
    # Use number of bins M
    M = ceiling(M)
    if (M == 1){
      res = getResMisOne(x,y);
      return(res);
    } else {
      xx =  seq(limits[1],limits[2],length.out=(1+M))
      N = length(xx);
      h = xx[2]-xx[1];
      midpoint = xx[1:(N-1)]+(h/2);
      
      zList = GetBins(x,y,xx,N);
      newy = zList$newy
      count = zList$count
      
      if (nonEmptyOnly){
        midpoint = midpoint[count > 0];
        newy = newy[count > 0];
        count = count[count > 0];
        M = length(midpoint);
      }      
            
      res = list(midpoint = midpoint, newy = newy, count = count, numBin = M, binWidth = h)
      return(res);  
    }
    
  }
}

 # Auxilary function 'GetResMisOne'
getResMisOne <-  function(x, y, h = diff(range(x))){
  r = h;
  M = 1;
  midpoint = r*0.5;
  count = length(x);
  newy = mean(y);
  zList = list(midpoint, newy, count, M, h)
  return( zList );
 
}
  
