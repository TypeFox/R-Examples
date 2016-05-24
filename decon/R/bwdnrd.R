bw.dnrd <- function(y,sig,error='normal')
  {
    homo = FALSE;
    sig = sig^2;
    if(length(sig)==1) homo=TRUE;
    if(!homo){
      if(length(y)!=length(sig))stop("Different length of 'y' and the variances.");
      if(any(is.na(y))||any(is.na(sig))){
        sele = (!is.na(y))&(!is.na(sig)); 
        y=y[sele];sig=sig[sele];
      }
      s2bar = mean(sig);sbar=sqrt(s2bar);
    }
    else{
      if(is.na(sig)||is.null(sig)) stop("SD(s) can not be empty!");
      s2bar=sig;sbar=sqrt(s2bar);}
    if(length(y)<3){stop("Data set is too small!");}
    else{n=length(y);}
    result=  switch(substr(tolower(error),1,3),
      lap = (5*sbar^4/n)^(1/9),
      nor = sbar*(log(n)/2)^{-.5},
      stop("This error type is not supported yet!")
      );
    return(result);
  }
