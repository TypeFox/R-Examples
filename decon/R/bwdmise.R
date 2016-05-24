bw.dmise <- function(y,sig,error='normal',kernel='support',grid=100,ub=2)
  {
    hsupp1 <- function(h1,s2bar,error,Rfx,n,grid=100,ub=2){
      switch(substr(tolower(error),1,3),
             lap = .C("SuppLap1",as.integer(n),as.double(Rfx),
               as.double(s2bar), as.double(h1),as.double(grid),as.double(ub))[[4]],
             nor = .C("SuppNorm1",as.integer(n),as.double(Rfx),
               as.double(s2bar), as.double(h1),as.double(grid),as.double(ub))[[4]],
             stop("This error type is not supported yet!") 
             );
    }
    hsupp2 <- function(h1,sig,error,Rfx,n,grid=100,ub=2){
      switch(substr(tolower(error),1,3),
             lap = .C("SuppLap2",as.integer(n),as.double(Rfx),
               as.double(sig), as.double(h1),as.double(grid),as.double(ub))[[4]],
             nor = .C("SuppNorm2",as.integer(n),as.double(Rfx),
               as.double(sig), as.double(h1),as.double(grid),as.double(ub))[[4]],
             stop("This error type is not supported yet!")
             );
    }
    hnorm1 <- function(h1,s2bar,error,Rfx,n,grid=100,ub=2){
      switch(substr(tolower(error),1,3),
             lap = .C("NormLap1",as.integer(n),as.double(Rfx),
               as.double(s2bar), as.double(h1),as.double(grid),as.double(ub))[[4]],
             nor = .C("NormNorm1",as.integer(n),as.double(Rfx),
               as.double(s2bar), as.double(h1),as.double(grid),as.double(ub))[[4]],
             stop("This error type is not supported yet!")
             );
    }
    hnorm2 <- function(h1,sig,error,Rfx,n,grid=100,ub=2){
      switch(substr(tolower(error),1,3),
             lap = .C("NormLap2",as.integer(n),as.double(Rfx),
               as.double(sig), as.double(h1),as.double(grid),as.double(ub))[[4]],
             nor = .C("NormNorm2",as.integer(n),as.double(Rfx),
               as.double(sig), as.double(h1),as.double(grid),as.double(ub))[[4]],
             stop("This error type is not supported yet!")
             );
    }
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
      s2y = var(y);
    }else{
      if(is.na(sig)||is.null(sig)) stop("Variance(s) can not be empty!");
      s2bar=sig;sbar=sqrt(s2bar);s2y = var(y);
    }
    if(s2y-s2bar<=0){
      stop("Berkson Model should be considered!");
    }else{
      Rfx = 0.09375*(s2y-s2bar)^(-2.5)*pi^(-.5); # R(f'')/4
    }
    if(length(y)<3){
      stop("Data set is too small!");
    }else{
      n=length(y);
    }
    h1 = bw.dnrd(y,sig,error=error);
    if(!is.numeric(grid)) grid=100;
    if(!is.numeric(ub)) ub=2;
    if(ub<0) ub=abs(ub);
    if(ub<1) up=1/ub;
    if(ub==1){warning("Invalid upper boundary of searcing grid!");ub=2}
    if(homo){
      result=switch(substr(tolower(kernel),1,4),
        supp = hsupp1(h1,s2bar,error,Rfx,n,grid=grid,ub=ub),
        norm = hnorm1(h1,s2bar,error,Rfx,n,grid=grid,ub=ub),
        stop("This kernel type is not supported yet!")
        );
    }else{
      result=switch(substr(tolower(kernel),1,4),
        supp = hsupp2(h1,sig,error,Rfx,n,grid=grid,ub=ub),
        norm = hnorm2(h1,sig,error,Rfx,n,grid=grid,ub=ub),
        stop("This kernel type is not supported yet!")
        );
    }
    return(result);
  }
