bw.dboot1 <- function(y,sig, h0='dnrd', error='normal',grid=100,ub=2)
  {
    homo = FALSE;    if(length(sig)==1) homo=TRUE;
	if (is.character(h0)) {
        h0 <- switch(tolower(h0), dnrd = bw.dnrd(y, sig, error = error), 
			stop("Unknown bandwidth rule!!!"))
    }
    if(!homo){
      if(length(y)!=length(sig))stop("Different length of 'y' and the SD(s).");
      if(any(is.na(y))||any(is.na(sig))){
        sele = (!is.na(y))&(!is.na(sig));
        y=y[sele];sig=sig[sele];
      }
    }else{
      if(is.na(sig)||is.null(sig)) stop("SD(s) can not be empty!");
    }
    if(length(y)<3){stop("The sample size is too small!");}else{n=length(y);}
    if(!is.numeric(grid)) grid=100;
    if(!is.numeric(ub)) ub=2;
    if(ub<0) ub=abs(ub);
    if(ub<1) up=1/ub;
    if(ub==1){warning("Invalid upper boundary of searcing grid!");ub=2}
    if(length(sig)==1){
      type=switch(substr(tolower(error),1,3),
        nor = 1,lap = 2,stop("This error type is not supported yet!"));
    }else{
      type=switch(substr(tolower(error),1,3),
        nor = 3,lap = 4,stop("This error type is not supported yet!"));
    }
    n=length(y);
    result=.C("bwBoot1",as.double(h0),as.integer(n),as.integer(type),
      as.double(y),as.double(sig),as.integer(grid),as.double(ub))[[1]];
    return(result);
  }
