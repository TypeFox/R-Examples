### Author of original python script of ARSER: Rendong Yang
### Email: cauyrd@gmail.com
### Associated literature: Rendong Yang and Zhen Su. Bioinformatics. 26(12):i168-74 (2010). 
### Website: http://bioinformatics.cau.edu.cn/ARSER/
### R script of ARSER: Gang Wu
### Email: wggucas@gmail.com
### Lab: John Hogenesch's lab in Perelman School of Medicine at University of Pennsylvania (http://hogeneschlab.org/)
###======================================================================================================================================
detrend_linear<-function(y)
{
  x=0:(length(y)-1);                                       ##'y' should be a vector type.
  C=cov(cbind(x,y))*(length(y)-1)/length(y);               ##'(length(y)-1)/length(y)' is used to keep the value as same at that calculated by 'cov(x,y,bias=1)' in the original colde of ARSER.
  b=C[1,2]/C[1,1];
  a=mean(y) - b*mean(x);
  return ( y- (b*x+a) );
}
###-----------------------------------------------------------------------------
###applies a Savitzky-Golay filter
savitzky_golay<-function(data, kernel=11, order=4)         
{
  ##input parameters:
  ##data => data as a vector type
  ##kernel => a positiv integer > 2*order giving the kernel size
  ##order => order of the polynomal
  ##returns smoothed data as a vector type
  ##-----------------------
  if ( (round(kernel)!=kernel) | (round(order)!=order) | (kernel<1) | (order<1) ) 
  { stop("The input values of kernel and order in 'savitzky_golay' should be positive integers."); }
  if (kernel %% 2 != 1)
  { stop("The input value of kernel in 'savitzky_golay' should be an odd number"); }
  if (kernel < order+2)
  { stop("The kernel value is to small for the polynomals, and it should be larger than (order+2)"); }
  ##-----------------------
  order_range=0:order;                                     #a second order polynomal has 3 coefficients
  half_window=(kernel-1)%/%2;
  b="b";
  for (k in (-half_window):(half_window) )
  {
    for (i in order_range)
    {
      if (b[1]=="b") {
       b=k**i;
       } else {
       b=c(b,k**i);
       } 
    }
  }
  b=matrix(b,nrow=half_window*2+1,ncol=length(order_range),byrow=TRUE);
  ##-----------------------
  ##'MPinv()' based on 'gnm' package in R can calculate the 'Moore-Penrose pseudo-inverse matrix' as 'numpy.linalg.pinv' in python
  #library(gnm);
  m=MPinv(b);
  m=m[1,];
  window_size=length(m);
  half_window=(window_size-1)%/%2; 
  offsets=(-half_window):half_window;                      #pre-compute the offset values for better performance
  offset_data=cbind(offsets,m);
  firstval=data[1];                                        #temporary data, extended with a mirror image to the left and right
  lastval=data[length(data)];
  leftpad=rep(0,half_window)+2*firstval;
  rightpad=rep(0,half_window)+2*lastval;
  leftchunk=data[2:(1+half_window)];
  leftpad=leftpad-rev(leftchunk);
  rightchunk=data[(length(data)-half_window):(length(data)-1)];
  rightpad=rightpad-rev(rightchunk);
  data=c(leftpad,data);
  data=c(data,rightpad);
  ##-----------------------
  smooth_data="smooth_data";
  for (i in half_window:(length(data) -half_window -1) )
  {
    value=0;
    for (j in 1:nrow(offset_data))
    { 
      offset=offset_data[j,1];
      weight=offset_data[j,2];
      value=value+weight*data[i+offset+1];
    }
    smooth_data=c(smooth_data,value);
  }
  smooth_data=smooth_data[2:length(smooth_data)];
  smooth_data=as.numeric(smooth_data);
  return (smooth_data);
}
###-----------------------------------------------------------------------------
###calculate model log-likelihood and two information criteria(AIC and BIC criterion values)                               
LL<-function(e,nobs,ncoef)                                 
{
  LL= -(nobs*1/2)*(1+log(2*pi)) - (nobs/2)*log((e%*%e)/nobs);
  aic= -2*LL/nobs + (2*ncoef/nobs);
  bic= -2*LL/nobs + (ncoef*log(nobs))/nobs;
  LL.out=c(LL,aic,bic);
  names(LL.out)<-c("ll","aic","bic");
  return (LL.out);
}
###-----------------------------------------------------------------------------
###multi-variate regression using OLS
OLS<-function(x, y, x_varnm='', y_varnm ='y')              
{
  ## y = dependent variable; y is a vector type;
  ## y_varnm = string with the variable label for y;y_varnm is a vector type containing one string.
  ## x = independent variables, note that a constant is added by default;x is a matrix type.
  ## x_varnm = string or list of variable labels for the columns of x; x_varnm is a vector type.
  ##-----------------------
  self.y=y;
  self.y_varnm=y_varnm;
  self.x=cbind(rep(1,nrow(x)),x);
  self.x_varnm=c("const",x_varnm);
  dimnames(self.x)[[2]]<-self.x_varnm;
  self.xT=t(self.x);
  self.nobs=length(self.y);
  self.ncoef=ncol(self.x);
  ##-----------------------
  self.inv_xx=solve( self.xT%*%self.x );
  xy=self.xT%*%self.y;
  self.b=self.inv_xx%*%xy;                                 #change 1-column matrix to vector type
  self.b=as.vector(self.b);
  self.df_e=self.nobs - self.ncoef;
  self.df_r = self.ncoef - 1;
  self.e=self.y - self.x%*%self.b;
  self.e=as.vector(self.e)                                 #change 1-column matrix to vector type
  self.sse=(self.e%*%self.e)/self.df_e;
  self.sse=as.vector(self.sse);                            #change 1-column matrix to vector type
  self.se=sqrt( diag(self.sse*self.inv_xx) );
  self.t=self.b/self.se;
  self.p=( 1 - pt(abs(self.t),self.df_e) )*2;
  e.adfac=(length(self.e)-1)/length(self.e);
  y.adfac=(length(self.y)-1)/length(self.y);
  self.R2=1-(var(self.e)*e.adfac)/(var(self.y)*y.adfac);   #warning: self.e and self.y should be vector type;                      
  self.R2adj=1-(1-self.R2)*((self.nobs-1)/(self.nobs-self.ncoef));
  self.F=(self.R2/self.df_r) / ((1-self.R2)/self.df_e);
  self.Fpv=1-pf(self.F,self.df_r, self.df_e);
  self.LL<-LL(e=self.e,nobs=self.nobs,ncoef=self.ncoef);
  OLS.out=list("b"=self.b,"R2"=self.R2,"R2adj"=self.R2adj,"F"=self.F,"Fpv"=self.Fpv,
               "ll"=self.LL["ll"],"aic"=self.LL["aic"],"bic"=self.LL["bic"],"t"=self.t,"tp"=self.p);
  return (OLS.out);
}
###======================================================================================================================================
###estimate possible cycling period of time series by AR spectral estimation
###is_filter=True, stands for using smoothed time-series values by HR
get_period<-function(per.x,per.dt_y,per.delta,per.pID="",is_filter=TRUE,ar_method="mle")
{
  self.x=per.x;
  self.dt_y=per.dt_y;
  self.delta=per.delta;
  num_freq_mese=500;
  set_order=24/self.delta;                                 #set the parameter of 'order' (the order of the AR model to be fitted) in spec.ar() 
  if (set_order == length(self.x))                         #If the expression profiles only cover one day[24/interval==length(time points)], 'order' is set to length(time points)/2;
  { set_order=length(self.x)/2; }
  filter_y<-try(savitzky_golay(self.dt_y),silent=TRUE);
  if (inherits(filter_y,"try-error"))
  {
    filter_y=savitzky_golay(self.dt_y, kernel=5, order=2);
  }
  flag=1;
  if (is_filter)
  {
    mese=try(spec.ar(filter_y,n.freq=num_freq_mese, order=set_order,method=ar_method,plot=FALSE),silent=TRUE);
    if ( inherits(mese,"try-error") )                      
    { 
      mese=NA;
      flag=NA;
    }
  } else {
    mese=try(spec.ar(self.dt_y,n.freq=num_freq_mese, order=set_order,method=ar_method,plot=FALSE),silent=TRUE);
    if (inherits(mese,"try-error"))                       
    { 
      mese=NA;
      flag=NA;
    }
  }
  if (is.na(flag))
  {
    periods=NA;
  } else {
    peaks_loc=c(NA,0);
    for (i in 2:(num_freq_mese-1))                         #search all the local peaks of maximum entropy spectrum
    {
      if ( ( mese$spec[i] > mese$spec[i+1] ) & ( mese$spec[i] > mese$spec[i-1] ) )
      {
        peaks_loc=rbind(peaks_loc,c(mese$spec[i],i) );
      }
    }
    if (is.vector(peaks_loc)) {
		periods=NA;
	} else {
		peaks_loc=peaks_loc[2:nrow(peaks_loc),];
		if ( is.vector(peaks_loc) ) 
		{ peaks_loc=matrix(peaks_loc,nrow=1,ncol=2); } 
		peaks_loc.sort=peaks_loc[1:nrow(peaks_loc),1];
		names(peaks_loc.sort)=peaks_loc[1:nrow(peaks_loc),2];
		peaks_loc.sort=sort(peaks_loc.sort,decreasing=TRUE);
		periods=try(1/mese$freq[as.numeric(names(peaks_loc.sort))]*self.delta,silent=TRUE);
		if ( inherits(periods,"try-error") )                   
		{  periods=NA;	}
	}
  }
  return (periods);
}
###-----------------------------------------------------------------------------
harmonic_regression<-function(har.x,har.dt_y,period)
{
  ##dt_y = mesor + sigma( A*cos(2*pi/T*x) + B*sin(2*pi/T*x) ) + error
  self.x=har.x;
  self.dt_y=har.dt_y;
  x=self.x;
  x_varnm_names=NA;
  for (T in period)
  {
    cosx=cos(2*pi/T*self.x);
    sinx=sin(2*pi/T*self.x);
    x=cbind(x,cosx,sinx);
    x_varnm_names=c(x_varnm_names,paste("cos",round(T,2),sep=""),paste("sin",round(T,2),sep=""));
  }
  x=x[,2:ncol(x)];
  x_varnm_names=x_varnm_names[2:length(x_varnm_names)];
  model = OLS(x=x, y=self.dt_y,x_varnm = x_varnm_names,y_varnm = 'y');
  return (model);
}
###-----------------------------------------------------------------------------
###evaluate the best model for each time series
evaluate<-function(eva.x,eva.y,eva.delta,eva.pID="",T_start=20, T_end=28, T_default=24,arsmethods=c("yule-walker","mle","burg"))           
{
  eva.dt_y=detrend_linear(eva.y);
  is_filter=c(TRUE, FALSE);
  ar_method=arsmethods;
  best_model=list("AIC"=1e6,"OLS"=NA,"period"=NA,"filter"=NA,"armethod"=NA);
  for (p1 in is_filter)
  {
    for (p2 in ar_method)
    {
      period=get_period(per.x=eva.x,per.dt_y=eva.dt_y,per.delta=eva.delta,per.pID=eva.pID,is_filter=p1,ar_method=p2);
      period=period[period >= T_start & period <= T_end];
      if (!length(period))
	  {
		p2="default";
		period=T_default;
	  }
	  if (length(period)==1)
      {
        if (is.na(period))
        {
          p2="default";
          period=T_default;
        }   
      }
      m=harmonic_regression(har.x=eva.x,har.dt_y=eva.dt_y,period=period);
      aic=m$aic;                                           #model selection by aic
      if (aic <= best_model$AIC)
      {
        best_model$AIC=aic;
        best_model$OLS=m;
        best_model$period=period;
        best_model$filter=p1;
        best_model$armethod=p2;
      }
    }
  }
  ##-----------------------
  m=best_model$OLS;                                        #record the best model parameters
  self.estimate_period = best_model$period;
  self.amplitude=NA;
  self.phase=NA;
  for (i in 1:length(self.estimate_period) )
  {
    phi=Arg(complex(real = m$b[2*i],imaginary = -m$b[2*i+1]));
    if (phi <= 1e-6) {
      self.phase=c(self.phase,abs(phi)/(2*pi)*self.estimate_period[i]);
    } else {
      self.phase=c(self.phase,(self.estimate_period[i] - phi/(2*pi)*self.estimate_period[i]) );
    }
    self.amplitude=c(self.amplitude,sqrt(m$b[2*i]**2 + m$b[2*i+1]**2) );
  }
  self.amplitude=self.amplitude[2:length(self.amplitude)];
  self.phase=self.phase[2:length(self.phase)];
  self.R2=m$R2;
  self.R2adj=m$R2adj;
  self.pvalue=m$Fpv;
  py.std=sd(eva.y)*sqrt( (length(eva.y)-1)/length(eva.y) );
  self.coef_var=py.std/mean(eva.y);                        #'py.std' is equal to 'np.std(eva.y)' in python.
  ##-----------------------
  eva.out=list("period"=self.estimate_period,"amplitude"=self.amplitude,"phase"=self.phase,"R2"=self.R2,"R2adj"=self.R2adj,
               "coefvar"=self.coef_var,"pvalue"=self.pvalue,"filter"=best_model$filter,"armethod"=best_model$armethod);
  return (eva.out);
}
###======================================================================================================================================
runARS <- function(indata,ARStime,minper=20,maxper=28, arsper=24, arsmet="", releaseNote=TRUE)
{
  #-----------------------
  if (releaseNote)  {
	cat("The ARS is in process from ",format(Sys.time(), "%X %m-%d-%Y"),"\n");
  }
  start=minper;
  end=maxper;
  pvalues <-NA
  #-----------------------
  self.data=indata;
  self.delta=ARStime[2]-ARStime[1];                
  time_points=ARStime;
  idorder<-dimnames(self.data)[[1]];
  dataM=as.matrix(self.data[,2:ncol(self.data)]);
  outID=as.character(self.data[,1]);
  names(outID)<-idorder;
  dimnames(dataM)<-list("r"=idorder,"c"=paste("T",self.delta*(0:(ncol(dataM)-1)),sep="" ) );
  #-----------------------
  expSD<-apply(dataM,1,sd);
  expMEAN<-apply(dataM,1,mean);
  constantID<-names(expSD[expSD == 0]);
  flagV<-rep(1,length(idorder));
  names(flagV)<-idorder;
  flagV[constantID]<-0;
  #-----------------------
  run_start=proc.time();
  set.seed(run_start["elapsed"]);
  header <- c("filter_type","ar_method","period_number","period","amplitude","phase","mean","R_square","R2_adjust","coef_var","pvalue");
  ars.outM <- header;
  ori.op <- options();
  #-----------------------
  ##try 'apply()' in latter version, it may improve the Computational Efficiency
  for (line in idorder )
  {
    if (flagV[line]) {
		d=evaluate(eva.x=time_points,eva.y=dataM[line,],eva.delta=self.delta,eva.pID=line,T_start=start, T_end=end, T_default=arsper,arsmethods=arsmet);
		ars.filter=0;
		if (d$filter)
		{ ars.filter=1; }
		ars.period.num=length(d$period);
		ars.period=paste(d$period,collapse=",");
		ars.amplitude=paste(d$amplitude,collapse=",");
		ars.phase=paste(d$phase,collapse=",");
		ars.out=c(ars.filter,d$armethod,ars.period.num,ars.period,ars.amplitude,ars.phase,mean(dataM[line,]),d$R2,d$R2adj,d$coefvar,d$pvalue);
		pvalues=c(pvalues,d$pvalue);
	} else {
		ars.out=c(rep(NA,4),0,NA,expMEAN[line],rep(NA,3),1);                           
		pvalues=c(pvalues,1);                                #assign it as '1' insead of 'NA' for avoiding error report in 'pi0.est' step
	}
    ars.outM=rbind(ars.outM,ars.out);
  }
  pvalues=pvalues[2:length(pvalues)];
  ars.outM=ars.outM[2:nrow(ars.outM),];
  options(ori.op);
  dimnames(ars.outM)[[1]]=idorder;
  names(pvalues)=idorder;
  #-----------------------
  #header=c("CirID",header,"qvalue","fdr_BH");               ## qvalue is not used in the new version of ARSER
  # pi0=pi0.est(pvalues);                                    #It will report error when analyzing expression profile with equal values among all time points.
  # if (pi0$p0 == 0)                                         #p0 is required for calculating qvalues by qvalue.cal(); p0 will define the largest q-value calculated by qvalue.cal
  # { pi0$p0=0.95; }                                         #If p0 is '0', the calculated q-value is '0' for all pvalues, thus change it to 0.95 is a stategy to reduce false positive?
  # qvalues=qvalue.cal(pvalues[idorder],pi0$p0);
  header=c("CycID",header,"fdr_BH");
  qvalues_BH=p.adjust(pvalues[idorder],"BH");
  ARSoutM=cbind(outID[idorder],ars.outM[idorder,],qvalues_BH);
  dimnames(ARSoutM)<-list("r"=idorder,"c"=header);
  if (releaseNote)  {
	cat("The analysis by ARS is finished at ",format(Sys.time(), "%X %m-%d-%Y"),"\n");
  }
  return(ARSoutM);
}
###======================================================================================================================================
