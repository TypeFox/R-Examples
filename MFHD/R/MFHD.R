derivcurves<-function(y1,method="bspline") fdata.deriv(fdata(y1),nderiv=1,method=method,class.out="fdata")$data
MFHD<-function(y1,y2,alpha=0.125,Beta=0.5,Time=NULL){
	y1<-data.matrix(y1)
	y2<-data.matrix(y2)
	n1<-nrow(y1)
	n2<-nrow(y2)
	T1<-ncol(y1)
	T2<-ncol(y2)
	no.y1<-sum(complete.cases(y1))
	no.y2<-sum(complete.cases(y2))
	if(no.y1!=n1)	stop("y1 contains missing cases.")
	if(no.y2!=n2)	stop("y2 contains missing cases.")
	if(n2!=n1)	stop("y1 and y2 do not have the same number of observations.")
	if(T2!=T1)	stop("y1 and y2 do not have the same number of measurements.")
	T<-T1
	if(!is.numeric(alpha) & !is.null(alpha))	stop("alpha should be NULL or a numeric.")
	if(!is.null(alpha)){
	      	if(alpha<=0 | alpha>1)  stop("If alpha is numeric, it should be in ]0-1].")
	}
     	if(!is.numeric(Beta))	stop("Beta should be a numeric.")
	if(Beta<0 | Beta>1)	stop("Beta should be in [0-1].")
	if(is.null(Time))	Time<-1:T
	if(!is.numeric(Time))	stop("Time should be a numeric vector.")
	if(length(Time)!=T)	stop("Time should contain T elements.")
	dTime<-diff(c(Time[1],Time,Time[T]),lag=2)
	if(min(dTime)<=0)	stop("Time intervals should all be strictly positive.")
	n<-nrow(y1)
	T<-ncol(y1)
	err.loc<-c()
	ybiv<-matrix(0,n,2)
	depths.time<-matrix(NA,n,T)
	weights<-matrix(1,1,T)
	disp<-matrix(0,2,T)
 	MFHDmedian=matrix(0,2,T)
	loc.outl=matrix(0,n,T);
	for (j in 1:T){
		ybiv[,1]<-y1[,j]
		ybiv[,2]<-y2[,j]	
		for (s in 1:n) depths.time[s,j]<-depth(ybiv[s,],ybiv)
		z<-NA
		z1<-matrix(NA,1,3)
		l<-0
		if(!is.null(alpha)){
			while(nrow(z1)<3){
				n0<-ceiling(n*alpha)-l
				if(n0<1)	stop("No valid depths region. Alpha needs to be changed.")
				z0<-try(isodepth(ybiv,mustdith=TRUE,dpth=n0,output=TRUE,factor=1)[[1]])
				if(is.numeric(z0))	z1<-unique(z0)
				l<-l+1
			}		
			w1<-try(deldir(z1[,1],z1[,2])$del.area)		
			if(is.numeric(w1)){	
				weights[j]<-w1
			} else {
				weights[j]<-0
			}
		}
		z<-try(MFHD_bagplot(ybiv[,1],ybiv[,2],verbose=FALSE,factor=3))
		if(is.list(z)){
			if(!is.null(z$outliers))	loc.outl[z$outliers,j]=1
			MFHDmedian[,j]=z$center
		} else {
			if(j>1){
				MFHDmedian[,j]<-MFHDmedian[,j-1]
				txt2<-". Using the MFHD median from the previous time period instead."
			} else {
				MFHDmedian[,j]<-colMedians(ybiv) 
				txt2<-". Using coordinate-wise medians instead."
			}
			err.loc<-c(err.loc,j)
			txt1<-"Could not compute Halfspace median for time point #"		
			print(paste0(txt1,j,txt2))
		}		
	}
	weights<-weights*dTime
	weights<-weights/sum(weights)
	depths<-depths.time%*%t(weights)
      	select<-which(depths>=quantile(depths,Beta))
	if(length(select)>1){
		region1<-y1[select,]
		disp[1,]<-rowDiffs(colRanges(region1))
	      	region2<-y2[select,]
		disp[2,]<-rowDiffs(colRanges(region2))
	}
	list(MFHDdepth=t(depths),MFHDmedian=MFHDmedian,weights=weights,disp=disp,loc.outl=loc.outl)
}
