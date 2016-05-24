#' Plot for preprocessed spectra using derivative, wavelet  an contiunuum removal transformation#'
#' @author Andrew Sila \email{asila@cgiar.org} and Tomislav Hengl \email{tom.hengl@wur.nl}

plot.trans <-function(x,...){
	raw<-x$raw;
	trans<-x$trans;
	tr<-x$transformation;
	dev.new(width=10,height=7);
	par(mfrow=c(2,1));
	waveb<-as.numeric(colnames(raw));
	plot(raw[1,]~waveb,type="l",ylim=c(min(raw),max(raw)),xlab="Wavebands",ylab="Absorption or Reflection",main="Raw spectra");
		for(i in 2:nrow(raw)){
		lines(raw[i,]~waveb);
		}
	if(tr!="wavelet transformed"){waveb<-as.numeric(colnames(trans));xl="Wavebands";yl="Absorption or Reflection"};
	if(tr=="wavelet transformed"){waveb<-c(1:128);xl="Wavelet coefficients from level 3";yl<-"Value wavelet coefficient"};
	if(tr=="derivative"){te<-"Derivative spectra"};
	if(tr=="continuum removed"){te<-"Continuum removed spectra"};
	if(tr=="wavelet transformed"){te<-"Wavelet transformed spectra"};
		plot(trans[1,]~waveb,type="l",ylim=c(min(trans),max(trans)),xlab=xl,ylab=yl,main=te);
		for(i in 2:nrow(raw)){
		lines(trans[i,]~waveb);
		}
	}
