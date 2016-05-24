RivivcA<-function(known.data,
		 impulse.data,
		 second.profile.data,
		 dose_iv=NULL,
		 dose_po=NULL,
		 mode="deconv",
		 explicit.interp=20,
		 implicit.interp=10,
		 optimization.maxit=200){

x<-known.data[,2]
known_time<-known.data[,1]
out_val<-list(regression=0,numeric=0)		 		 

if (mode=="deconv") {
#deconvolution mode
wynik<-NumDeconv(impulse.data,
		    second.profile.data,
		    dose_iv=dose_iv,
		    dose_po=dose_po,
		    deconv.timescale=known_time,
		    explicit.interpolation=explicit.interp,
		    implicit.interpolation=implicit.interp,
		    optim.maxit=optimization.maxit)		    
}else if (mode=="conv"){
#convolution_mode
wynik<-NumConv(impulse.data,
		    second.profile.data,
		    conv.timescale=known_time,
		    explicit.interpolation=explicit.interp
		    )


}

y<-wynik$par[,2]

regr<-lm(y~x)

out_val$regression<-regr
out_val$numeric<-wynik

return(out_val)

}
