magmap<-
function(data, lo=0, hi=1, flip=FALSE, range=c(0,2/3), type='quan', stretch='lin', stretchscale=1, bad=NA, clip=''){
if(stretch=='log'){
	good=is.na(data)==FALSE & is.nan(data)==FALSE & is.infinite(data)==FALSE & is.null(data)==FALSE & data>0
	if(length(which(good))==0){stop('There is no numeric data with a value greater than 0!')}
	}else{
	good=is.na(data)==FALSE & is.nan(data)==FALSE & is.infinite(data)==FALSE & is.null(data)==FALSE
	if(length(which(good))==0){stop('There is no numeric data!')}
	}
if(type=='quan'){
	lo=quantile(data[good],lo)
	hi=quantile(data[good],hi)
}
if(type=='num'){
	lo=lo
	hi=hi
}
if(type=='sig'){
	lo=quantile(data[good],pnorm(lo))
	hi=quantile(data[good],pnorm(hi))
}
if(type=='rank'){
	lo=1
	hi=length(data[good])
	data[order(data[good])]=lo:hi
}
if(stretch=='log' & lo==0){stop('lo value is 0 and stretch=\'log\'- this is not allowed!')}
if(stretch=='log' & hi==0){stop('hi value is 0 and stretch=\'log\'- this is not allowed!')}
if(lo>hi){stop('lo>hi is not allowed')}
if(lo==hi){data=rep((range[2]+range[1])/2,length(data))}
if(lo<hi){
	if(stretch=='log'){
	  lo=log10(lo)
	  hi=log10(hi)
    data=suppressWarnings(log10(data))
	  #if(is.finite(bad)){bad=log10(bad)}
	}
  if(stretch=='atan'){
    lo=atan(lo*stretchscale)
    hi=atan(hi*stretchscale)
    data=atan(data*stretchscale)
    #if(is.finite(bad)){bad=atan(bad*stretchscale)}
  }
	if(stretch=='asinh'){
	  lo=asinh(lo*stretchscale)
	  hi=asinh(hi*stretchscale)
    data=asinh(data*stretchscale)
	  #if(is.finite(bad)){bad=asinh(bad*stretchscale)}
	}
  losel=data<lo; hisel=data>hi
	data[losel]=lo; data[hisel]=hi
	data=data-lo
	data=range[1]+(data*(range[2]-range[1])/(hi-lo))
	if(flip){data=range[2]-data+range[1]}
  if(clip=='NA'){data[losel]=NA;data[hisel]=NA}
	#if(is.finite(bad)){
  #  bad=bad-lo
  #  bad=range[1]+(bad*(range[2]-range[1])/(hi-lo))
  #  if(bad<range[1]){bad=range[1]}
  #  if(bad>range[2]){bad=range[2]}
  #  if(flip){bad=range[2]-bad+range[1]}
	#}
}
data[! good]=bad
return(list(map=data,datalim=c(lo,hi),maplim=range,loclip=length(which(data[good]==range[1]))/length(data[good]),hiclip=length(which(data[good]==range[2]))/length(data[good])))
}

