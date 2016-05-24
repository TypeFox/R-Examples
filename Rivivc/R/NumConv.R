NumConv<-function(impulse.matrix,
		    input.matrix,
		    conv.timescale=NULL,
		    explicit.interpolation=1000
		    )
{
# require(signal)
# require(compiler)

##################################################
##prepare variables

res.out<-list(par=0,par_explicit=0)

accuracy<-explicit.interpolation

input_orig<-input.matrix

impulse_orig<-impulse.matrix

#check convolution timescale
if (is.null(conv.timescale)){
  time_orig<-input_orig[,1]
}else{
  if(length(conv.timescale)>2){
    time_orig<-conv.timescale
    }else{
      time_orig<-seq(from=conv.timescale[1],to=conv.timescale[2],by=((conv.timescale[2]-conv.timescale[1])/accuracy))
  }
}

time_imp_orig<-impulse_orig[,1]

lk_row_orig<-nrow(input_orig)

time<-seq(from=time_orig[1],to=time_orig[lk_row_orig],by=((time_orig[lk_row_orig]-time_orig[1])/accuracy))

lk_row1<-length(time)

input1<-pchip(time_orig,input_orig[,2],time)

impulse1<-pchip(time_imp_orig,impulse_orig[,2],time)

##################################################

PKconvolution <- function(input,impulse,lk_row){
# compute convolution for unequal time steps
# based on F. Langenbucher / European Journal of Pharmaceutics and Biopharmaceutics 56 (2003) 429–437


convolution<-vector(length=lk_row)


for(i in 1:lk_row) {
  pom<-0
    for(k in 1:i){
      

      input1<-input[i-k+1]
      
      input2<-0
      if(i!=k){	
       input2<-input[i-k]}


      if (k==1){
      pom<-pom+(impulse[k])/2*(input1-input2)
      }else{
	    pom<-pom+(impulse[k]+impulse[k-1])/2*(input1-input2)
	  }      
    }
    convolution[i]<-pom
}
return(convolution)
}

PKconvolution_comp<-cmpfun(PKconvolution)

convol_res<-PKconvolution_comp(input1,impulse1,lk_row1)

interpolated<-pchip(time,convol_res,time_orig)


out_values<-matrix(nrow=length(interpolated),ncol=2)

out_values[,1]<-time_orig
out_values[,2]<-interpolated

res.out$par<-out_values


out_values<-matrix(nrow=length(time),ncol=2)

out_values[,1]<-time
out_values[,2]<-convol_res

res.out$par_explicit<-out_values

res.out

}
#end of NumConv