NumDeconv<-function(impulse.matrix,
		    resp.matrix,
		    dose_iv=NULL,
		    dose_po=NULL,
		    deconv.timescale=NULL,
		    explicit.interpolation=20,
		    implicit.interpolation=10,
		    optim.maxit=200
  )
{

##prepare variables

# require(signal)
# require(compiler)

res.out<-list(par=0,par_explicit=0,par_implicit=0)

iv_multiplication_factor<-1

if (!is.null(dose_iv)&& !is.null(dose_po)) {
  iv_multiplication_factor<-dose_po/dose_iv
}

accuracy<-explicit.interpolation
multipl_2<-implicit.interpolation
maxit_optim<-optim.maxit

## explicit interpolation


impulse_orig<-impulse.matrix

resp_orig<-resp.matrix

time_resp_orig<-resp_orig[,1]

time_imp_orig<-impulse_orig[,1]

if (iv_multiplication_factor != 1) {
  
for(i in 1:nrow(impulse_orig)) {
      impulse_orig[i,2]<-impulse_orig[i,2]*iv_multiplication_factor
    }
}
#check deconvolution timescale
if (is.null(deconv.timescale)){
time_orig<-time_imp_orig
}else{
  if(length(deconv.timescale)>2){
    time_orig<-deconv.timescale
    }else{
      time_orig<-seq(from=deconv.timescale[1],to=deconv.timescale[2],by=((deconv.timescale[2]-deconv.timescale[1])/accuracy))
  }
}
lk_row_orig<-length(time_orig)

time<-seq(from=time_orig[1],to=time_orig[lk_row_orig],by=((time_orig[lk_row_orig]-time_orig[1])/accuracy))

lk_row1<-length(time)

#input_interp<-pchip(time_orig,input_orig[,2],time)

impulse_interp<-pchip(time_imp_orig,impulse_orig[,2],time)

resp_interp<-pchip(time_resp_orig,resp_orig[,2],time)

##setting up implicit interpolation
time_2<-seq(from=time_orig[1],to=time_orig[lk_row_orig],by=((time_orig[lk_row_orig]-time_orig[1])/(multipl_2*accuracy)))

resp_interp_2<-pchip(time_resp_orig,resp_orig[,2],time_2)

impulse_interp_2<-pchip(time_imp_orig,impulse_orig[,2],time_2)

lk_row_2<-length(time_2)

##################################################

error_function<-function(input_data){

MSE<-function(vect1,vect2){
result=0
rows_no<-length(vect1)
for(i in 1:rows_no) result<-result +(vect1[i]-vect2[i])^2    
   result<-(result/rows_no)
  return(result)
}

PKconvolution <- function(input,impulse,lk_row){


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
    convolution[i]<-as.numeric(pom)
}
return(convolution)
}
  
  error<-1000

   input_data_2<-pchip(time,input_data,time_2)
  
   try(convol_res<-PKconvolution(input_data_2,impulse_interp_2,lk_row_2))
 
   try(error<-MSE(convol_res,resp_interp_2))

  return(error)  

}

##initialize the input vector values
input_discovered<-vector(length=lk_row1)

multipl<-1/lk_row1

for(i in 1:lk_row1) {
input_discovered[i]<-i*multipl
}


print("START")

error_function_comp<-cmpfun(error_function)

## optim with optim(SANN)

try(fit1 <- optim(
	input_discovered,
	error_function_comp,    
	method="BFGS",
	control=list(trace=1,maxit=maxit_optim)
))
summary(fit1)


#input_discovered<-smooth(fit1$par)
input_discovered<-fit1$par

interpolated<-pchip(time,input_discovered,time_orig)

interpolated_2<-pchip(time,input_discovered,time_2)

out_values<-matrix(nrow=length(interpolated),ncol=2)

out_values_1<-matrix(nrow=lk_row1,ncol=2)

out_values_2<-matrix(nrow=lk_row_2,ncol=2)

#Orig_scale
for(i in 1:lk_row_orig) {
#cat(time_orig[i]," ",interpolated[i],"\n")
out_values[i,1]<-time_orig[i]
out_values[i,2]<-interpolated[i]
}

#Explicit interpolation
for(i in 1:lk_row1) {
out_values_1[i,1]<-time[i]
out_values_1[i,2]<-input_discovered[i]
}


#Implicit interpolation
for(i in 1:lk_row_2) {
#cat(time_2[i]," ",interpolated_2[i],"\n")
out_values_2[i,1]<-time_2[i]
out_values_2[i,2]<-interpolated_2[i]
}


res.out$par <-out_values
res.out$par_explicit <-out_values_1
res.out$par_implicit <-out_values_2

res.out

}
#end NumDdeconv
