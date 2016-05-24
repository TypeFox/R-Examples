test_dim<-function(S,yv,k,link=1,disc=0,difl=0,multi0=1:J,multi1,tol=10^-10,disp=FALSE){
	
# disp = to display log-likelihood evolution step by step

# test_dim<-function(S,yv,k,link=1,disc=0,difl=0,multi0=1:J,multi1,tol=10^-10,disp=FALSE)
    J = dim(S)[2]
	out0 = est_multi_poly(S,yv,k,link=link,disc=disc,difl=difl,multi=multi0,tol=tol,disp=disp)  # unidimensional model
	out1 = est_multi_poly(S,yv,k,link=link,disc=disc,difl=difl,multi=multi1,tol=tol,disp=disp)  # bidimensional model
	dev = 2*(out1$lk-out0$lk)
	df = out1$np-out0$np
	pv = 1-pchisq(dev,df)
	out = list(out0=out0,out1=out1,dev=dev,df=df,pv=pv,table=table,call=match.call())	
	class(out) = "test_dim"
	out

}