sa <- function(par, fn, method=c("FAST"), ...,xval=NULL){
	#example:
	#example_model1<-function(p,x){
	#         return(p[1]*x+p[2]*(1-x))
	#}
	#par <-  matrix(c(0,0,0,1,2,2),ncol=2)
	# sa(par=par, fn=example_model1, x=0.5)
	# sa(par=matrix(c(0,0,0,1,1,1),ncol=2), fn=sens.model1, x=seq(0,1,by=0.1), xval=seq(0,1,by=0.1))

	stopifnot(mode(par)=="numeric")
	stopifnot(NCOL(par)==2)
	match.arg(method, several.ok=FALSE)
	stopifnot(mode(fn)=="function")
      if(method=="FAST"){
	     paras<-fast_parameters(minimum=par[,1],maximum=par[,2])
	     model_results <- apply(paras, 1, fn, ...)
	     model_results
	     if(!is.null(dim(model_results)) & is.null(xval)){
		     xval=1:NROW(model_results)
             }
	     if(is.null(xval)){
		     sensitivity <- sensitivity(x=model_results, numberf=NROW(par), make.plot=FALSE)
	     } else {
		     sensitivity <- sensitivity_rep(data=model_results, direction=1, numberf=NROW(par), xval=xval)
             }
	     return(sensitivity)

	      }
}

