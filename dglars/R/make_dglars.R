make_dglars <- function(object,control){
	np <- object$np
	nav <-object$nav
	beta <- make_coef(object)
	action <- make_action(beta[-1,],np)
	if(control$algorithm=="pc"){
		ru <- make_ru_pc(object)
		out <- list(call=NULL,family=NULL,np=np,beta=beta,ru=ru,dev=object$dev[1:np],df=object$df[1:np],
					g=object$g_seq[1:np],X=object$X,y=object$y,action=action,conv=object$conv,control=control)
		
	}
	if(control$algorithm=="ccd"){
		out <- list(call=NULL,family=NULL,np=np,beta=beta,dev=object$dev[1:np],df=object$df[1:np],
					g=object$g_seq[1:np],X=object$X,y=object$y,action=action,conv=object$conv,control=control)
	}
	class(out) <- "dglars"
	out
}
