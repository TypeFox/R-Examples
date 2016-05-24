KERNEL <- function(X,type,dim=2,Qcomp="spher"){
    if(dim!=2) stop("Function KERNEL only defined for bivariate cases... check input data")
	
    if(type=="gaus") return(KBivN(X))
    if(type=="quar") return(KBivQ(X,Qcomp))
    
}
