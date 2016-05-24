
MCGHD <- function(data=NULL, gpar0=NULL, G=2, max.iter=100, eps=1e-2,  label=NULL, method="km",scale=TRUE,nr=10, modelSel="AIC" ) {
	data=as.matrix(data)
    if( scale==TRUE){
        data=scale(data)}
    	pcol=ncol(data)
        #if (nrow(data)<((G-1)+G*(4*pcol+2+pcol*(pcol-1)/2)))stop('G is too big, number of parameters > n')
	if (is.null(data)) stop('data is null')
	if (nrow(data) == 1) stop('nrow(data) is equal to 1')
	#if (ncol(data) == 1) stop('ncol(data) is equal to 1; This function currently only works with multivariate data p > 1')
	if (any(is.na(data))) stop('No NAs allowed.')
	if (is.null(G)) stop('G is NULL')
	#if ( G < 1) stop('G is not a positive integer')
	if (  max.iter< 1) stop('max.iter is not a positive integer')
	
    if(modelSel=="BIC"){

    bico=-Inf
    t=length(G)
    BIC=matrix(NA,t,1)
    cont=0
	for(b in 1:t){
    mo=try(MainMCGHD(data=data, gpar0=gpar0, G=G[b], max.iter, eps,  label, method,nr=nr),silent = TRUE)
    cont=cont+1
    if(is.list(mo)){
        bicn=mo$BIC
    BIC[cont]=bicn}
    else{bicn=-Inf
        BIC[cont]=NA}
    if(bicn>bico){
        bico=bicn
        sg=G[b]
        model=mo
    }
    }
     val=list(index=BIC,model=model)
    cat("The best model (BIC) for the range of  components used is  G = ", sg,".\nThe BIC for this model is ", bico,".",sep="")
        return(val)}
   else if(modelSel=="ICL"){
        
        iclo=-Inf
        t=length(G)
        ICL=matrix(NA,t,1)
        cont=0
        for(b in 1:t){
            mo=try(MainMCGHD(data=data, gpar0=gpar0, G=G[b], max.iter, eps,  label, method,nr=nr),silent = TRUE)
            cont=cont+1
            if(is.list(mo)){
                icln=mo$ICL
                ICL[cont]=icln}
            else{icln=-Inf
                ICL[cont]=NA}
            if(icln>iclo){
               iclo=icln
                sg=G[b]
                model=mo
            }
        }
        val=list(index=ICL,model=model)
        cat("The best model (ICL) for the range of  components used is  G = ", sg,".\nThe ICL for this model is ", iclo,".",sep="")
        return(val)}
   else if(modelSel=="AIC3"){
       
       iclo=-Inf
       t=length(G)
       AIC3=matrix(NA,t,1)
       cont=0
       for(b in 1:t){
           mo=try(MainMCGHD(data=data, gpar0=gpar0, G=G[b], max.iter, eps,  label, method,nr=nr),silent = TRUE)
           cont=cont+1
           if(is.list(mo)){
               icln=mo$AIC3
               AIC3[cont]=icln}
           else{icln=-Inf
               AIC3[cont]=NA}
           if(icln>iclo){
               iclo=icln
               sg=G[b]
               model=mo
           }
       }
       val=list(index=AIC3,model=model)
       cat("The best model (AIC3) for the range of  components used is  G = ", sg,".\nThe AIC3 for this model is ", iclo,".",sep="")
       return(val)}
   
   else {
       
       iclo=-Inf
       t=length(G)
       AIC=matrix(NA,t,1)
       cont=0
       for(b in 1:t){
           mo=try(MainMCGHD(data=data, gpar0=gpar0, G=G[b], max.iter, eps,  label, method,nr=nr),silent = TRUE)
           cont=cont+1
           if(is.list(mo)){
               icln=mo$AIC
               AIC[cont]=icln}
           else{icln=-Inf
               AIC[cont]=NA}
           if(icln>iclo){
               iclo=icln
               sg=G[b]
               model=mo
           }
       }
       val=list(index=AIC,model=model)
       cat("The best model (AIC) for the range of  components used is  G = ", sg,".\nThe AIC for this model is ", iclo,".",sep="")
       return(val)}
	}

