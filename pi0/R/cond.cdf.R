cond.cdf=function(p.eval,ncp,test=c("t","z"),alternative=c("two.sided","less","greater"),
    df=if(test=="z")Inf else df,keep.cdf=NULL,suppressWarnings=TRUE) #,save=FALSE)
{
# conditional cdf of the p-value evaluated at p.eval, given ncp. 
# The output dim is c( length(p.eval), length(ncp) )

    p.eval=as.vector(p.eval)
    ncp=as.vector(ncp)
    Npeval = length(p.eval)  ;
    Nncp=length(ncp)
    MAXncp=ceiling(max(ncp))
    test=match.arg(test)
    if(test=='z')df=Inf
    alternative=match.arg(alternative)

    cdfp.prefix=paste("Npeval_",Npeval,"__MAXncp",MAXncp,"__Nncp_",Nncp,"__test_",test,
                     "__alternative_",alternative,"__df_",df, sep='')

    #exist.hiddenEnv=exists(".pi0cdfp", envir=globalenv())
	if(is.logical(keep.cdf)){
		if(isTRUE(keep.cdf)){
			warning("As of version 1.4-0, 'keep.cdf=TRUE' is no longer supported to comply with new CRAN policies. Please supply an environment")
		}
		keep.cdf= NULL
	}
	if(!is.null(keep.cdf) && !is.environment(keep.cdf)){
		stop('"keep.cdf" needs to be an environment when it is non-null')
	}
    #if(exist.hiddenEnv) {	
    if(!is.null(keep.cdf)  ) {	
        ## .pi0cdfp is a hidden environment in the global environment. 
        ## objects in this environment should be conditional cdf of p-values 
        ## given a vector of non-centrality parameters (ncp) 
        ## this section of code check if the desired result was already 
        ## kept in the environment or not


        #hiddenEnv=get('.pi0cdfp', envir=globalenv())
        hiddenEnv=keep.cdf
        existing.cdfps=grep(paste("^",cdfp.prefix,sep=''),ls(envir=hiddenEnv),value=TRUE)
        if(length(existing.cdfps)>0){
            check.attr=function(x) {
                x=get(x,envir=hiddenEnv)
                if(isTRUE(all.equal(attr(x,"p.eval"),p.eval)) &&
                   isTRUE(all.equal(attr(x,"ncp"),ncp)) &&
                   isTRUE(all.equal(attr(x,"test"),test)) &&
                   isTRUE(all.equal(attr(x,"alternative"),alternative)) &&
                   isTRUE(all.equal(attr(x,"df"),df))
                ) TRUE else FALSE
            }
            check.rslt=sapply(existing.cdfps,check.attr)
            if(sum(check.rslt)==1)
                return(get(existing.cdfps[[which(check.rslt)]],envir=hiddenEnv))
        }
    }

#    add.cdfp.to.global=function(x){
#	if( "cdfp"%in%ls(envir=globalenv()) ){
#		cdfp=get("cdfp",envir=globalenv())
#		if(!("list"%in%class(cdfp)))
#			cdfp=list(cdfp)
#    		n=length(cdfp)
#		cdfp[[n+1]]=x
#	}else cdfp=x
#	assign("cdfp",cdfp,envir=globalenv())
#    }

#    fname=paste("cdfp","Npeval",Npeval,test,alternative,"df",
#        ifelse(is.finite(df),df,"Inf"),"RData",sep=".")		
#        ### FIXME: the length(ncp) needs to be in the filename.
#    if(fname%in%dir()){
#        load(fname);cat(fname,"loaded",fill=TRUE)
#        attr(cdfp,"Npeval")=Npeval
#        attr(cdfp,"test")=test
#        attr(cdfp,"alternative")=alternative
#        attr(cdfp,"df")=df
#	if(keep.cdf)add.cdfp.to.global(cdfp)
#        return(cdfp)
#    }

    if(isTRUE(suppressWarnings)){
        old.warning.level=options('warn')
        options(warn=-1)
		on.exit(options(old.warning.level))
    }
    if(test=="z"){
        if(alternative=="greater"){
            z =  tcrossprod(qnorm(p.eval, lower.tail=FALSE),rep(1,Nncp)) ;
            p =  pnorm(z,tcrossprod(rep(1,Npeval),ncp),lower.tail=FALSE ) ;
        }else if (alternative=="less") {
            z =  tcrossprod(qnorm(p.eval),rep(1,Nncp)) ;
            p =  pnorm(z,tcrossprod(rep(1,Npeval),ncp)) ;
        }else if (alternative=="two.sided"){
            z =  tcrossprod(abs(qnorm( 1 - p.eval/2)), rep(1,Nncp));
            ncpmat=tcrossprod(rep(1,Npeval),ncp)
            p = 1-pnorm(z,ncpmat)+pnorm(-z,ncpmat)  ;
        }else {
            stop("unsupported alternative type.")
        }
    }else if(test=='t'){
        if(alternative=="greater"){
            tt=tcrossprod(qt(p.eval,df,lower.tail=FALSE),rep(1,Nncp))
            p=pt(tt,df, tcrossprod(rep(1,length(p.eval)),ncp),lower.tail=FALSE)
        }else if (alternative=="less") {
            tt = tcrossprod(qt(p.eval,df),rep(1,length(ncp))) ;
            p =  pt(tt, df, tcrossprod(rep(1,length(p.eval)),ncp)  ) ;
        }else if (alternative=="two.sided"){
            tt=tcrossprod(abs(qt(1- p.eval/2 ,df)),rep(1,length(ncp)))
            ncpmat=tcrossprod(rep(1,Npeval),ncp)
            p=1-pt(tt,df,ncpmat)+pt(-tt,df,ncpmat)
        }else {
            stop("unsupported alternative type.")
        }
    }else{
            stop('test type not supported yet.')
    }
#    if(isTRUE(suppressWarnings)){
#        options(warn=old.warning.level$warn)
#    }

    attr(p,"p.eval")=p.eval
    attr(p,"ncp")=ncp
    attr(p,"test")=test
    attr(p,"alternative")=alternative
    attr(p,"df")=df

	#if(keep.cdf){
	if(!is.null(keep.cdf)) {
        #if(!exist.hiddenEnv) hiddenEnv=new.env()	
		hiddenEnv = keep.cdf

        #if(exist.hiddenEnv && length(existing.cdfps)>0){
		if(length(existing.cdfps)>0){
            rslts=sapply(strsplit(existing.cdfps,"__"),tail,n=1)
            rslt.number=max(as.numeric(sapply(strsplit(rslts,"_"),tail,n=1)))+1
        }else rslt.number=1
        assign(paste(cdfp.prefix,"__rslt_",rslt.number,sep=''),p,envir=hiddenEnv)
        #if(!exist.hiddenEnv)assign(".pi0cdfp",hiddenEnv,envir=globalenv())
    }
                
        
#    if(is.logical(save) & !save) return(p)
#    if (is.logical(save)) save=getwd()
#    cdfp=p
#    save(cdfp,file=paste(sub("[/\\\\]$","",save),"/",fname,sep=""))
#    return(cdfp)
    p
}
