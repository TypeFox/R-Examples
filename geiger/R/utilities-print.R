print.constraint.m=function(x, printlen=3, ...){
	cat("matrix representation of unique and shared transitions \n")
	tt=table(x)
	if(any(tt==1)){
		cat("\tunique transition classes:", paste(sort(as.integer(names(tt[tt==1]))), collapse=", "), "\n")
	}
	if(any(tt>1)){
		cat("\tshared transition classes:", paste(sort(as.integer(names(tt[tt>1]))), collapse=", "),"\n")
	}
	cat("\n")
	class(x)="matrix"
	print(x)
}

print.rbm=function (x, printlen = 3, ...)
{
    cat("likelihood function for relaxed-rates univariate continuous trait evolution\n")
	aa=names(argn(x))
    cat("\targument names:", paste(aa, collapse = ", "))
	cat("\n\n")
	f=x
	attributes(f)=NULL
	cat("definition:\n")
	print(f)
	
	cat("\n\nNOTE: 'rates' vector should be given in order given by the function 'argn()'")
}


print.gprior=function(x, printlen = 3, ...){
    y=x
    attributes(y)=NULL
    print(y)
}

print.auteurRAW <- function(x, printlen=3, ...){
	cat("raw reversible-jump Markov Chain Monte Carlo (MCMC) output:\n")
	
	cat("\t", paste(names(x), collapse=", "))
}

print.rjmcmc <- function (x, printlen=3, ...)
{
    cat("reversible-jump Markov Chain Monte Carlo (MCMC) output:\n\tstart =", start(x),
    "\n\tend =", end(x), "\n\tthinning interval =", thin(x), "\n")
	
	nc=ncol(x)
    cat("\ncolnames:\n")
    if (nc > 2*printlen) {
        cat(paste(paste(colnames(x)[1:printlen], collapse = ", "), ", ..., ", paste(colnames(x)[(nc-printlen+1):nc], collapse = ", "), "\n", sep = ""))
    } else {
		cat(paste(colnames(x), collapse=", "))
	}
}

print.rjmcmcmc <- function (x, printlen=3, ...)
{
    cat("combined reversible-jump Markov Chain Monte Carlo (MCMC) output:\n\tstart =", start(x),
    "\n\tend =", end(x), "\n\tthinning interval =", thin(x), "\n\tchains =", attr(x,"nrun"), "\n")
	
	nc=ncol(x)
    cat("\ncolnames:\n")
    if (nc > 2*printlen) {
        cat(paste(paste(colnames(x)[1:printlen], collapse = ", "), ", ..., ", paste(colnames(x)[(nc-printlen+1):nc], collapse = ", "), "\n", sep = ""))
    } else {
		cat(paste(colnames(x), collapse=", "))
	}
}

print.mcmc.list <- function(x, printlen=3, ...)
{
    cat("combined reversible-jump Markov Chain Monte Carlo (MCMC) output:\n\tstart =", start(x),
    "\n\tend =", end(x), "\n\tthinning interval =", thin(x), "\n\tchains =", length(x), "\n")
}

print.transformer=function (x, printlen = 3, ...) 
{
    cat("function for tree transformation\n")
    cat("\targument names:", paste(argn(x), collapse = ", "))		 
	cat("\n\n")
	f=x
	attributes(f)=NULL
	cat("definition:\n")
	print(f)
}

## GENERIC
print.glnL=function(x, ...){
    print(as.numeric(x))
}


## GENERIC
print.bm=function (x, printlen = 3, ...) 
{
    cat("likelihood function for univariate continuous trait evolution\n")
    cat("\targument names:", paste(argn(x), collapse = ", "))		 
	cat("\n\n")
	f=x
	attributes(f)=NULL
	cat("definition:\n")
	print(f)
}

## GENERIC
print.mkn=
function (x, printlen = 3, ...) 
{
    cat("likelihood function for univariate discrete trait evolution\n")
    cat("\targument names:", paste(argn(x), collapse = ", "))	
    if(!is.null(al<-attr(x, "levels"))) {
        fmt=.format.levels.print(length(al))
        cat("\n\n\tmapping\n\t\t", paste(sprintf(fmt, 1:length(al)), al, collapse="\n\t\t", sep=": "), sep="")
    }
	cat("\n\n")
	f=x
	attributes(f)=NULL
	cat("definition:\n")
	print(f)
}


print.gfit=function(x, format=c("default", "oldestyle"), ...){
    
    format=match.arg(format, c("default", "oldestyle"))
    
    lik=x$lik
    if("bm"%in%class(lik)) type="continuous" else if("mkn"%in%class(lik)) type="discrete" else stop("unresolvable class of likelihood function")
        
    if(format=="oldestyle"){
        
        if(type=="discrete"){
          res=.oldestyle.gfit.discrete(x)
          res=list(res)
          names(res)="Trait1"
          return(res)
        }
        
        if(type=="continuous"){
            res=.oldestyle.gfit.continuous(x)
            res=list(res)
            names(res)="Trait1"
            return(res)
        }
        
    } else {
        solnfreq=function(x, tol = .Machine$double.eps^0.5){
            ll=logLik(x)
            aa=abs(x$res[,"lnL"]-ll)<=tol
            sum(aa)/length(aa)
        }
        
        
        ## OVERVIEW
        cat(paste("GEIGER-fitted comparative model of", type, "data", sep=" "))

        lik=x$lik
        att=attributes(lik)
        
        ## Q MATRIX
        if(type=="discrete"){
            Q=.Qmatrix.from.gfit(x)
            rownames(Q)=paste("    ", rownames(Q), sep="")

            cat("\n fitted Q matrix:\n")
            print(Q, ...)
        
            model=att$transform
            if(model!="none"){
                pars=att$argn[!att$trns]
            } else {
                pars=NULL
            }
        }
        
        if(type=="continuous"){
            model=att$model
            pars=c(argn(lik),"z0")
        }
        
        ## MODEL PARAMETERS
        if(model!="none"){
            cat("\n fitted ", sQuote(model), " model ", ifelse(length(pars)==1, "parameter:\n\t", "parameters:\n\t"), paste(names(x$opt[pars]), sprintf("%f", x$opt[pars]), sep=" = ", collapse="\n\t"), "\n", sep="")
        }
        
        ## SUMMARY
        rr=x$opt[c("lnL", "aic", "aicc")]
        names(rr)=c("log-likelihood", "AIC", "AICc")

        cat("\n model summary:\n\t", paste(names(rr), sprintf("%f",rr), sep=" = ", collapse="\n\t"), "\n\tfree parameters = ", x$opt$k, "\n", sep="")
        
        
        ## CONVERGENCE
        cat("\nConvergence diagnostics:\n\toptimization iterations = ",nrow(x$res),"\n\tfailed iterations = ", sum(rownames(x$res)=="FAIL"), "\n\tfrequency of best fit = ", sprintf("%.2f", solnfreq(x)), "\n", sep="")
        
        
        ## OBJECT SUMMARY
        cat("\n object summary:\n\t'lik' -- likelihood function\n\t'bnd' -- bounds for likelihood search\n\t'res' -- optimization iteration summary\n\t'opt' -- maximum likelihood parameter estimates\n")
    }
}

print.gfits=function(x, format=c("default", "oldestyle"), ...){
     format=match.arg(format, c("default", "oldestyle"))
    
    if(format=="oldestyle"){
    	lik=x[[1]]$lik
        if("bm"%in%class(lik)) type="continuous" else if("mkn"%in%class(lik)) type="discrete" else stop("unresolvable class of likelihood function")
        
        if(type=="discrete"){
          res=lapply(x, .oldestyle.gfit.discrete)
          names(res)=names(x)
          return(res)
         }
        
        if(type=="continuous"){
          res=lapply(x, .oldestyle.gfit.continuous)
          names(res)=names(x)
          return(res)
        }
    } else {
       print.default(x)
    }
}

.oldestyle.gfit.continuous=function(x){
    lik=x$lik
    arg=argn(lik)
    model=attr(lik, "model")
    oddoldemodels=c("kappa", "white")
    if(model%in%"white"){
        par=structure(c("z0", "sigsq"), names=c("mean", "nv"))      
    } else {
        par=structure(c("sigsq","lambda","kappa","delta","alpha","a","z0","drift"), names=c("beta","lambda","kappa","delta","alpha","a","mean","mu"))
    }
    
    zpar=x$opt[arg]
    names(zpar)=ifelse(names(zpar)%in%par, names(par), names(zpar))
    
    lnl=logLik(x)
    rr=x$res
    conv=as.numeric(rr[,"convergence"][min(which(rr[,"lnL"]==lnl))])

    res=c(lnl=lnl, zpar, convergence=conv, x$opt[c("aic", "aicc", "k")]) 
    return(res)
}

.Qmatrix.from.gfit=function(x){
    lik=x$lik
    
    numberize=function(x){
        y=gsub("q", "", x)
        sp=(nn<-nchar(y))/2
        as.numeric(c(substring(y, 1, sp), substring(y, sp+1, nn)))
    }
        
    att=attributes(lik)
    att$k=length(att$levels)
    Qmat=matrix(0, att$k, att$k)

    nms=att$argn[att$trns]
    other=att$argn[!att$trns]

    if("constrained"%in%class(lik)){
        cpars=x$opt[argn(lik)]
        apars=names(lik(unlist(cpars), pars.only=TRUE))
        nms=apars[!apars%in%other]
    }
    trns=x$opt[nms]
    for(i in 1:length(trns)){
        nm=names(trns)[i]
        idx=numberize(nm)
        Qmat[idx[1], idx[2]]=trns[[i]]
    }
    diag(Qmat)=-rowSums(Qmat)
    rownames(Qmat)<-colnames(Qmat)<-levels(lik)
    Qmat   
}


.oldestyle.gfit.discrete=function(x){

           lik=x$lik
           numberize=function(x){
                y=gsub("q", "", x)
                sp=(nn<-nchar(y))/2
                as.numeric(c(substring(y, 1, sp), substring(y, sp+1, nn)))
            }
        
            att=attributes(lik)
            att$k=length(att$levels)

            nms=att$argn[att$trns]
            other=att$argn[!att$trns]

            if(length(other)) tt=x$opt[[other]] else tt=NULL
            
            rr=x$res
            lnl=logLik(x)
            conv=as.numeric(rr[,"convergence"][min(which(rr[,"lnL"]==lnl))])
            msg=ifelse(conv==0, "R thinks that this is the right answer.", "Warning: may not have converged to a proper solution.")
            res=list(lnl=lnl, q=.Qmatrix.from.gfit(x), treeParam=tt, message=msg)
            res=res[!sapply(res, is.null)]
            return(res)
}