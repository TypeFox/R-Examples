svycdf<-function(formula,design,na.rm=TRUE,...) UseMethod("svycdf",design)

svycdf.default<-function(formula, design,na.rm=TRUE,...){
    if (inherits(formula, "formula")) 
        x <- model.frame(formula, model.frame(design), na.action = na.pass)
    else if (typeof(formula) %in% c("expression", "symbol")) 
        x <- eval(formula, model.frame(design, na.action = na.pass))
    else x<-formula
    if (na.rm) {
        nas <- rowSums(is.na(x))
        x <- x[nas == 0, , drop = FALSE]
    }
    rval<-vector("list",ncol(x))
    names(rval)<-names(x)
    for(i in 1:ncol(x)){
    		xx<-x[,i]
    		w <- weights(design,type="sampling")[nas==0]
    		oo<-order(xx)
    		cum.w<-cumsum(w[oo])/sum(w)
    		cdf <- approxfun( xx[oo],cum.w, method = "constant", 
         	   yleft =0, yright =1,ties="max")

    		class(cdf)<-"stepfun"
                call.i<-match.call()
                call.i$formula<-as.formula(paste("~",names(x)[i]))
    		attr(cdf,"call")<-call.i
    		rval[[names(x)[i]]]<-cdf
	}
    class(rval)<-"svycdf"
    cc<-sys.call()
    cc[[1]]<-as.name(.Generic)
    attr(rval,"call")<-cc
    rval
}


print.svycdf<-function(x,...){
	cat("Weighted ECDFs: ")
	print(attr(x,"call"))
	invisible(x)
	}

plot.svycdf<-function(x,xlab=NULL,...){
  if(is.null(xlab)) 
    xlab<-names(x)
  else if (length(xlab)==1)
    xlab<-rep(xlab,length(names(x)))
  
  for (i in 1:length(x)) plot(x[[i]], xlab =xlab[i], ...)
}
