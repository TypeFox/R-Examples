#'@describeIn ztable
#'
ztable.matrix=function(x,digits=NULL,...){
    result=data.frame(x,stringsAsFactors=FALSE)
    colnames(result)=colnames(x)
    out=ztable(result,...)
    out
}

#'@describeIn ztable
#'
ztable.lm=function(x,digits=NULL,...){
    result=data.frame(summary(x)$coeff)
    colnames(result)=c("Estimate","Std. Error","t value","Pr(>|t|)")
    h=deparse(x$call)
    #h=gsub("~","$\\sim$",h,fixed=TRUE)
    h=paste("Call: ",h,sep="")
    attr(result,"footer")=h
    if (is.null(digits)) mydigits=c(1,4,4,2,4)
    else mydigits=digits
    out=ztable_sub(result,digits=mydigits,...)
    out
}

#'@describeIn ztable
#'
ztable.fitdistr=function(x,digits=NULL,...){
    if(is.null(digits)) mydigits=3
    else mydigits=digits

    result=rbind(x$estimate,x$sd)
    rownames(result)=c("estimate","sd")
    result=data.frame(result)

    string=paste("N=",x$n,", The log-likelihood=",round(x$loglik,2),sep="")
    attr(result,"footer")=string
    out=ztable_sub(result,digits=mydigits)
    out
}

#'@describeIn ztable
#'
ztable.nls=function(x,digits=NULL,...){
    result=data.frame(summary(x)$coeff)
    colnames(result)=c("Estimate","Std. Error","t value","Pr(>|t|)")

    s=deparse(formula(x))
    h1=paste("  model: ", s,"\n",sep="")
    h2=paste("  data: ", deparse(x$data),"\n", sep = "")
    h=c("Nonlinear regression model\n",h1,h2)
    attr(result,"heading")=h
    if (is.null(digits)) mydigits=c(1,4,4,2,4)
    else mydigits=digits
    out=ztable_sub(result,digits=mydigits,...)
    out
}


#'@describeIn ztable
#'
ztable.aov=function(x,digits=NULL,...){
    result=summary(x)[[1]]
    if(!is.null(x$call)){
        h=deparse(x$call)
        h=paste("Call: ",h,sep="")
        attr(result,"footer")=h
    }
    if (is.null(digits)) mydigits=c(1,0,2,2,2,4)
    else mydigits=digits
    out=ztable_sub(result,digits=mydigits,...)
    out
}

#'@describeIn ztable
#'
ztable.anova=function(x,digits=NULL,...){
    result=data.frame(x)
    colnames(result)=colnames(x)
    if(is.null(digits)) {
        if(ncol(x)==4) mydigits=c(1,0,2,0,2)
        else if (ncol(x)==5) mydigits=c(1,0,2,2,2,4)
        else mydigits=c(1,1,2,1,2,2,4)
    }
    else mydigits=digits
    #attr(result,"heading")=attr(x,"heading")
    h=c()
    if(!is.null(attr(x,"heading"))) {
        heading=attr(x,"heading")
        for(i in 1:length(heading)) {
            h=c(h,unlist(strsplit(heading[i],"\n")))
        }
    }
    attr(result,"heading")=h
    if(!is.null(x$call)){
        h=deparse(x$call)
        h=paste("Call: ",h,sep="")
        attr(result,"footer")=h
    }
    out=ztable_sub(result,digits=mydigits,...)
    out
}

#'@describeIn ztable
#'
ztable.glm=function(x,digits=NULL,...){
    a=summary(x)$coeff
    b=data.frame(a)
    colnames(b)=colnames(a)

    suppressMessages(d<-confint(x))
    OR=data.frame(exp(coef(x)),exp(d))
    OR=round(OR,4)
    OR=cbind(OR,round(summary(x)$coefficient[,4],4))
    #result=na.omit(result)
    colnames(OR)=c("OR","lcl","ucl","p")
    i=apply(OR,1,function(x) any(is.na(x)))
    OR[i,c(1,2,3)]=NA
    out=cbind(b,OR[c(1,2,3)])
    h=deparse(x$call)
    if(length(h)==1) h=paste("Call: ",h,sep="")
    else if(length(h)==2) h=paste("Call: ",h[1],h[2],sep="")
    attr(out,"footer")=h

    if (is.null(digits)) mydigits=c(1,4,4,2,4,2,2,2)
    else mydigits=digits
    out=ztable_sub(out,digits=mydigits,...)
    out
}

#'@describeIn ztable
#'
ztable.coxph=function(x,digits=NULL,...){
    a=summary(x)$coeff
    b=summary(x)$conf.int
    result=cbind(b[,c(1,3,4)],a[,c(3,4,5)])
    result=data.frame(result)
    h=deparse(x$call)
    if(length(h)==1) h=paste("Call: ",h,sep="")
    else if(length(h)==2) h=paste("Call: ",h[1],h[2],sep="")

    attr(result,"footer")=h
    colnames(result)=c("HR","lcl", "ucl", "se(coef)","z","Pr(>|z|)")
    if (is.null(digits)) mydigits=c(0,3,3,3,3,3,4)
    else mydigits=digits
    out=ztable_sub(result,digits=mydigits,...)
    out
}

#'@describeIn ztable
#'
ztable.prcomp=function(x,digits=NULL,...){
    result=data.frame(x$rotation)
    colnames(result)=colnames(x$rotation)
    attr(result,"heading") <- "Rotation:"
    if(!is.null(x$call)){
        h=deparse(x$call)
        h=paste("Call: ",h,sep="")
        attr(result,"footer")=h
    }
    if (is.null(digits)) mydigits=c(1,4,4,4,4)
    else mydigits=digits
    out=ztable_sub(result,digits=mydigits,...)
    out
}

#'@describeIn ztable
#'
ztable.summary.prcomp=function(x,digits=NULL,...){
    result=data.frame(x$importance)
    colnames(result)=colnames(x$importance)
    attr(result,"heading") <- "Importance of components:"
    if(!is.null(x$call)){
        h=deparse(x$call)
        h=paste("Call: ",h,sep="")
        attr(result,"footer")=h
    }
    if (is.null(digits)) mydigits=c(1,4,4,4,4)
    else mydigits=digits
    out=ztable_sub(result,digits=mydigits,...)
    out
}

