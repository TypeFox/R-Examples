############################################################################
# parameterization 3: b,c,d,e,f -> b,c,d,g,h, where g is the inflection point, and h is the hill slope at the inflection point



# vectorized functions for converting parameters between parameterizations
# param can be a vector or a matrix where each row is a parameter value
cla2gh=function(param){
    is.v=FALSE
    if(is.vector(param)) {
        is.v=TRUE
        tmp=substr(names(param),1,1)
        param=matrix(param, nrow=1)
        colnames(param)=tmp
    } else {
        colnames(param)=substr(colnames(param),1,1)
    }
    
    if(ncol(param)==4) {
        is.4pl=TRUE
        param=cbind(param, "f"=rep(1,nrow(param)))
    } else {
        is.4pl=FALSE
    }
    
    if(!"b" %in% colnames(param) & "logmb"%in%colnames(param)) b=unname(-exp(param[,"logmb"])) else b=param[,"b"]
    if(!"e" %in% colnames(param) & "loge"%in%colnames(param)) e=unname(exp(param[,"loge"])) else e=param[,"e"]
    if(!"f" %in% colnames(param) & "logf"%in%colnames(param)) f=unname(exp(param[,"logf"])) else f=param[,"f"]    
    c=param[,"c"]; d=param[,"d"]; 
    g=log(e)-(1/b)*log(f)
    h=-b*(d-c)/(1+(1/f))^(f+1)
    if(is.v) {
        res=c(c,d,g=unname(g),h=unname(h))
        if (!is.4pl) res=c(res, f) 
    } else {
        res=cbind(c,d,g=unname(g),h=unname(h))
        if (!is.4pl) res=cbind(res, f) 
    }
    res
}

gh2cla=function(param){
    is.v=FALSE
    if(is.vector(param)) {
        is.v=TRUE
        tmp=substr(names(param),1,1)
        param=matrix(param, nrow=1)
        colnames(param)=tmp
    }
    
    if(ncol(param)==4) {
        is.4pl=TRUE
        param=cbind(param, "f"=rep(1,nrow(param)))
    } else {
        is.4pl=FALSE
    }
    
    c=param[,"c"]; d=param[,"d"]; g=param[,"g"]
    if(!"f" %in% colnames(param) & "logf"%in%colnames(param)) f=unname(exp(param[,"logf"])) else f=param[,"f"]    
    if(!"h" %in% colnames(param) & "logh"%in%colnames(param)) h=unname(exp(param[,"logh"])) else h=param[,"h"]    
    
    b=-(h/(d-c))*(1+(1/f))^{f+1}
    e=exp(g+(1/b)*log(f))
    if(is.v) {
        res=c(b=unname(b),c,d,e=unname(e))
        if (!is.4pl) res=c(res, f)
    } else {
        res=cbind(b=unname(b),c,d,e=unname(e))
        if (!is.4pl) res=cbind(res, f)
    }
    res
}

cla2ed50=function(param){
    is.v=FALSE
    if(is.vector(param)) {
        is.v=TRUE
        tmp=substr(names(param),1,1)
        param=matrix(param, nrow=1)
        colnames(param)=tmp
    }else colnames(param)=substr(colnames(param),1,1)
    
    if(ncol(param)==4) {
        is.4pl=TRUE
        param=cbind(param, "f"=rep(1,nrow(param)))
    } else {
        is.4pl=FALSE
    }
    
    b=param[,"b"]; c=param[,"c"]; d=param[,"d"]; e=param[,"e"]; f=param[,"f"]    
    tao=e*(2^{1/f}-1)^{1/b}
    if(is.v) {
        res=c(b,c,d,logtao=unname(log(tao)))
        if (!is.4pl) res=c(res, f) 
    } else {
        res=cbind(b,c,d,logtao=unname(log(tao)))
        if (!is.4pl) res=cbind(res, f) 
    }
    res
}

cla2ed50b=function(param){
    is.v=FALSE
    if(is.vector(param)) {
        is.v=TRUE
        tmp=substr(names(param),1,1)
        param=matrix(param, nrow=1)
        colnames(param)=tmp
    }else colnames(param)=substr(colnames(param),1,1)
    
    if(ncol(param)==4) {
        is.4pl=TRUE
        param=cbind(param, "f"=rep(1,nrow(param)))
    } else {
        is.4pl=FALSE
    }
    
    b=param[,"b"]; c=param[,"c"]; d=param[,"d"]; e=param[,"e"]; f=param[,"f"]    
    tao=e*(2^{1/f}-1)^{1/b}
    h=-b*(d-c)/(1+(1/f))^(f+1)
    names(h)="h"
    if(is.v) {
        res=c(c,d,logtao=unname(log(tao)),h)
        if (!is.4pl) res=c(res, f) 
    } else {
        res=cbind(c,d,logtao=unname(log(tao)),h)
        if (!is.4pl) res=cbind(res, f) 
    }
    res
}

ed502cla=function(param){
    is.v=FALSE
    if(is.vector(param)) {
        is.v=TRUE
#        tmp=substr(names(param),1,1)
        tmp=names(param)
        param=matrix(param, nrow=1)
        colnames(param)=tmp
    }
    
    if(ncol(param)==4) param=cbind(param, "f"=rep(1,nrow(param)))
    
    c=param[,"c"]; d=param[,"d"]; b=param[,"b"]; f=param[,"f"]; logtao=param[,"logtao"]
    if(!"f" %in% colnames(param) & "logf"%in%colnames(param)) f=unname(exp(param[,"logf"])) else f=param[,"f"]    
    if(!"b" %in% colnames(param) & "logmb"%in%colnames(param)) b=unname(-exp(param[,"logmb"])) else b=param[,"b"]
    
    loge=logtao-(1/b)*log(2^{1/f}-1)
    e=exp(loge)
    if(is.v) c(b=unname(b),c,d,e=unname(e),f) else cbind(b=unname(b),c,d,e=unname(e),f)
}




ed50b2cla=function(param){
    is.v=FALSE
    if(is.vector(param)) {
        is.v=TRUE
#        tmp=substr(names(param),1,1)
        tmp=names(param)
        param=matrix(param, nrow=1)
        colnames(param)=tmp
    }
    
    if(ncol(param)==4) param=cbind(param, "f"=rep(1,nrow(param)))
    
    c=param[,"c"]; d=param[,"d"]; h=param[,"h"]; f=param[,"f"]; logtao=param[,"logtao"]
    if(!"f" %in% colnames(param) & "logf"%in%colnames(param)) f=unname(exp(param[,"logf"])) else f=param[,"f"]    
    if(!"h" %in% colnames(param) & "logh"%in%colnames(param)) h=unname(exp(param[,"logh"])) else h=param[,"h"]    
    
    b=-h/((d-c)/(1+(1/f))^(f+1))
    loge=logtao-(1/b)*log(2^{1/f}-1)
    e=exp(loge)
    if(is.v) c(b=unname(b),c,d,e=unname(e),f) else cbind(b=unname(b),c,d,e=unname(e),f)
}




# x can be a single number or a vector
# return a matrix, where each row corresponds to an x
vpl1.deriv = function(x,param){
    b=param["b"]; c=param["c"]; d=param["d"]; e=param["e"]; f=param["f"]    
    u=(x/e)^{b}
    y=c+((d-c)/((1+u)^{f}))
    res=matrix(
        c(  "c" = 1-(1/((1+u)^{f})), 
            "d" = (1/((1+u)^{f})), 
            "e" = (y-c)*((b*f)/e)*(u/(1+u)),
            "loge" = (y-c)*(b*f)*(u/(1+u)),
            "f" = -(y-c)*log(1+u),
            "b" = -(y-c)*(f/b)*(u/(1+u))*log(u)
        ), 
        nrow=length(x)
    )
    colnames(res)=c("c","d","e","loge","f","b")
    res
}
## test
#vpl1.deriv(c(1,100), coef(fit))

# return a list of derivative functions
vpl1.deriv.func = function(param){
    b=param["b"]; c=param["c"]; d=param["d"]; e=param["e"]; f=param["f"]    
    f.c=function(x) {
        u=(x/e)^{b}; y=c+((d-c)/((1+u)^{f}))
        unname(1-(1/((1+u)^{f})))
    }
    f.d=function(x) {
        u=(x/e)^{b}; y=c+((d-c)/((1+u)^{f}))
        unname((1/((1+u)^{f})))
    }
    f.e=function(x) {
        u=(x/e)^{b}; y=c+((d-c)/((1+u)^{f}))
        unname((y-c)*((b*f)/e)*(u/(1+u)))
    }
    f.loge=function(x) {
        u=(x/e)^{b}; y=c+((d-c)/((1+u)^{f}))
        unname((y-c)*b*f*(u/(1+u)))
    }
    f.f=function(x) {
        u=(x/e)^{b}; y=c+((d-c)/((1+u)^{f}))
        unname(-(y-c)*log(1+u))
    }
    f.b=function(x) {
        u=(x/e)^{b}; y=c+((d-c)/((1+u)^{f}))
        unname(-(y-c)*(f/b)*(u/(1+u))*log(u))
    }
    list("c"=f.c, "d"=f.d, "e"=f.e, "loge"=f.loge, "f"=f.f, "b"=f.b)
}


############################################################################
# parameterization 2: b,c,d,e,f -> b,c,d,tao,f, where tao is the ED50

# x can be a single number or a vector
# return a matrix, where each row corresponds to an x
vpl2.deriv = function(x,param){
    c=param["c"]; d=param["d"]; logtao=param["logtao"]; b=param["b"]; f=param["f"]    
    t=log(x); u=exp(b*(t-logtao)+log(2^{1/f}-1)); y=c+((d-c)/((1+u)^{f}))
    matrix(
        c(  "c" = 1-(1/((1+u)^{f})), 
            "d" = (1/((1+u)^{f})), 
            "logtao" = b*f*(y-c)*(u/(1+u)),
            "b" = -f*(t-logtao)*(y-c)*(u/(1+u)),
            "f" = ((log(2))/f)*((2^{1/f})/(2^{1/f}-1))*(y-c)*(u/(1+u))
        ), 
        nrow=length(x)
    )
}
## test
#vpl2.deriv(c(1,100), theta)

# return a list of derivative functions
vpl2.deriv.func = function(param){
    c=param["c"]; d=param["d"]; logtao=param["logtao"]; b=param["b"]; f=param["f"]    
    f.c=function(x) {
        t=log(x); u=exp(b*(t-logtao)+log(2^{1/f}-1)); y=c+((d-c)/((1+u)^{f}))
        unname(1-(1/((1+u)^{f})))
    }
    f.d=function(x) {
        t=log(x); u=exp(b*(t-logtao)+log(2^{1/f}-1)); y=c+((d-c)/((1+u)^{f}))
        unname((1/((1+u)^{f})))
    }
    f.logtao=function(x) {
        t=log(x); u=exp(b*(t-logtao)+log(2^{1/f}-1)); y=c+((d-c)/((1+u)^{f}))
        unname(b*f*(y-c)*(u/(1+u)))
    }
    f.b=function(x) {
        t=log(x); u=exp(b*(t-logtao)+log(2^{1/f}-1)); y=c+((d-c)/((1+u)^{f}))
        unname(-f*(t-logtao)*(y-c)*(u/(1+u)))
    }
    f.f=function(x) {
        t=log(x); u=exp(b*(t-logtao)+log(2^{1/f}-1)); y=c+((d-c)/((1+u)^{f}))
        unname(((log(2))/f)*((2^{1/f})/(2^{1/f}-1))*(y-c)*(u/(1+u)))
    }
    list("c"=f.c, "d"=f.d, "logtao"=f.logtao, "b"=f.b, "f"=f.f)
}



# x can be a single number or a vector
# return a matrix, where each row corresponds to an x
vpl3.deriv = function(x,param){
    c=param["c"]; d=param["d"]; g=param["g"]; h=param["h"]; f=param["f"]    
    t=log(x); u=(1/f)*exp(-(h/(d-c))*(1+(1/f))^{f+1}*(t-g)); y=c+((d-c)/((1+u)^{f}))
    matrix(
        c(  "c" = 1-(1/((1+u)^{f})), 
            "d" = (1/((1+u)^{f})), 
            "g" = -h*((y-c)/(d-c))*f*(1+(1/f))^{f+1}*(u/(1+u)),
            "h" = (t-g)*((y-c)/(d-c))*f*(1+(1/f))^{f+1}*(u/(1+u)),
            "f" = -(y-c)*{log(1+u)-(u/(1+u))*(log(u)+log(f)+1)+f*(u/(1+u))*(log(u)+log(f))*log(1+(1/f))}
        ), 
        nrow=length(x)
    )
}
## test
#vpl3.deriv(c(1,100), theta)

# return a list of derivative functions
vpl3.deriv.func = function(param){
    c=param["c"]; d=param["d"]; g=param["g"]; h=param["h"]; f=param["f"]    
    f.c=function(x) {
        t=log(x); u=(1/f)*exp(-(h/(d-c))*(1+(1/f))^{f+1}*(t-g)); y=c+((d-c)/((1+u)^{f}))
        unname(1-(1/((1+u)^{f})))
    }
    f.d=function(x) {
        t=log(x); u=(1/f)*exp(-(h/(d-c))*(1+(1/f))^{f+1}*(t-g)); y=c+((d-c)/((1+u)^{f}))
        unname((1/((1+u)^{f})))
    }
    f.g=function(x) {
        t=log(x); u=(1/f)*exp(-(h/(d-c))*(1+(1/f))^{f+1}*(t-g)); y=c+((d-c)/((1+u)^{f}))
        unname(-h*((y-c)/(d-c))*f*(1+(1/f))^{f+1}*(u/(1+u)))
    }
    f.h=function(x) {
        t=log(x); u=(1/f)*exp(-(h/(d-c))*(1+(1/f))^{f+1}*(t-g)); y=c+((d-c)/((1+u)^{f}))
        unname((t-g)*((y-c)/(d-c))*f*(1+(1/f))^{f+1}*(u/(1+u)))
    }
    f.f=function(x) {
        t=log(x); u=(1/f)*exp(-(h/(d-c))*(1+(1/f))^{f+1}*(t-g)); y=c+((d-c)/((1+u)^{f}))
        unname(-(y-c)*{log(1+u)-(u/(1+u))*(log(u)+log(f)+1)+f*(u/(1+u))*(log(u)+log(f))*log(1+(1/f))})
    }
    list("c"=f.c, "d"=f.d, "g"=f.g, "h"=f.h, "f"=f.f)
}

# ylim=NULL; col=NULL; lty=NULL; lwd=1; plot.legend=FALSE; add=FALSE; legend=NULL; main=NULL # default 
lines5PL=function(param, xlim, ...) plot5PL(param, xlim, add=TRUE, ...)
plot5PL=function(param, xlim, ylim=NULL, col=NULL, lty=NULL, lwd=1, plot.legend=FALSE, add=FALSE, legend=NULL, main=NULL, xlab=NULL, ylab=NULL, xaxt="s", 
    yaxis.log.scale=FALSE, expy=FALSE, logy=FALSE) {
    
    if (!is.matrix(param) & !is.vector(param)) stop("param has to be vector or matrix")
    if(is.vector(param)) {
        tmp=names(param)
        param=matrix(param, nrow=1)
        dimnames(param)[[2]]=tmp
    }
    
    tlim=log(xlim)
    t.1=seq(tlim[1], tlim[2], length=1000)
    if (is.null(ylim)) ylim=range(apply(param,1,function(x) FivePL.t(t.1, x)))
    if (expy) ylim=exp(ylim)
    if (logy) ylim=log(ylim)
    if (!add) {
        plot(1,1,xlim=xlim, ylim=ylim, type="n", xlab=ifelse(is.null(xlab),"t",xlab), ylab=ifelse(is.null(ylab),"y",ylab), main=main, xaxt=xaxt, log=ifelse(yaxis.log.scale,"xy","x"))
#        plot(1,1,xlim=tlim, ylim=ylim, type="n", xlab=ifelse(is.null(xlab),"x",xlab), ylab=ifelse(is.null(ylab),"y",ylab), main=main, xaxt="n", log=log)
#        axis(side=1, at=seq(xlim[1], xlim[2], length=10), labels=round(exp(seq(xlim[1], xlim[2], length=10)),1))
    }
    if (is.null(col)) col=1:nrow(param) else if(length(col)==1) col=rep(col, nrow(param))
    if (is.null(lty)) lty=rep(1, nrow(param))
    for (i in 1:nrow(param)) {
        yy=FivePL.t(t.1, param[i,])
        if (expy) yy=exp(yy)
        if (logy) yy=log(yy)
        lines(exp(t.1), yy, col=col[i], lty=lty[i], lwd=lwd)
        #lines(t.1, yy, col=col[i], lty=lty[i], lwd=lwd)
    }
    
    if(is.null(legend)) legend=1:nrow(param)
    if (!add & plot.legend) mylegend(x=9, legend=legend, col=col, lty=lty)
}
