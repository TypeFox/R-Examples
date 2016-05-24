# param is a named numeric vector or a matrix.
# if param is a vector, it specifies a 5PL curve. 
# if param is a matrix, each row of it specifies a 5PL curve
# param can be in either the classical parameterization or the g-h parameterization
# param can come from drm fit or bcrm fit
# the ed50 parameterization is not currently supported
# return a list with 5 elements: b, c, d, e, f
get.curve.param.list=function(param){
    
    if (is.matrix(param)) {
        # do nothing
    } else if (is.vector(param)) {
        param=matrix(param, nrow=1, dimnames=list(NULL, names(param))) # as.matrix turns a vector into a column
    } else {
        stop ("param is not matrix or vector")
    }
    
    tmp=substr(colnames(param),1,1)
    tmp["logtao"==colnames(param)]="logtao" # logtao can not be shortened to l
    colnames(param)=tmp

    if (!"c" %in% colnames(param)) stop("param does not have d")
    if (!"d" %in% colnames(param)) stop("param does not have d")
    
    if("b" %in% tmp & "e" %in% tmp) { 
        # classical parameterization
    } else if ("g" %in% tmp & "h" %in% tmp) {
        # gh parameterization
        param=gh2cla(param) 
    } else if ("logtao" %in% tmp & "b" %in% tmp) {
        # ED50 parameterization
        param=ed502cla(param) 
    } else if ("logtao" %in% tmp & ("h" %in% tmp | "logh" %in% tmp)) {
        # ED50 parameterization
        param=ed50b2cla(param) 
    } else {
        stop("not gh not cla not ED50: "%+%concatList(colnames(param),","))
    }
    
    b=param[,"b"]; c=param[,"c"]; d=param[,"d"]; e=param[,"e"]
    names(b)=NULL; names(c)=NULL; names(d)=NULL; names(e)=NULL; # so that the outcome is not named incorrectly
    if (any(e<0)) stop("e cannot be negative")
    res=list(b=b, c=c, d=d, e=e)
    
    if ("f" %in% colnames(param)) {
        f=param[,"f"]    
        names(f)=NULL
        if (any(f<0)) stop("f cannot be negative")
        res=c(res, f=list(f))
    } else {
        res=c(res, f=list(rep(1, length(e))))
    }
        
    res
}


############################################################################
# classical parameterization


# following functions are vectorized for t, but not for param
FivePL.t=function (t,param) {
    param.list=get.curve.param.list(param)
    if (!(length(t)==1 | length(param.list$c)==1 | length(t)==length(param.list$c))) stop ("t or x and param not match in length")
    
    out=with(param.list, (d-c)/{1+exp(b*t-b*log(e))}^f+c)
    names(out)=rep("y", length(out))
    out
}

FivePL.t.func=function (param) {
    param.list=get.curve.param.list(param)
    return (with(param.list, function (t) (d-c)/{1+exp(b*t-b*log(e))}^f+c))
}

FivePL.x=function (x,param) {
    FivePL.t(log(x), param)
}


# returns the concentration for a given set of y, if out of bound, returns either 0 or Inf
FivePL.t.inv = function (y,param) {
    
    param.list=get.curve.param.list(param)
    if (!(length(y)==1 | length(param.list$c)==1 | length(y)==length(param.list$c))) stop ("y and param not match in length")
    
    with(param.list, {
        out=suppressWarnings(1/b*log( ((d-c)/(y-c))^(1/f) - 1 ) + log(e))
        # deal with NaN, need to do this b/c of the need to differentiate between out of left bound and out of right bound, 
        # which is needed when a replacement value will be used
        out = ifelse(y>c & y<d, out, {
            ifelse(y<c, -Inf, Inf)*ifelse(b<0,1,-1)
        })
        names(out)=rep("x", length(out))
        out
    })
    
}

FivePL.t.inv.func=function (param) {
    
    param.list=get.curve.param.list(param)
    
    with(param.list, {
        function (y) {            
            out=suppressWarnings(1/b*log( ((d-c)/(y-c))^(1/f) - 1 ) + log(e))
            # deal with NaN, need to do this b/c of the need to differentiate between out of left bound and out of right bound, 
            # which is needed when a replacement value will be used
            out = ifelse(y>c & y<d, out, {
                ifelse(y<c, -Inf, Inf)*ifelse(b<0,1,-1)
            })
            names(out)=rep("x", length(out))
            out
        }
    })
    
}


FivePL.x.inv = function (y,param) exp(FivePL.t.inv(y,param))
    
FivePL.x.inv.func=function (param) {
    function (y) {            
        exp(FivePL.t.inv.func(param)(y))
    }    
}




FourPL.x <- function (x, param) {
    # convert param to a matrix if necessary
    if (is.matrix(param)) {
        # do nothing
    } else if (is.vector(param)) {
        param=matrix(param, nrow=1, dimnames=list(NULL, substr(names(param),1,1))) # as.matrix turns a vector into a column
    } else {
        stop ("param is not matrix or vector")
    }    
    # add f=1
    param=cbind(param, f=1)
    
    FivePL.x(x,param)        
}

FourPL.x.inv <- function (y, param) {
    # convert param to a matrix if necessary
    if (is.matrix(param)) {
        # do nothing
    } else if (is.vector(param)) {
        param=matrix(param, nrow=1, dimnames=list(NULL, substr(names(param),1,1))) # as.matrix turns a vector into a column
    } else {
        stop ("param is not matrix or vector")
    }    
    # add f=1
    param=cbind(param, f=1)
    
    FivePL.x.inv(y,param)
}

FourPL.t.func=function (param) {
    # convert param to a matrix if necessary
    if (is.matrix(param)) {
        # do nothing
    } else if (is.vector(param)) {
        param=matrix(param, nrow=1, dimnames=list(NULL, substr(names(param),1,1))) # as.matrix turns a vector into a column
    } else {
        stop ("param is not matrix or vector")
    }    
    # add f=1
    param=cbind(param, f=1)
    FivePL.t.func(param)
}


ED5PL = function (param, tao) {
    
    names(param)=substr(names(param),1,1)
    b=param["b"]; c=param["c"]; d=param["d"]; e=param["e"]; f=param["f"]; g=param["g"]; h=param["h"]; logtao=param["logtao"]
    if(is.na(b) | is.na(e)) {
        if(!is.na(g) & !is.na(h)) {
            pp=gh2cla(c(c,d,g,h=h,f=f)) 
            b=pp["b"]
            e=pp["e"]
        } else if (!is.na(logtao)) {
            pp=ed502cla(c(c,d,b=b,logtao=unname(logtao),f=f))            
            e=pp["e"]
        } else  {
            stop("parameterization not recognized")
        }
    } 
    
    unname(e*(tao^{-1/f}-1)^{1/b})
}


#Treat out of bound concentration estimates
#y: a number. The readout
#p: a vector of number. Parameters for a 5pl/4pl curve.
#t.range: a vector of two numbers. The range of log standard samples concentrations.
#If y is less than lower asymptote, return t.range[1]+log(1/2), i.e. log of half of smallest standard samples concentration.
#If y is higher than upper asymptote, return t.range[2], i.e. log of largest standard samples concentration
treat.out.of.bound=function(y, p, t.range){
    if (y<get.curve.param.list(p)$c) {
        t.0=t.range[1]+log(1/2) # half of the smallest standard concentration
    } else if (y>get.curve.param.list(p)$d) {
        t.0=t.range[2] # the largest standard concentration
    } else stop("treat.out.of.bound: this cannot be right.")
    t.0
}




###########################################################################################################
# simulate one curve, return FI, a vectorized function
simulate1curve=function(param, t, sd.e=0.1, expy=TRUE, gamma=0) {
    if (expy) {
        .mean=FivePL.t(t, param)
        y = rnorm (n=length(.mean), mean=.mean, sd=sd.e*.mean^(gamma/2))
        exp(y) 
    } else {
        .mean=(FivePL.t(t, param)) 
        y = rnorm (n=length(.mean), mean=.mean, sd=sd.e*.mean^(gamma/2))
        y
    }
}



###########################################################################################################
# summary measure for comparing two curves

get.abs.dev = function(p1, p2, t.range, y.range) {
    if (!"f" %in% names(p1)) p1=c(p1, f=1)
    if (!"f" %in% names(p2)) p2=c(p2, f=1)

    f.0  =FivePL.t.inv.func(p1)
    f.hat=FivePL.t.inv.func(p2)

    integrate( function(y) {
        t.0 =   f.0  (y)
        t.0 = sapply (1:length(y), function (i) { # as y is a vector and treat.out.of.bound is not vectorized, we need to loop through the lenght of y
            if (is.nan(t.0[i])) treat.out.of.bound (y[i], p1, t.range) else t.0[i]
        })
        t.hat = f.hat(y)
        t.hat = sapply (1:length(y), function (i) {
            if (is.nan(t.hat[i]) | Inf==abs(t.hat[i])) treat.out.of.bound (y[i], p2, t.range) else t.hat[i]
        })
        abs(t.hat - t.0)
    }, lower=y.range[1], upper=y.range[2], subdivisions=1000 )$value/(y.range[2]-y.range[1])
}

# get area between two curves betwee two curves
get.abc = function(p1, p2, t.range) {
    
    if (!"f" %in% names(p1)) p1=c(p1, f=1)
    if (!"f" %in% names(p2)) p2=c(p2, f=1)
    
    f1=FivePL.t.func( p1 )
    f0=FivePL.t.func( p2 )
    
    integrate( function(t) abs(f1(t)-f0(t)) , lower=t.range[1], upper=t.range[2], subdivisions=1000 )$value/(t.range[2]-t.range[1])
    
}

# get S1 betwee two curves
get.S1 = function(p1, p2, t.range) {
    if (!"f" %in% names(p1)) p1=c(p1, f=1)
    if (!"f" %in% names(p2)) p2=c(p2, f=1)
    
    f1=FivePL.t.func( p1 )
    f0=FivePL.t.func( p2 )
    
    integrate( function(t) (f1(t)-f0(t))^2 , lower=t.range[1], upper=t.range[2], subdivisions=1000 )$value/(t.range[2]-t.range[1])
}

# get S2 betwee two curves, percent bias
get.S2 = function(p1, p2, t.range) {
    if (!"f" %in% names(p1)) p1=c(p1, f=1)
    if (!"f" %in% names(p2)) p2=c(p2, f=1)
    
    f0=FivePL.t.func( p1 )
    f1.inv=FivePL.x.inv.func(p2)
    
    integrate( function(t) {
        x0= exp(t)
        x1= f1.inv(f0(t))
        abs(x1-x0)/x0 * 100
    }, lower=t.range[1], upper=t.range[2], subdivisions=1000 )$value/(t.range[2]-t.range[1])
}
