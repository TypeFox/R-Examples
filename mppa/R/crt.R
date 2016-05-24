mkF <- function(x, start=0, end=1, adjust=1, disallow.zero=TRUE){
    f = density(x, adjust=adjust, from=start, to=end);
    if (disallow.zero) {m=min(f$y[f$y>1e-6]); f$y=pmax(f$y,m)}
    dx=(f$x)[2]-(f$x)[1]
    Fy=cumsum(f$y)*dx
    approxfun(f$x, Fy, yleft=min(Fy), yright=max(Fy))
}

htest <- function(p, T, data.name, alternative, method, crts, TMTval=NA, TMTi=NA){
  names(T)="T"
  RVAL <- list(statistic = T, p.value = p, alternative = alternative, 
               method = method, data.name = data.name, t=crts, TMTval=TMTval, TMTi=TMTi)
  class(RVAL) <- "htest"
  return(RVAL)
}

corrtest <- function(A, B, start=0, end=1, F=NULL, alternative="causes", method="timeout", transform=FALSE, usebeforeA1=FALSE, maxtau=NA, careful=TRUE){
    data.name=paste(deparse(substitute(A)), deparse(substitute(B)))
    if (careful){
        if (length(A) == 0 | length(B)==0){stop("A or B is empty")}
        if (is.unsorted(A)){A=sort(A)}
        if (is.unsorted(B)){B=sort(B)}
        if (B[1] < start | A[1] < start | B[length(B)]>end | A[length(A)]>end) stop("A and/or B are not within [start,end]")        
        if (length(B)>1 & min(diff(B))==0) warning("B has duplicate entries, results could be spurious")
        if (!is.null(F)){
            trueStart=ifelse(usebeforeA1, start, A[1])
            testPoints = c(trueStart,sort(runif(10, min=trueStart, max=end)), end)
            values = sapply(testPoints, function(point) F(point))
            if (is.unsorted(values)){
                stop("the supplied F must be non-decreasing over [start, end] if usebeforeA1 and [A[1], end] otherwise")                
            }            
            FA=F(A); FB=F(B); if (((length(FA)>1)&min(diff(FA))==0) | (length(FB)>1 & min(diff(FB))==0)) warning("F(A) or F(B) has duplicate entries (F is not strictly increasing), results could be spurious")
            if (sum(duplicated(c(FA[!duplicated(FA)], FB[!duplicated(FB)])))>=1){
                warning("F(a)=F(b) for some a,b (F is not strictly increasing) returning (probably spuriously) p=0")
            }
        }
    }
    if (!usebeforeA1 & length(B[B>=A[1]])==0){return(htest(p=1, T=0, data.name=data.name, alternative=alternative, method=method, crts=B))}
    if (transform){
        A=F(A); B=F(B); Fstart=F(start); Fend=F(end); A = (A-Fstart)/(Fend-Fstart); B = (B - Fstart)/(Fend-Fstart); F=NA; start=0; end=1
    }
    if (usebeforeA1 & method=="timeout"&(alternative=="inhibits" | alternative=="causes")){
        ##don't want the timeout to wrap around
        if (is.na(maxtau)) maxtau=Inf
        maxtau=min(maxtau, max(diff(c(A,end))))
    }

    if (alternative=="causes" | alternative=="inhibits"){
        if (!is.null(F)){
            Astar=c(A,end)    
            crts=sapply(B[B>=A[1]], function(b){
                tau=b-max(A[A<=b])
                sum(sapply(1:length(A), function(i){
                    max(0,F(min(Astar[i+1],A[i]+tau))-F(A[i]))
                }))
            })
            if (method=="timeout" & !is.na(maxtau)){
                maxtaut= sum(sapply(1:length(A), function(i){
                    max(0,F(min(Astar[i+1],A[i]+maxtau))-F(A[i]))
                }))
            }
            R=F(end)-F(A[1])
            if (usebeforeA1){
                crts=c(R-F(start)+F(B[B<A[1]]), crts);
                total=F(end)-F(start);
                crts=crts/total
            }
            else{
                crts=crts/R; if (method=="timeout" & !is.na(maxtau)){maxtaut=maxtaut/R}
            }
        }
        else {
            ##Can go a little faster
            Bstar=B[B>=A[1]]
            rt=sapply(Bstar, function(b){b - max(A[A<=b])})
            Deltas=diff(c(A,end))
            crts=sapply(1:length(Bstar), function(i){sum(pmin(Deltas, rt[i]))})
            if (method=="timeout" & !is.na(maxtau)){
                maxtaut=sum(pmin(Deltas, maxtau))
            }
            if (usebeforeA1) {
                crts=c(end-A[1]+B[B<A[1]]-start, crts); crts=crts/(end-start)
            }
            else {
                crts=crts/(end-A[1]); if (method=="timeout" & !is.na(maxtau)){maxtaut=maxtaut/(end-A[1])}
            }
        }
    }    
    else {##alternatives are correlation or anticorrelation
        if (length(A)==1) mids=c(start, end)        
        else mids=c(start, (A[1:(length(A)-1)]+A[2:length(A)])/2, end)
        if (!is.null(F)){
            crts=sapply(B, function(b){
                tau=min(abs(b-A))
                sum(sapply(1:length(A), function(i){
                    F(min(mids[i+1], A[i]+tau)) - F(max(mids[i], A[i]-tau))
                }))
            })
            if (method=="timeout"&!is.na(maxtau)){
                maxtaut=sum(sapply(1:length(A), function(i){
                    F(min(mids[i+1], A[i]+maxtau)) - F(max(mids[i], A[i]-maxtau))
                }))
            }
            crts=crts/(F(end)-F(start))
      }
        else {
            Deltas = diff(sort(c(start,A,mids, end)))
            rt=sapply(B, function(b){min(abs(A-b))})
            crts=sapply(rt, function(x){sum(pmin(Deltas,x))})
            if (method=="timeout"&!is.na(maxtau)){
                maxtaut=sum(pmin(Deltas,maxtau))
            }
            crts=crts/(end-start)
        }
    }
    if (alternative=="inhibits" | alternative=="anticorrelation") crts=1-crts

    if (min(crts)==0) return(htest(p=0, T=Inf, data.name=data.name, alternative=alternative, method=method, crts=crts))
    TMTval=NA; TMTi=NA
    if (method=="timeout") {
        maxtaut=ifelse(is.na(maxtau), 1.1, maxtaut)
        stat=TMT.test(crts, maxtau=maxtaut)
        p=stat$p; T = stat$lk; TMTi=stat$TMTi
        if (!is.na(TMTi)){
            Bi=order(crts)[TMTi]+(ifelse(usebeforeA1, 0, sum(B<A[1])))
            TMTval = B[Bi]-max(A[A<=B[Bi]])
        }
    }
    if (method=="simes"){
        pT=simes.test(crts, returnstat=TRUE); p=pT[1]; T=pT[2]
    }    
    if (method=="fisher"){
        pT=F.test(crts, returnstat=TRUE); p=pT[1]; T=pT[2]
    }
    htest(p,T,data.name,alternative,method,crts, TMTval, TMTi)
}

## A=c(.1, .2, .3)
## B=c(.09, .22, .31)

## res=corrtest(A,B)
## corrtest(A,B, usebeforeA1=TRUE); corrtest(A,B, usebeforeA1=TRUE)$t
## corrtest(A,B, method="simes", usebeforeA1=TRUE); corrtest(A,B, method="simes", usebeforeA1=TRUE)$t
## res$t

## a


## events = list(c(0.1, 0.2, 0.3), c(.4), c(.05, .11, .21, .32), c(.06, 0.22, .33))
## stack = list(c(1,2,3), c(3,4)) #means we think 1 and 2 cause 3, and 3 causes 4
## a = ttransform(events, stack)
## stack = list(c(1,2,3), c(1,3,4)) #means we think 1 and 2 cause 3, and 1 and 3 cause 4
## b = ttransform(events, stack)
## print(a)
## print(b)

## events = list(c(.1, .2, .3), c(.05, .11, .21, .31), c(.06, .12, .22, .32))
## stack = list(c(1,2), c(2,3))
## e1 = ttransform(events, stack, multipass=TRUE)
## e2 = ttransform(events, stack, multipass=FALSE)
## print(e1)
## print(e2)
