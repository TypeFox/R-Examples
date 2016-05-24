## Functions to Reproduce tests from Rohmel and Kieser (2013, Stat in Med, 2335-2348)
##      We always use odds ratios equal to 
##        OR<-((threshold+delta)*(1-threshold))/(threshold*(1-threshold-delta))
##        where delta=difference margin and threshold is the threshold for 
##        the control rate for switching 
##        between odds ratio and difference margins 
## NiM1/T1 = fmecTest with type="threshold" (except we use OR=1.71, see above)
## NiM3/T2 = fmecTest with type="max" 
## NiM3/T3 = fmecExact with type="max"
## NiM3/T4 = brkTest

FarrMann<-function(X1,X2,M1,M2,delta0,alternative=c("less","greater","two.sided"),fudgeDigits=8){
        alt<-match.arg(alternative)
        ## Got this from Farrrington and Manning, 
        ## Stat in Med 1990, 1447-1454
        ## they use Theta1-Theta2=delta
        ## instead of Theta2-Theta1=delta
        ## so switch
        i<-X2
        N1<-M2
        j<-X1
        N2<-M1        
        P1hat<- i/N1
        P2hat<-j/N2
        #if (delta0==0){
        #    p1d<-P1hat
        #    p2d<-P2hat
        #} else {
        ### get MLEs (p1d and p2d) from Appendix of Farrington
        ### and Manning, Stat in Med, 1990, 1447-1454
        theta<-N2/N1
        a<-1+theta
        s0<- delta0
        b<- -1*(1+theta+P1hat+theta*P2hat+s0*(theta+2) )
        cc<- s0^2 + s0*(2*P1hat+theta+1) + P1hat + theta*P2hat
        d<- -P1hat*s0*(1+s0)
        v<- b^3/(3*a)^3 - (b*cc)/(6*a^2) + d/(2*a)
        u<- sign(v)*(b^2/(3*a)^2 - cc/(3*a) )^(1/2)
        temp<-v/u^3
        ## define 0/0=1 to avoid NaNs messing things up
        temp[v==0 & u==0]<-1
        temp<-pmax(-1,temp)
        temp<-pmin(1,temp)
        w<-(1/3)*(pi+acos(temp) )
        p1d<- 2*u*cos(w) - b/(3*a)
        p1d<-pmax(s0,p1d)
        p1d<-pmin(1,p1d)
        p2d<-p1d-s0
        p1d<-round(p1d,fudgeDigits)
        p2d<-round(p2d,fudgeDigits)
        #}
        ### Now do Z 
        ### round std.err and numerator to set fuzz to 0
        std.err<-   round(    ( (p1d*(1-p1d))/N1 + 
           (p2d*(1-p2d))/N2 )^(1/2), fudgeDigits)
        numerator<- round( (P1hat-P2hat-delta0), fudgeDigits)
        #if (any(is.na(std.err)) | any(std.err<0)) browser()     
        Z<- numerator/std.err
        ## set no effect to zero even as std.err goes to zero
        Z[numerator==0 & std.err==0]<-0
        #
        #out<-list(Z=Z,p.value=pnorm(Z))
    # one-sided p-value
    pval<-pnorm(Z)
    if (alt=="less"){
        pval<-pval
    }  else if (alt=="greater"){
        pval<-1-pval
    } else if (alt=="two.sided"){
        pval<-min(1,2*pval)
    }
    pval
}



EC<-function(X1,X2,M1,M2,or=1.71,alternative=c("less","greater","two.sided")){
    alt<-match.arg(alternative)
    fisher.test(matrix(c(X2,M2-X2,X1,M1-X1),2,2),or=or,alternative=alt)$p.value
}


fmecTest<-function(x1,n1,x2,n2,threshold=0.2,delta=0.1,
    alternative=c("less","greater"),
    type=c("max","switch")){
    alt<-match.arg(alternative)
    OR<-((threshold+delta)*(1-threshold))/(threshold*(1-threshold-delta))
    pEC<-EC(x1,x2,n1,n2,or=OR,alternative=alt)
    pFM<-FarrMann(x1,x2,n1,n2,delta,alternative=alt)
    #list(pEC=pEC,pFM=pFM)
    type<-match.arg(type)
    if (type=="max"){
        pval<-max(pEC,pFM)
    } else if (type=="switch"){
        pval<-ifelse(x1/n1<threshold,pEC,pFM)
    } else { stop("type should be either 'max' or 'threshold' ") }
    METHOD<-paste0("Odds Ratio/Difference Test \n type=",type)
    dname<-paste0("x1=",x1," n1=",n1," x2=",x2," n2=",n2)
    statistic<-c(threshold,delta,OR)
    names(statistic)<-c("threshold","difference","odds ratio")
    null.value<-delta
    names(null.value)<-"difference in proportions"


    out<-list(statistic=statistic,data.name=dname,method=METHOD,
          p.value=pval,alternative=alt,null.value=null.value)
    class(out)<-"htest"
    out
}


fmecExact<-function(x1,n1,x2,n2,threshold=0.2,delta=0.1,
    alternative=c("less","greater"),
    type=c("max","switch"),
    ngrid=1000){
    alternative<-match.arg(alternative)
    OR<-((threshold+delta)*(1-threshold))/(threshold*(1-threshold-delta))
    # See equation 9 of Rohmel and Kieser (2013, Stat in Med, 2335-2348
    Tmat<-matrix(NA,n1+1,n2+1)
    type<-match.arg(type)
    for (i in 0:n1){
        for (j in 0:n2){
            Tmat[i+1,j+1]<- fmecTest(x1=i,n1=n1,x2=j,n2=n2,
               threshold=threshold,delta=delta,
               alternative=alternative,type=type)$p.value
        }
     }
     Delta<-delta
     NIM<-function(p){ nimDiffOR(p,q=threshold,delta=Delta) }  
     pp<- 0:ngrid/ngrid
     qq<- NIM(pp)
     pvalVector<-rep(NA,ngrid+1)
     for (j in 1:(ngrid+1)){
         ProbMat<- matrix(dbinom(0:n1,n1,pp[j]),n1+1,1) %*% 
                   matrix(dbinom(0:n2,n2,qq[j]),1,n2+1)
         pvalVector[j]<- sum( ProbMat[Tmat <= Tmat[x1+1,x2+1]] ) 
     }
     p.value<- max(pvalVector)
    METHOD<-paste0("Exact Odds Ratio/Difference Test \n type=",type)
    dname<-paste0("x1=",x1," n1=",n1," x2=",x2," n2=",n2)
    statistic<-c(threshold,delta,OR)
    names(statistic)<-c("threshold","difference","odds ratio")
    functionName<-deparse(substitute(g))
    null.value<-delta
    names(null.value)<-"difference in proportions"

    out<-list(statistic=statistic,data.name=dname,method=METHOD,null.value=null.value,
       alternative=alternative,p.value=p.value)
    class(out)<-"htest"
    out
}


#fmecTest(6,10,2,12,alternative="less",type="max")
#fmecExact(6,10,2,12,alternative="less",type="max")

getij<-function(R){
    n1<-nrow(R)-1
    n2<-ncol(R)-1
    imat<-matrix(rep(1:(n1+1),n2+1),n1+1,n2+1)
    jmat<-matrix(rep(1:(n2+1),each=n1+1),n1+1,n2+1)
    ivec<-as.vector(imat[R])
    jvec<-as.vector(jmat[R])
    list(i=ivec,j=jvec)
}



findPower<-function(ir,jr,n1,n2,psearch,qsearch){
    # find power when F1[i]=psearch[i] and F2[i]=qsearch[i]
    p<-psearch
    q<-qsearch
    pow<-rep(NA,length(p))
    X1<-0:n1
    X2<-0:n2
    for (i in 1:length(p)){
        f1<-dbinom(X1,n1,p[i])
        f2<-dbinom(X2,n2,q[i])
        pow[i]<-sum(f1[ir]*f2[jr])
      }
    pow
}


getieje<-function(ir,jr){
    ui<-sort(unique(ir))
    uj<-sort(unique(jr))
    iout<-jout<-rep(NA,length(ui)+length(uj))
    for (j in 1:length(uj)){
        jout[j]<-uj[j]
        iout[j]<-min(ir[jr==uj[j]])
    }
    b<-length(uj)
    for (i in 1:length(ui)){
        iout[b+i]<-ui[i]
        jout[b+i]<-max(jr[ir==ui[i]])
    }
    o<-order(iout,jout)
    iout<-iout[o]
    jout<-jout[o]
    idiff<-c(1,diff(iout))
    jdiff<-c(1,diff(jout))
    extras<- (idiff==0 & jdiff==0)
    iout<-iout[!extras]
    jout<-jout[!extras]
    list(ie=iout,je=jout)
}


getPossibleR<-function(ie,je,n1,n2,decreasei=TRUE,increasej=TRUE){
    #ie<-inlist$ie
    #je<-inlist$je
    #n1<-inlist$n1
    #n2<-inlist$n2
    # input ie=i index on the edge of the rejection region
    #       je=j index on the edge of the rejection region
    if (!decreasei | !increasej) stop("need to rewrite program to handle increasing i or decreasing j")
    if (max(je)==n2+1){
       jnew<-min(je):max(je)
    } else {
       jnew<-min(je):(max(je)+1)
    }
    inew<-rep(NA,length(jnew))
    for (a in 1:length(jnew)){
        imatch<- ie[je==jnew[a]]
        if (length(imatch)==0){
            # if no values of i match je=jnew[a]
            # then this is added on to the end
            # so inew[a] must be the largest value
            inew[a]<-n1+1
        } else {
            inew[a]<- min(imatch)-1
            if (inew[a]<1){
                jnew[a]<-NA
            }
        }    

        if (a>1 && inew[a]==inew[a-1]){
            # do not want to keep larger jnew at same inew level
            # do not set inew[a]<-NA, so we can run the if statement 
            # for the next iteration
            jnew[a]<-NA
        }     
    }
    keep<- !is.na(inew) & !is.na(jnew)
    iout<-inew[keep]
    jout<-jnew[keep]
    list(inew=iout,jnew=jout)
}


getpadd<-function(i,j,n1,n2,psearch,qsearch){
    dbinom(i-1,n1,psearch)*dbinom(j-1,n2,qsearch)
}


getR<-function(i,j,n1,n2){
    ivec<-rep(1:(n1+1),n2+1)
    jvec<-rep(1:(n2+1),each=n1+1)
    Rvec<-rep(FALSE,length(ivec))
    ni<-length(i)
    for (h in 1:ni){
        Rvec[ivec==i[h] & jvec==j[h]]<-TRUE
    }
    R<-matrix(Rvec,n1+1,n2+1)
    R
}


getfij<- function(n1,n2,p,q){
    matrix(dbinom(0:n1,n1,p),n1+1,1) %*% matrix(dbinom(0:n2,n2,q),1,n2+1)
}

findPowerR<-function(R,g,psearch=(0:1000)/1000){
    n1<-nrow(R)-1
    n2<-ncol(R)-1
    p<-psearch
    q<-g(p)
    pow<-rep(NA,length(p))
    for (i in 1:length(p)){
        fij<-getfij(n1,n2,p[i],q[i])
        check<-sum(fij)
        if (abs(check-1)>0.0001) stop("check fij function")
        pow[i]<-sum(fij[R])
    }
    d<-data.frame(p=p,q=q,power=pow)
    #d[pow==max(pow),]
    #list(d=d,maxpower=max(pow))
    #max(pow)
    d
}

brkCalc<-function(n1,n2,threshold,delta,g=nimDiffOR,alpha=0.025,alphastar=0.001,ngrid=1000){

    T2<-matrix(NA,n1+1,n2+1)
    for (i in 0:n1){
        for (j in 0:n2){
            T2[i+1,j+1]<- fmecTest(i,n1,j,n2,threshold=threshold,
                    delta=delta,
                    alternative="less",type="max")$p.value
       }
    }



    psearch<-0:ngrid/ngrid
    qsearch<-g(psearch,delta=delta,q=threshold)
    n1<-nrow(T2)-1
    n2<-ncol(T2)-1
    # get a current rejection region
    R<-matrix(FALSE,nrow(T2),ncol(T2))
    R[T2<=alphastar]<-TRUE
    # get indices for rejection region in R, then get p-value vector
    currij<-getij(R)
    ci<-currij$i
    cj<-currij$j
    cp<-findPower(ci,cj,n1,n2,psearch,qsearch)
    prepmat<-prei<-prej<-NULL
    # create list of input for iterating function, stepi
    # prepmat, prei, and prej are used later, set to NULL now

    # since computers have problems with ties, round pvalues to rdig digits
    rdig<- (-1)*ceiling(log10(alphastar)) + 3
    PVAL<-matrix(NA,nrow(T2),ncol(T2))
    ## do not calculate p-value for each values in the starting rejection 
    ## region, just put the upper bound for every value
    PVAL[R]<- round(max(cp),rdig)
    ## SYMB gives the symbol in the expression: p.value SYMB PVAL
    ## so the starting rejection region will have: p.value <= PVAL
    SYMB<-matrix(NA,nrow(T2),ncol(T2))
    SYMB[R]<-"<="
    ## later SYMB will be "=" as we add values to rejection region, 
    ##  and at the end, all values that fail to reject will be ">"
    inList<-list(ci=ci,cj=cj,cp=cp,
        prepmat=NULL,
        prei=NULL,
        prej=NULL,p.value=round(max(cp),rdig),
        PVAL=PVAL, SYMB=SYMB)

    stepi<-function(L){
        ci<-L$ci; cj<-L$cj; cp<-L$cp; 
        prepmat<-L$prepmat; prei<-L$prei ; prej<-L$prej
        PVAL<-L$PVAL
        SYMB<-L$SYMB
        # get indices at the edge of rejection region of R
        ije<-getieje(ci,cj)
        # get indeces for possible next additions to rejection region
        newije<-getPossibleR(ije$ie,ije$je,n1,n2)
        newie<-newije$inew
        newje<-newije$jnew
        paddmat<-matrix(NA,length(newie),length(psearch))
        newpe<-rep(NA,length(newie))
        # find additional point that increases p-value smallest
        # amount, ties go to smaller j (new treatment failures)
        # if padd has already been calculated, do not repeat the calculation
        if (is.null(prepmat)){
            for (h in 1:length(newie)){          
                paddmat[h,]<-  getpadd(newie[h],newje[h],
                                    n1,n2,psearch,qsearch)
                newpe[h]<- max(cp+paddmat[h,])
            }
        } else {
            for (h in 1:length(newie)){
                I<-prei==newie[h] & prej==newje[h]

                if (any(I)){        
                    paddmat[h,]<- prepmat[I,]
                } else {
                    paddmat[h,]<-  getpadd(newie[h],newje[h],n1,n2,psearch,qsearch)
                }
                newpe[h]<- max(cp+paddmat[h,])
            }
        }
        # since computers have problems with ties, round
        newpe<-round(newpe,rdig)
        # pick out index with smallest p-value, ties pick smaller j
        minp<-min(newpe)
        jatminp<- newje[newpe==minp]
        bestj<- min(jatminp)
        pick<-newje==bestj & newpe==minp
        besti<- newie[pick]
        ## June 19, 2015: change second if statement to <1, instead of ==1
        if (length(newie[!pick])==1){
             prepmat<-matrix(paddmat[!pick,],nrow=1)
        } else if (length(newie[!pick])<1){
             prepmat<-NULL
        } else {
             prepmat<-paddmat[!pick,]
        }
        PVAL[besti,bestj]<- minp
        SYMB[besti,bestj]<-"="
        outList<-list(ci=c(ci,besti),
            cj=c(cj,bestj),
            cp=cp + paddmat[pick,],
            prepmat=prepmat,
            prei=newie[!pick],
            prej=newje[!pick],p.value=minp, PVAL=PVAL, SYMB=SYMB)
        outList
    }
    # stepi= gets next rejection point, given current and precalculated paddmat
    p.value<-inList$p.value
    iter<-0
    while(p.value<=alpha){
        iter<-iter+1
        outList<-stepi(inList)
        p.value<-outList$p.value
        #cat("iter=",iter)
        #cat(" p.value=",p.value,"\n")
        if (p.value<=alpha) inList<-outList
    }
    ## get SYMB and PVAL from inList
    SYMB<-inList$SYMB
    SYMB[is.na(SYMB)]<-">"
    PVAL<-inList$PVAL
    PVAL[is.na(PVAL)]<-alpha
    #out<-c(list(R=getR(inList$ci,inList$cj,n1,n2)),
    #  inList,
    # list(next.pvalue=p.value))

    matname<-list(paste0("x1=",0:n1),paste0("x2=",0:n2))
    R<-getR(inList$ci,inList$cj,n1,n2)
    attr(R,"sig.level")<-alpha
    PVALUES<-matrix(paste0("p",SYMB,PVAL),n1+1,n2+1)
    dimnames(R)<-dimnames(PVALUES)<-dimnames(PVAL)<-dimnames(SYMB)<-matname

    out2<-list(R=R,
      PVALbounds=PVAL, 
      PVALsymbols=SYMB,
      PVALUES=PVALUES)

     out2
}


brkControl<-function(alpha=0.025,alphastar=0.001,ngrid=1000){
    list(alpha=alpha,alphastar=alphastar,ngrid=ngrid)
}

brkTest<-function(x1,n1,x2,n2,threshold=0.2,delta=0.1,
        control=brkControl()){
    # get pvalues from fmecTest(type="max") to get 
    # starting rejection region
    if (n1<5 | n2<5) stop("n1 and n2 must be at least 5")
    bout<-brkCalc(n1,n2,threshold=threshold,delta=delta,g=nimDiffOR,alpha=control$alpha,alphastar=control$alphastar,
                ngrid=control$ngrid)

    symb<-bout$PVALsymbols[x1+1,x2+1] 
    if (symb==">"){ 
        p.value<- 1
    } else if (symb=="=" | symb=="<="){
        p.value<-bout$PVALbounds[x1+1,x2+1]
    }
    METHOD<-paste0("Barnard-Rohmel-Kieser Test")
    dname<-paste0("x1=",x1," n1=",n1," x2=",x2," n2=",n2)
    OR<-((threshold+delta)*(1-threshold))/(threshold*(1-threshold-delta))
    statistic<-c(threshold,delta,OR,x1,x2)
    names(statistic)<-c("threshold","difference","odds ratio","x1","x2")
    #functionName<-deparse(substitute(g))
    out<-list(FullResults=bout,
            statistic=statistic,data.name=dname,method=METHOD,p.value=p.value)
    class(out)<-c("brk","htest")
    out
}


print.brk<-function(x, digits=getOption("digits"), prefix="\t",...){
    cat("\n")
    cat(strwrap(x$method, prefix = prefix), sep = "\n")
    cat("\n")
    cat("data:  ", x$data.name, "\n", sep = "")

    out<- character() 
    out <- c(out, paste(names(x$statistic[1:3]), "=", format(signif(x$statistic[1:3], 
            max(1L, digits - 2L)))))
    cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
 
    x1<- x$statistic["x1"]
    x2<- x$statistic["x2"]

    cat("\n")
    cat("One-sided p-value:")
    cat("\n")
    cat(x$FullResults$PVALUES[x1+1,x2+1])
    cat("\n")
    if (x$FullResults$PVALsymbols[x1+1,x2+1]!="="){
       cat("Note: to save computational time, only bound on p-value calculated.")
       cat("\n")
    }
    cat("\n")
    cat("For rejection region and p-values for any possible")
    cat("\n")
    cat("result with these sample sizes, save output as x,")
    cat("\n")
    cat("and see the list x$FullResults")
    cat("\n")
    invisible(x)
}


#x<-brkTest(3,6,0,5)

