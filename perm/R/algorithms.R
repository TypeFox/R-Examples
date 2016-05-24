### add p.conf.level to all .exact.mc
### need to add it to the permXX also

`twosample.pclt` <-
function(scores,group){
    tab<-table(group,scores)
    m<-sum(tab[2,])
    n<-length(scores)
    ## note that twosample.pclt is not called directly, but from within permTS
    ## so we know that the first group has group=1 and the second has group=0
    ##
    ## we add the following in case 
    ## it is called directly, then it gives a correct answer even if group is 
    ## defined with characters
    ## note "1"==1 is TRUE
    ## (note it even works if group is continuous numeric since if y<-rnorm(1), then as.character(y)==y is TRUE)
    Grp1<-dimnames(tab)[[1]][2]
    grp<-rep(0,n)
    grp[group==Grp1]<-1
    T0<-sum( scores*grp) 
    SSE.scores<- sum( (scores - mean(scores))^2 )  
    SSE.grp<- sum( (grp - mean(grp))^2 )  

    Z<-sqrt(n-1)*(T0 - n*mean(scores)*mean(grp))/
                  sqrt(SSE.scores*SSE.grp)
    p.lte<-pnorm(Z)
    p.gte<-1-pnorm(Z)
    ## in this case both twosided pvalues will always be equal
    p.twosidedAbs<- 1-pchisq(Z^2,1)
    p.values<-c(p.twosided=min(1,2*min(p.lte,p.gte)),p.lte=p.lte,p.gte=p.gte,p.twosidedAbs=p.twosidedAbs)
    out<-list(p.values=p.values,Z=Z)
    out
}

`twosample.exact.network` <-
function(scores,group,digits=12){
    tab<-table(group,scores)
    M<-tab[2,]
    y<-tab[1,]
    ## T is a vector of totals for each unique score
    T<-M+y
    ## V is a matching vector of the values of those unique scores
    V<-as.numeric(dimnames(tab)[[2]])
    ## test statistic is sum of scores in one group
    ## an equivalent test statistic (i.e., gives same p-values) 
    ## is the difference in group means of the scores
    T0<- sum( V*M )
    k<-length(M)
    m<-sum(M)
    n<-sum(T)-m
    N<-sum(T)
    ## for notation and terminology see agresti, mehta, and patel (1990, JASA, 453-458)
    wj<-0
    ## create vectors too big on purpose, growing vectors takes more time
    stage<-Wj<-rep(NA,k*N)
    stage[1]<-Wj[1]<-0
    cnt<-2
    ## first go through and define all nodes
    for (j in 0:(k-1)){
        temp.wj<-min(pmax(wj,m-N+sum(T[1:(j+1)]))):
            max(pmin(m,wj+T[j+1]))
        n.wj<-length(temp.wj)
        stage[cnt:(cnt+n.wj-1)]<-rep(j+1,n.wj)
        Wj[cnt:(cnt+n.wj-1)]<-temp.wj
        wj<-temp.wj
        cnt<-cnt+n.wj
    }
    ## not delete parts of stage and Wj not used
    stage<-stage[!is.na(stage)]
    Wj<-Wj[!is.na(Wj)]
    nNodes<-length(stage)
    Nodes<-1:nNodes
    Edges<-probL<-rankL<-rep(NA,nNodes*(max(T)+1))
    cnt<-1
    edgesize<-cnt.edge<-rep(NA,nNodes)
    cnt.stage<-rep(1,k+1)
    cnt.edge[1]<-1
    for (i in 1:(nNodes-1)){
            wj<-Wj[i]  
            temp<-max(wj,m-N+sum(T[1:(stage[i]+1)])):
                min(m,wj+T[stage[i]+1])
            edgesize[i]<-length(temp)
            ## get number of nodes to go with (i,temp)
            pick<-rep(FALSE,nNodes)
            for (h in 1:length(temp)){
                pick[stage==stage[i]+1 & Wj==temp[h]]<-TRUE
            }
            cnt2<-cnt+edgesize[i]
            Edges[cnt:(cnt2-1)]<-Nodes[pick]
            probL[cnt:(cnt2-1)]<- choose(T[stage[i]+1],temp-Wj[i])
            rankL[cnt:(cnt2-1)]<-(temp-Wj[i])*V[stage[i]+1]
            cnt<-cnt2
            cnt.edge[i+1]<-cnt
            if (stage[i+1]-stage[i]>0) cnt.stage[stage[i+1]+1]<-cnt2
    }

    #g<-list(stage=stage,Wj=Wj,Edges=Edges,wgts=wgts,probL=probL,rankL=rankL)

        thk<-rankL[cnt.stage[1]:(cnt.stage[2]-1)]
        Chk<-probL[cnt.stage[1]:(cnt.stage[2]-1)]
        nodehk<-Edges[cnt.stage[1]:(cnt.stage[2]-1)]
        for (i in 2:k){
            ## most of the time that this function takes comes from the getcnt function
            ## to do: move that function to C code, or speed it up some other way
            CNT<-getcnt(nodehk,cnt.edge,edgesize)
            newthk<-rep(thk,times=edgesize[nodehk]) + 
                rankL[CNT]
            newChk<-rep(Chk,times=edgesize[nodehk])*
                probL[CNT]
            newnodehk<-Edges[CNT]
            thk<-newthk
            Chk<-newChk
            nodehk<-newnodehk
        }   

        out<-calcPvals(thk,T0,digits,denom=choose(N,m),wgts=Chk)
        out
}

`twosample.exact.mc` <-
function(scores,group,alternative="two.sided",nmc=10^4-1,seed=1234321,digits=12,p.conf.level=.99,setSEED=TRUE){
    t0<-sum(scores*group)
    N <- nmc
    if (setSEED) set.seed(seed)
    ti <- rep(NA, N)
    for (i in 1:N) {
        ti[i] <- sum(scores*sample(group))
    }    
    out<-calcPvalsMC(ti,t0,digits,alternative,FALSE,p.conf.level)
    out
}

`twosample.exact.ce` <-
function(scores,group,cm=NULL,digits=12){
    tab<-table(group,scores)
    m<-sum(tab[2,])
    n<-length(scores)
    if (is.null(cm)) cm<-chooseMatrix(n,m)

    Grp1<-as.character(dimnames(tab)[[1]][2])
    grp<-rep(0,n)
    grp[as.character(group)==Grp1]<-1
    Tj<- cm %*% scores
    T0<-sum( scores*grp) 
    out<-calcPvals(Tj,T0,digits)
    out
}
calcPvals<-function(Tj,T0,digits,denom=NULL,wgts=NULL,twosidedTstat=FALSE){
    if (is.null(denom)) nTj<-length(Tj)
    if (is.null(wgts)) wgts<-rep(1,nTj)
    ## subtract off mean, needed only for p.twosidedAbs
    denom<-sum(wgts)
    mu<-sum(wgts*Tj)/denom
    Tj<-Tj-mu
    T0<-T0-mu
    ## to avoid problems with computer rounding error
    Tj<-signif(Tj,digits)
    T0<-signif(T0,digits)
    if (twosidedTstat){
        p.values<-c(p.twosided=sum(wgts[Tj>=T0])/denom,
                  p.equal=sum(wgts[Tj==T0])/denom)
    } else {
        p.lte<-sum(wgts[Tj<=T0])/denom
        p.gte<-sum(wgts[Tj>=T0])/denom
        p.equal<-sum(wgts[Tj==T0])/denom
        ## p.twosidedAbs is alternative="two.sided" with control=permControl(tsmethod="abs")
        ## p.twosided is alternative="two.sided" with control=permControl(tsmethod="central")
        p.twosidedAbs<-sum(wgts[abs(Tj)>=abs(T0)])/denom
        p.values<-c(p.twosided=min(1,2*min(p.lte,p.gte)),p.twosidedAbs=p.twosidedAbs,
            p.lte=p.lte,p.gte=p.gte,p.equal=p.equal)
    }
   out<-list(p.values=p.values)
   out
}
calcPvalsMC<-function(Tj,T0,digits,alternative=NULL,twosidedTstat=FALSE,p.conf.level=.99){
    N<-length(Tj)
    wgts<-rep(1,N)
    ## subtract off mean, needed only for p.twosidedAbs
    mu<-mean(Tj)
    Tj<-Tj-mu
    T0<-T0-mu
    ## to avoid problems with computer rounding error
    Tj<-signif(Tj,digits)
    T0<-signif(T0,digits)

    S.lte <- sum(wgts[Tj <= T0])
    S.gte <- sum(wgts[Tj >= T0])
    S.abs <- sum(wgts[abs(Tj) >= abs(T0)])
  
    ## always count the observed value, this ensures valid p-values 
    ## see Fay, Kim and Hachey, 2007, JCGS, 946-967, eq 5.3
    p.lte <- (S.lte + 1)/(N + 1)
    p.gte <- (S.gte + 1)/(N + 1)
    p.equal<-  (length((1:N)[Tj == T0])+1)/(N+1)
    p.twosidedAbs<- (S.abs+1)/(N+1)

    if (twosidedTstat & !is.null(alternative)) warning("twosidedTstat=TRUE so alternative ignored")

    if (twosidedTstat){
        p.conf.int<-pCI(S.gte,N,p.conf.level)
        p.values<-c(p.twosided=p.gte,p.equal=p.equal)
    } else {
        if (alternative=="less"){ 
            p.conf.int<-pCI(S.lte,N,p.conf.level)
        } else if (alternative=="greater"){
            p.conf.int<-pCI(S.gte,N,p.conf.level)
        } else if (alternative=="two.sidedAbs"){
            p.conf.int<-pCI(S.abs,N,p.conf.level)
        } else if (alternative=="two.sided"){
            if (S.lte<S.gte){
                p.conf.int<-2*pCI(S.lte,N,p.conf.level)
            } else p.conf.int<-2*pCI(S.gte,N,p.conf.level)
            if (p.conf.int[2]>1) p.conf.int[2]<-1
        } 
         p.values<-c(p.twosided=min(1,2*min(p.lte,p.gte)),p.twosidedAbs=p.twosidedAbs,
            p.lte=p.lte,p.gte=p.gte,p.equal=p.equal)
    }
    out<-list(p.values=p.values,p.conf.int=p.conf.int)
    out
}

pCI<-function(S,N,p.conf.level=.99){
    p.alpha<-1-p.conf.level
    # p<-(S+1)/(N+1)
    ci.lower<-ifelse(S == 0,0,qbeta(p.alpha/2, S, N - S + 1))
    ci.upper<-ifelse(S==N,1,qbeta(1 - p.alpha/2, S + 1, N - S))
    pvalue.ci <- c(ci.lower, ci.upper)
    attr(pvalue.ci, "conf.level") <- p.conf.level
    pvalue.ci
}


`ksample.pclt` <-
function(scores,group){
    ### (see e.g., fay and shih, 1998, JASA, p. 389, eq 4)
    if (is.factor(group)) group<-as.character(group)
    ug<-unique(group)
    ng<-length(ug)
    n<-length(scores)
    ## standardize scores so that they sum to zero
    scores<-scores - mean(scores)
    SSE.scores<- sum( scores^2 )  
    N<-rep(NA,ng)
    names(N)<-ug
    mean.scores<-N
    for (j in 1:ng){
        mean.scores[j]<-mean(scores[group==ug[j]])
        N[j]<-length(scores[group==ug[j]])
    }
    chisq.value<- ((n-1)/SSE.scores)*sum( N*(mean.scores^2) )
    p.twosided <- 1 - pchisq(chisq.value, ng - 1)
    p.values<-c(p.twosided=p.twosided)
    out<-list(p.values=p.values,chisq.value=chisq.value,df=ng-1)
    out
}
`ksample.exact.mc` <-
function(scores,group,nmc=10^4-1,seed=1234321,digits=12,p.conf.level=.99,setSEED=TRUE){
    scores<-scores-mean(scores)
    if (is.factor(group)) group<-as.character(group)
    ug<-unique(group)
    ng<-length(ug)
    calcTestStat<-function(x,g,Ng=ng,Ug=ug){
        mean.scores<-N<-rep(NA,Ng)
        for (j in 1:Ng){
            mean.scores[j]<-mean(x[g==Ug[j]])
            N[j]<-length(x[g==Ug[j]])
        }
        sum( N*(mean.scores^2) )
    }

    t0<-calcTestStat(scores,group)
    N <- nmc
    if (setSEED) set.seed(seed)
    ti <- rep(NA, N)
    for (i in 1:N) {
        ti[i] <- calcTestStat(scores,sample(group,replace=FALSE))
    }    
    out<-calcPvalsMC(ti,t0,digits,NULL,TRUE,p.conf.level)
    out
}


`trend.pclt` <-
function(scores,group){
    if (!is.numeric(group)) stop("for trend tests group must be numeric")
    T0<-sum( scores*group) 
    SSE.scores<- sum( (scores - mean(scores))^2 )  
    SSE.grp<- sum( (group - mean(group))^2 )  
    n<-length(scores)

    Z<-sqrt(n-1)*(T0 - n*mean(scores)*mean(group))/
                  sqrt(SSE.scores*SSE.grp)
    p.lte<-pnorm(Z)
    p.gte<-1-pnorm(Z)
    p.twosided<-2*min(p.lte,p.gte)
    p.values<-c(p.twosided=p.twosided,p.twosidedAbs=p.twosided,p.lte=p.lte,p.gte=p.gte)
    out<-list(p.values=p.values,Z=Z)
    out
}


`trend.exact.mc` <-
function(scores,group,alternative="two.sided",nmc=10^3-1,seed=1234321,digits=12,p.conf.level=.99,setSEED=TRUE){
    if (!is.numeric(group)) stop("for trend tests group must be numeric")

    calcTestStat<-function(x,g){ sum( x*g) }

    t0<-calcTestStat(scores,group)
    N <- nmc
    if (setSEED) set.seed(seed)
    ti <- rep(NA, N)
    for (i in 1:N) {
        ti[i] <- calcTestStat(scores,sample(group,replace=FALSE))
    }   
    out<-calcPvalsMC(ti,t0,digits,alternative,FALSE,p.conf.level) 
    out
}


