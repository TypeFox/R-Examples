##  From the R-help e-mail by Ted Harding: http://tolstoy.newcastle.edu.au/R/e2/help/07/03/12853.html
##  See also http://tolstoy.newcastle.edu.au/R/help/05/05/4254.html
pava <- function(x, wt = rep(1, length(x)))
{
    n <- length(x)
    if (n <= 1) return(x)
    lvlsets <- 1:n
    repeat 
    {
        viol <- (as.vector(diff(x)) < 0)
        if (!(any(viol))) break
        i <- min((1:(n-1))[viol])
    
        lvl1 <- lvlsets[i]
        lvl2 <- lvlsets[i+1]
        ilvl <- ( (lvlsets == lvl1) | (lvlsets == lvl2) )
    
        x[ilvl] <- sum(x[ilvl] * wt[ilvl]) / sum(wt[ilvl])     
        lvlsets[ilvl] <- lvl1
    }
    x
} 


# from Cai Tianxi
## count how many YY's are smaller or equal to yy
N.L.E <- function(yy, YY)  ## sum I(YY <= yy[i])
{
   rank(c(yy+1e-8,YY))[1:length(yy)] - rank(yy)  ### add a small pertubation to avoid ties when calculating rank
}
N.L <- function(yy, YY)  ## sum I(YY < yy[i])
{
   rank(c(yy-1e-8,YY))[1:length(yy)] - rank(yy)  ### add a small pertubation to avoid ties when calculating rank
}
N.G.E <- function(yy, YY)  ## sum I(YY >= yy[i])
{
   length(YY)-(rank(c(yy-1e-8,YY))[1:length(yy)] - rank(yy))
}

# get first prinipal component
pr.1=function(x){
    x.s=scale(x)
    pr.s = prcomp(x.s)
    out=c(x.s %*% pr.s$rotation[,1])
    attr(out, "rotation")=pr.s$rotation[,1]
    out
}

# compute covariability as defined in Clarke 1995 page 2271, second formula
# x is a vector of characters, so is y
covariability=function(x,y){
    tmp=table(x)
    p.x=data.frame(tmp/(sum(tmp)))
    tmp=table(y)
    p.y=data.frame(tmp/(sum(tmp)))
    
    # tabulate all pairwise combination of aa
    pair = paste(x,y,sep="")
    tmp=table(pair)
    dat=data.frame(tmp/(sum(tmp)))
    
    ps=numeric(nrow(dat))
    for (i in 1:nrow(dat)) {
        row.=dat[i,]
        a2=as.character(row.[1,1])
        p.i=p.x$Freq[p.x$x==substr(a2,1,1)]
        p.j=p.y$Freq[p.y$y==substr(a2,2,2)]
        p.ij=row.[1,2]
        ps[i]=p.ij**2*log(p.ij/p.i/p.j)        
    }
    sum(ps)
}
## test
#dat=readFastaFile ("D:/1CHMM/domains/SET/SETpfamseed/seq/SETpfamseed_aligned.fasta")
#seqs.mat.a=mysapply(dat, s2c)
#covar=matrix(0,ncol(seqs.mat.a),ncol(seqs.mat.a))
#for (i in 1:(ncol(seqs.mat.a)-1)){
#    for (j in (i+1):ncol(seqs.mat.a)){
#        covar[i,j]<-covar[j,i]<-covariability(seqs.mat.a[,i], seqs.mat.a[,j])
#    }
#}
#sort(round(covar,3), decre=TRUE)[1:100]

# output format:
#+1 1:0.708333 2:1 3:1 4:-0.320755 5:-0.105023 6:-1 7:1 8:-0.419847 9:-1 10:-0.225806 12:1 13:-1 
# y: vector of 1 and non-1
# z: matrix or data.frame of covariates
# file.name: e.g. test. There should be no file extension
# ws: vector of weights
write.svm.input=function(y, z, file.name="svm_input", ws=NULL){
    rows=sapply(1:length(y), function(i){
        paste(c(ifelse(y[i]==1,"+1","-1"),(1:ncol(z))%+%":"%+%z[i,]), collapse=" ")
    })
    write(paste(rows,collapse="\n"), file = file.name%+%".dat", append = FALSE, sep = " ")        
    # write weight file
    if (!is.null(ws)) write(paste(ws,collapse="\n"), file = file.name%+%".wgt", append = FALSE, sep = " ")    
}

# m is data frame or matrix, each row is an observation
# w is a named vector, the names may not match m, may not have same number even
# return a vector
weighted.ave=function(m, w){
    if(is.data.frame(m)){
        nam=names(m)
        m=as.matrix(m)
    } else {
        nam=colnames(m)
    }
    ret=w[match(nam, names(w))]
    ret[is.na(ret)]=0
    names(ret)=nam
    c(m %*% ret )
}

# st.fit is an object returned by sn::st.mle, skew normal
get.st.mean=function(st.fit) {
    p=st.fit$dp
    
    xi=p["location"];
    alpha=p["shape"];
    df=p["df"]
    #omega=sd(X);
    omega=p["scale"]
    delta=alpha/sqrt(1+alpha^2);
    mu=delta*sqrt(df/pi)*gamma(0.5*(df-1))/gamma(0.5*df);
    
    # the 4 first moments
    moment1=xi+omega*mu;
    moment2=xi^2 + 2*omega*mu*xi + omega^2 * df/(df-2);
    moment3=xi^3 + 3*omega*mu*xi^2 + 3*omega^2*df/(df-2)*xi +
    omega^3*mu*(3-delta^2)*df/(df-3);
    moment4=xi^4 + 4*omega*mu*xi^3 + 6*omega^2*df/(df-2)*xi^2 +
    4*omega^3*mu*(3-delta^2)*df/(df-3)*xi+omega^4*3*df^2/((df-2)*(df-4));
    
    # the 4 useful measures
    mean=moment1;
    var=moment2;
    skew=moment3/var^(3/2);
    kurt=moment4/var^2 - 3;
    
    c(mean)
}


# returns sandwich estimator of variance matrix
# from Thomas Lumley
infjack.glm<-function(glm.obj,groups){
    umat<-estfun.glm(glm.obj)
    usum<-rowsum(umat,groups,reorder=FALSE)
    modelv<-summary(glm.obj)$cov.unscaled
    modelv%*%(t(usum)%*%usum)%*%modelv
}
jack.glm<-function(glm.obj,groups){
    umat<-jackvalues(glm.obj)
    usum<-rowsum(umat,groups,reorder=FALSE)
    t(usum)%*%usum*(nrow(umat)-1)/nrow(umat)
}
jackvalues<-function(glm.obj){
    db<-lm.influence(glm.obj)$coef
    t(t(db)-apply(db,2,mean))
}   
estfun.glm<-function(glm.obj){
    if (is.matrix(glm.obj$x)) 
        xmat<-glm.obj$x
    else {
        mf<-model.frame(glm.obj)
        xmat<-model.matrix(terms(glm.obj),mf)       
    }
    residuals(glm.obj,"working")*glm.obj$weights*xmat
}


#iterated sum, like diff
summ=function(x) {
    x[-1]+x[-length(x)]
}

list_args <- Vectorize( function(a,b) c( as.list(a), as.list(b) ), SIMPLIFY = FALSE)
make_args_mtx <- function( alist ) Reduce(function(x, y) outer(x, y, list_args), alist)
multi.outer <- function(f, ... ) {
  args <- make_args_mtx(list(...))
  apply(args, 1:length(dim(args)), function(a) do.call(f, a[[1]] ) )
}


get.sim.res = function(foldername, verbose=FALSE) {
    
    if (!requireNamespace("abind")) {print("abind does not load successfully"); return (NULL) }
    
    ## get file names
    tmp=list.files(path=foldername, pattern="batch[0-9]+.*.Rdata")
    if (length(tmp)==0) {
        # no files yet
        stop ("no files\n")
    } else {
        fileNames = foldername %+%"/"%+%tmp
    }
    cat("number of files: "%+%length(fileNames),"\n")
    
    ## read files
    res=lapply(fileNames, function(x) {load(file=x); res})
    res.all = do.call(abind::abind, res)
    names(dimnames(res.all))=names(dimnames(res[[1]]))
    
    if(verbose) {
        cat("res.all:\n")
        print(str(res.all))        
        cat("\nFirst dimension col names: \n")
        print(dimnames(res.all)[[1]])
    }
    
    res.all    
}
