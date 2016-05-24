Fclust <- function (X)
{
if (missing(X))
stop("The data set must be given")
if (is.null(X))
stop("The data set X is empty")
X=as.matrix(X)
if (any(is.na(X)))
stop("The data set X must not contain NA values")
if (!is.numeric(X)) 
stop("The data set X is not a numeric data.frame or matrix")
cat(" ",fill=TRUE)
cat("WELCOME to the interactive Fclust program",fill=TRUE)
cat("Warning: If you insert an object of mode CHARACTER when not requested, an error occurs and the program stops!")
cat(" ",fill=TRUE)
n=nrow(X)
startU=NULL
k=NULL
m=NULL
ent=NULL
vp=NULL
RS=NULL
conv=NULL
maxit=NULL
stand=NULL
delta=NULL
ana=1
ok=0
param=0
while (ana==1)
    {
    while (ok!=1)
        {
        cat(" ",fill=TRUE)
        cat("Specify the clustering algorithm: ",fill=TRUE)
        cat(" 1: Fuzzy k-means (function FKM) ",fill=TRUE)
        cat(" 2: Fuzzy k-means with noise cluster (function FKM.noise) ",fill=TRUE)
        cat(" 3: Fuzzy k-means with entropy regularization (function FKM.ent) ",fill=TRUE)
        cat(" 4: Fuzzy k-means with entropy regularization and noise cluster (function FKM.ent.noise) ",fill=TRUE)
        cat(" 5: Gustafson-Kessel-like Fuzzy k-means (function FKM.gk) ",fill=TRUE)
        cat(" 6: Gustafson-Kessel-like Fuzzy k-means with noise cluster (function FKM.gk.noise) ",fill=TRUE)
        cat(" 7: Gustafson-Kessel-like Fuzzy k-means with entropy regularization (function FKM.gk.ent) ",fill=TRUE)
        cat(" 8: Gustafson-Kessel-like Fuzzy k-means with entropy regularization and noise cluster (function FKM.gk.ent.noise) ",fill=TRUE)
        cat(" 9: Fuzzy k-medoids (function FKM.med) ",fill=TRUE)
        cat("10: Fuzzy k-medoids with noise cluster (function FKM.med.noise) ",fill=TRUE)
        cat("11: Fuzzy k-means with polynomial fuzzifier (function FKM.pf) ",fill=TRUE)
        cat("12: Fuzzy k-means with polynomial fuzzifier and noise cluster (function FKM.pf.noise) ",fill=TRUE)
        model=scan("",n=1) 
        if (any(model==1:12,na.rm=TRUE))
            ok=1
        else	
            {
            cat(" ",fill=TRUE)
            cat("Error! Make a proper choice for the clustering algorithm ",fill=TRUE)
            }
        } 
        if ((is.null(k)) || (param!=1))
            {
            cat(" ",fill=TRUE)
            cat("Specify the number of cluster k (default =2): ",fill=TRUE)
            k=scan("",n=1) 
            if (length(k)==0)
                {
                k=2
                }
            if ((k>ceiling(n/2)) || (k<2))
                {
                k=2
                cat("The number of clusters k must be an integer in {2, 3, ..., ceiling(n/2)}: the default value k=2 will be used ",fill=TRUE)
                }
            if (k%%ceiling(k)>0)  
                {
                k=ceiling(k)
                cat("The number of clusters k must be an integer in {2, 3, ..., ceiling(nrow(X)/2)}: the value ceiling(k) will be used ",fill=TRUE)
                }
            }
       if (any(model==c(1,2,5,6,9,10)))
            {
            if ((is.null(m)) || (param!=1))
                {
                cat(" ",fill=TRUE)
                cat("Specify the parameter of fuzziness m: (default =2)",fill=TRUE)
                m=scan("",n=1) 
                if (length(m)==0)
                    {
                    m=2
                    }
                if (m<=1) 
                    {
                    m=2
                    cat("The parameter of fuzziness m must be >1: the default value m=2 will be used ",fill=TRUE)
                    }
                }
            }
       if (any(model==c(3,4,7,8)))
            {
            if ((is.null(ent)) || (param!=1))
                {
                cat(" ",fill=TRUE)
                cat("Specify the degree of fuzzy entropy ent (default =1): ",fill=TRUE)
                ent=scan("",n=1) 
                if (length(ent)==0)
                    {
                    ent=1
                    }
                if (ent<=0)
                    {
                    ent=1
                    cat("The degree of fuzzy entropy ent must be >0: the default value ent=1 will be used ",fill=TRUE)
                    }
                }
            }
       if (any(model==c(11,12)))
            {
            if ((is.null(b)) || (param!=1))
                {
                cat(" ",fill=TRUE)
                cat("Specify the parameter of the polynomial fuzzifier beta (default =0.5): ",fill=TRUE)
				b=scan("",n=1) 
                if (length(b)==0)
                    {
                    b=0.5
                    }
                if ((b>1) || (b<0)) 
                    {
                    b=0.5
                    cat("The parameter of the polynomial fuzzifier beta must be in [0,1]: the default value beta=0.5 will be used ",fill=TRUE)
                    }
                }
            }
        if (any(model==c(3,4,5,6)))
            {
            if ((is.null(vp)) || (param!=1))
                {
                cat(" ",fill=TRUE)
                cat("Specify the volume parameter vp (default =rep(1,k)): ",fill=TRUE)
                vp=scan("",n=k) 
                if (length(vp)==0)
                    {
                    vp=rep(1,k)
                    }
                if (min(vp)<=0)
                    {
                    vp=rep(1,k)
                    cat("The volume parameter vp must be >0: the default value vp=rep(1,k) will be used ",fill=TRUE)
                    }
                }
            }
        if ((is.null(RS)) || (param!=1))
            {
            cat(" ",fill=TRUE)
            cat("Specify the number of random starts RS (default =1): ",fill=TRUE)
            RS=scan("",n=1) 
            if (length(RS)==0)
                {
                RS=1
                }
            if (RS<1)
                {
                cat("The number of starts RS must be an integer >=1: the default value RS=1 will be used ",fill=TRUE)
                RS=1
                } 
            if (RS%%ceiling(RS)>0)
                {
                cat("The number of starts RS  must be an integer >=1: the value ceiling(RS) will be used ",fill=TRUE)
                RS=ceiling(RS)
                }
            }
        if ((is.null(conv)) || (param!=1))
            {
            cat(" ",fill=TRUE)
            cat("Specify the convergence criterion conv (default =1e-9): ",fill=TRUE)
            conv=scan("",n=1) 
            if (length(conv)==0)
                {
                conv=1e-9
                }
            if (conv<=0) 
                {
                cat("The convergence criterion conv must be a (small) value >0: the default value conv=1e-9 will be used ",fill=TRUE)
                conv=1e-9
                }
            }
        if ((is.null(maxit)) || (param!=1))
            {
            cat(" ",fill=TRUE)
            cat("Specify the maximum number of iterations maxit (default =1e+6): ",fill=TRUE)
            maxit=scan("",n=1) 
            if (length(maxit)==0)
                {
                maxit=1e+6
                }
            if (maxit<=0)
                {
                cat("The maximum number of iterations maxit must be an integer >0: the default value maxit=1e+6 will be used ",fill=TRUE)
                maxit=1e+6
                } 
            if (maxit%%ceiling(maxit)>0)
                {
                cat("The maximum number of iterations maxit must be an integer >0: the value ceiling(maxit) will be used ",fill=TRUE)
                maxit=1e+6
                }
            }
        if ((is.null(stand)) || (param!=1))
            {
            cat(" ",fill=TRUE)
            cat("If you want to standardize the dataset before running the clustering algorithm, specify 1: ",fill=TRUE)
            stand=scan("",n=1) 
            if (length(stand)==0)
                {
                stand=0
                }
            }
        if (model==c(2,4,6,8,10,12))
            {
            if ((is.null(delta)) || (param!=1))
                {
                cat(" ",fill=TRUE)
                cat("Specify the noise distance delta (see ?FKM.noise for the default value): ",fill=TRUE)
                delta=scan("",n=1) 
                if (length(delta)==0)
                    {
                    Dd=matrix(0,nrow=n,ncol=k) 
                    Hd=FKM(X,k,m,RS=1,stand,conv=1e-6)$H
                    for (i in 1:n) 
                        {
                        for (c in 1:k) 
                            {
                            Dd[i,c]=sum((X[i,]-Hd[c,])^2)
                            }
                        }
                    delta=mean(Dd)
                    }
                if (delta<0)
                    {
                    cat("The noise distance delta must be non negative: the default value (see ?FKM.noise) will be used ",fill=TRUE)
                    Dd=matrix(0,nrow=n,ncol=k) 
                    Hd=FKM(X,k,m,RS=1,stand,conv=1e-6)$H
                    for (i in 1:n) 
                        {
                       for (c in 1:k) 
                            {
                            Dd[i,c]=sum((X[i,]-Hd[c,])^2)
                            }
                        }
                    delta=mean(Dd)
                    } 
                }
            }
        if (model==1)
           {
            clust=FKM(X,k,m,RS,stand,startU,conv,maxit)
            }
        if (model==2)
            {
            clust=FKM.noise(X,k,m,delta,RS,stand,startU,conv,maxit)
            }
        if (model==3)
            {
            clust=FKM.ent(X,k,ent,RS,stand,startU,conv,maxit)
            }
        if (model==4)
            {
          clust=FKM.ent.noise(X,k,ent,delta,RS,stand,startU,conv,maxit)
            }
        if (model==5)
            {
            clust=FKM.gk(X,k,m,vp,RS,stand,startU,conv,maxit)
            }
        if (model==6)
            {
          clust=FKM.gk.noise(X,k,m,vp,delta,RS,stand,startU,conv,maxit)
            }
        if (model==7)
        {
            clust=FKM.gk.ent(X,k,ent,vp,RS,stand,startU,conv,maxit)
        }
        if (model==8)
        {
          clust=FKM.gk.ent.noise(X,k,ent,vp,delta,RS,stand,startU,conv,maxit)
        }
        if (model==9)
        {
          clust=FKM.med(X,k,m,RS,stand,startU,conv,maxit)
        }
        if (model==10)
        {
          clust=FKM.med.noise(X,k,m,delta,RS,stand,startU,conv,maxit)
        }
        if (model==11)
        {
		   clust=FKM.pf(X,k,b,RS,stand,startU,conv,maxit)
        }
        if (model==12)
        {
          clust=FKM.pf.noise(X,k,b,delta,RS,stand,startU,conv,maxit)
        }
        cat(" ",fill=TRUE)
        cat("Some preliminary results are displayed: ",fill=TRUE)
        cat(" ",fill=TRUE)
        cat("Membership degree matrix ",fill=TRUE)
        print(round(clust$U,digits=2))
        cat(" ",fill=TRUE)
        cat("Prototype matrix ",fill=TRUE)
        print(round(clust$H,digits=2))
        cat(" ",fill=TRUE)
        cat("If you want to perform a new cluster analysis, specify 1: ",fill=TRUE)
        ana=scan("",n=1) 
        if (length(ana)==0)
            {
            ana=0
            }
        if (ana==1)
            {
            cat(" ",fill=TRUE)
            cat("If you want to run the same clustering algorithm, specify 1: ",fill=TRUE)
            ok=scan("",n=1) 
            if (length(ok)==0)
                {
                ok=0
                }
            if (ok!=1)
                {
                cat(" ",fill=TRUE)
                cat("You chose to run a new clustering algorithm",fill=TRUE)
                cat("If you want to use the same parameter values, specify 1: ",fill=TRUE)
                param=scan("",n=1) 
                if (length(param)==0)
                    {
                    param=0
                    }
                }
            else
                param=0
            if (param==0)
                {
                startU=NULL
                k=NULL
                m=NULL
                ent=NULL
                vp=NULL
                RS=NULL
                conv=NULL
                maxit=NULL
                stand=NULL
                delta=NULL
                }
            }
        }
    cat(" ",fill=TRUE)
    cat("The results of the last analysis are displayed or saved ",fill=TRUE)
    return(clust)
}