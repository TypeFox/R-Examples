NPtest<-function(obj, n=NULL, method="T1", ...){
#npt<-function(obj, n=NULL, method="T1", ...){
#-------------------------------------------------
# changes: 2011-12-06
#-------------------------------------------------
# in ... koennen ausser den spezifikationen fuer die
# einzelnen statistiken wie z.B. idx=, stat=, etc.
# nun zusaetzlich step=, burn_in=, und seed=
# angegeben werden
# ausserdem mit RSinfo=TRUE wird eine summary des
# RaschSampler output objects geprintet
#-------------------------------------------------
# Aenderung im Code:
# - alle methods bei switch haben jetzt ... in der
#   argument liste
# - check des input objekts
#-------------------------------------------------


#   require(RaschSampler)   # removed as RaschSampler is a part of eRm since 0.15-0

   dots<-as.list(substitute(list(...)))[-1]
   nn<-names(dots)
   for (i in seq(along=dots)) assign(nn[i],dots[[i]])

   if(!exists("burn_in", inherits = FALSE)) burn_in <- 256
   if(!("step" %in% nn)) step<-32
   if(!exists("seed", inherits = FALSE)) seed<-0
   if(is.null(n)) n <- 500

   if(is.matrix(obj) || is.data.frame(obj)){ # input is datamatrix -  RaschSampler object is generated
      if (!all(obj %in% 0:1)) stop("Data matrix must be binary, NAs not allowed")
      itscor<-colSums(obj) # rh 2011-03-03
      itcol<-(itscor==0|itscor==nrow(obj))
      if (any(itcol)){
        cat("The following columns in the data show complete 0/full responses: \n")
        cat((1:ncol(obj))[itcol],sep=", ")
        cat("\n")
        stop("NPtest using these items is meaningless. Delete them first!")
      }
      obj<-rsampler(obj,rsctrl(burn_in=burn_in, n_eff=n, step=step, seed=seed))
#browser()

   } else if(class(obj)!="RSmpl"){
        stop("Input object must be data matrix/data frame or output from RaschSampler")
   }

   if(exists("RSinfo", inherits = FALSE)) if(get("RSinfo")) summary(obj)

  switch(method,
    "T1"    = T1(obj, ...),
    "T1l"   = T1l(obj, ...),
    "T1m"   = T1m(obj, ...),
    "Tmd"   = Tmd(obj, ...),
    "T2"    = T2(obj, ...),
    "T2m"   = T2m(obj, ...),
    "T4"    = T4(obj, ...),
#    "T7"    = T7(obj, ...),
#    "T7a"   = T7a(obj, ...),
    "T10"   = T10(obj, ...),
    "T11"   = T11(obj, ...),
    "Tpbis" = Tpbis(obj, ...),
    "MLoef" = MLoef.x(obj, ...)
  )
}

MLoef.x<-function(rsobj, splitcr=NULL, ...){
     # user function
     MLexact<-function(X,splitcr){
       rmod<-RM(X)
       LR<-MLoef(rmod,splitcr)$LR
       LR
     }
     #if(!exists("splitcr", inherits = FALSE)) splitcr="median"
     if(is.null(splitcr)) splitcr="median"
     res <- rstats(rsextrobj(rsobj, 2), MLexact, splitcr)

     rmod<-RM(rsextrmat(rsobj,1))                     # MLoef for original data
     MLres<-MLoef(rmod,splitcr)
     class(MLres)<-c(class(MLres),"MLx")              # for printing without blank line
     res1<-MLres$LR

     n_eff<-rsobj$n_eff                         # number of simulated matrices
     res<-unlist(res)
     prop<-sum((res[1:n_eff]>=res1)/n_eff)

     result<-list(MLres=MLres, n_eff=n_eff, prop=prop, MLoefvec=res) # MLobj
     class(result)<-"MLobj"
     result
}

Tpbis <- function(rsobj, idxt=NULL, idxs=NULL, ...){ # fixed 2013-08-09
  Tpbis.stat <- function(x){
    rb <- rowSums(x[, idxs, drop = FALSE])     # all raw scores
    t  <- x[, idxt]                            # dichotomous item
    r  <- tapply(rb, t, sum, simplify = FALSE) # raw scores by item; simplify = FALSE to be on the safe side
    n1 <- sum(t)                               # n_1 = sum of raw scores with t == 1
    n0 <- sum(1 - t)                           # n_0 = sum of raw scores with t == 0
    return(n0 * r[[2L]][1L] - n1*r[[1L]][1L])  # n_0 * sum(r_1) - n_1 * sum(r_0)
  }

  if(is.null(idxs)) stop("No item(s) for subscale  specified (use idxs!)")
  if(is.null(idxt)) stop("No test item for testing against subscale specified (use idx!)")
  li1 <- length(idxt)
  li2 <- length(idxs)
  k   <- rsobj$k
  if(li1 > 1L ||li2 >= k || (li1 + li2) > k || any(idxt %in% idxs) || any(c(idxt,idxs) > k)){
    stop("Subscale and/or test item incorrectly specified.")
  }

  n_eff <- rsobj$n_eff                   # number of simulated matrices
  n_tot <- rsobj$n_tot                   # number of simulated matrices

  res     <- rstats(rsobj, Tpbis.stat)              # calculates statistic for each matrix
  corrvec <- do.call(cbind, lapply(res, as.vector)) # converts result list to matrix

  prop <- sum(corrvec[2L:(n_tot)] <= corrvec[1L]) / n_eff   # T(A_s) >= T(A_0)

  # Tpbisobj
  result <- list("n_eff"    = n_eff,
                 "prop"     = prop,
                 "idxt"     = idxt,
                 "idxs"     = idxs,
                 "Tpbisvec" = corrvec)
  class(result)<-"Tpbisobj"
  return(result)
}

Tmd<-function(rsobj, idx1=NULL, idx2=NULL, ...){
     Tmd.stat<-function(x){
        r1<-rowSums(x[,idx1, drop=FALSE])
        r2<-rowSums(x[,idx2, drop=FALSE])
        corr<-cor(r1,r2)
        corr
     }

     if(is.null(idx1))
         stop("No item(s) for subscale 1 specified (use idx1!)")
     if(is.null(idx2))
         stop("No item(s) for subscale 2 specified (use idx2!)")
     li1<-length(idx1)
     li2<-length(idx2)
     k<-rsobj$k
     if(li1>=k ||li2>=k || li1+li2>k || any(idx1 %in% idx2))
         stop("Subscale(s) incorrectly specified.")

     n_eff<-rsobj$n_eff                         # number of simulated matrices
     n_tot<-rsobj$n_tot                         # number of simulated matrices

     res<-rstats(rsobj,Tmd.stat)               # calculates statistic for each matrix
     corrvec<-do.call(cbind, lapply(res,as.vector)) # converts result list to matrix

     prop<-sum(corrvec[2:(n_tot)]<=corrvec[1])/n_eff

     result<-list(n_eff=n_eff, prop=prop, idx1=idx1, idx2=idx2, Tmdvec=corrvec)   # Tmdobj
     class(result)<-"Tmdobj"
     result
}


T1m<-function(rsobj, ...){
     T1mstat<-function(x){      # calculates statistic T1m
        unlist(lapply(1:(k-1),function(i) lapply((i+1):k, function(j) sum(x[,i]==x[,j]))))
     }
     n_eff<-rsobj$n_eff                         # number of simulated matrices
     n_tot<-rsobj$n_tot                         # number of simulated matrices
     k<-rsobj$k                                 # number of columns of matrices

     res<-rstats(rsobj,T1mstat)                  # calculates statistic for each matrix

     res<-do.call(cbind, lapply(res,as.vector)) # converts result list to matrix
     T1mvec<-apply(res, 1, function(x) sum(x[2:(n_tot)]<=x[1])/n_eff)
     T1mmat<-matrix(,k,k)
     T1mmat[lower.tri(T1mmat)] <- T1mvec           # lower triangular matrix of p-values
     result<-list(n_eff=n_eff, prop=T1mvec, T1mmat=T1mmat) # T1mobj
     class(result)<-"T1mobj"
     result
}

T1<-function(rsobj, ...){
     T1stat<-function(x){      # calculates statistic T1
        unlist(lapply(1:(k-1),function(i) lapply((i+1):k, function(j) sum(x[,i]==x[,j]))))
     }
     n_eff<-rsobj$n_eff                         # number of simulated matrices
     n_tot<-rsobj$n_tot                         # number of simulated matrices
     k<-rsobj$k                                 # number of columns of matrices

     res<-rstats(rsobj,T1stat)                  # calculates statistic for each matrix

     res<-do.call(cbind, lapply(res,as.vector)) # converts result list to matrix
     T1vec<-apply(res, 1, function(x) sum(x[2:(n_tot)]>=x[1])/n_eff)
     T1mat<-matrix(,k,k)
     T1mat[lower.tri(T1mat)] <- T1vec           # lower triangular matrix of p-values
     result<-list(n_eff=n_eff, prop=T1vec, T1mat=T1mat) # T1obj
     class(result)<-"T1obj"
     result
}

T1l<-function(rsobj, ...){
     T1lstat<-function(x){      # calculates statistic T1
        unlist(lapply(1:(k-1),function(i) lapply((i+1):k, function(j) sum(x[,i] & x[,j]))))
     }
     n_eff<-rsobj$n_eff                         # number of simulated matrices
     n_tot<-rsobj$n_tot                         # number of simulated matrices
     k<-rsobj$k                                 # number of columns of matrices

     res<-rstats(rsobj,T1lstat)                  # calculates statistic for each matrix

     res<-do.call(cbind, lapply(res,as.vector)) # converts result list to matrix
     T1lvec<-apply(res, 1, function(x) sum(x[2:(n_tot)]>=x[1])/n_eff)
     T1lmat<-matrix(,k,k)
     T1lmat[lower.tri(T1lmat)] <- T1lvec           # lower triangular matrix of p-values
     result<-list(n_eff=n_eff, prop=T1lvec, T1lmat=T1lmat) # T1obj
     class(result)<-"T1lobj"
     result
}
T2<-function(rsobj,idx=NULL,stat="var", ...){

     T2.Var.stat<-function(x){       # calculates statistic T2
        var(rowSums(x[,idx, drop=FALSE]))
     }
     T2.MAD1.stat<-function(x){       # calculates statistic T2
        y<-rowSums(x[,idx, drop=FALSE])           # mean absolute deviation
        mean(abs(y-mean(y)))
     }
     T2.MAD2.stat<-function(x){       # calculates statistic T2
        mad(rowSums(x[,idx, drop=FALSE]),constant=1) # unscaled median absolute deviation
     }
     T2.Range.stat<-function(x){     # calculates statistic T2
        diff(range(rowSums(x[,idx, drop=FALSE])))
     }
     n<-rsobj$n
     n_eff<-rsobj$n_eff
     k<-rsobj$k                      # number of columns of matrices
     if(is.null(idx))
         stop("No item(s) for subscale specified (use idx!)")
     res<-switch(stat,
          "var"=rstats(rsobj,T2.Var.stat),
          "mad1"=rstats(rsobj,T2.MAD1.stat),
          "mad2"=rstats(rsobj,T2.MAD2.stat),
          "range"=rstats(rsobj,T2.Range.stat),
          stop("stat must be one of \"var\", \"mad1\", \"mad2\", \"range\"")
     )
     res<-unlist(res)
     prop<-sum(res[2:(n_eff+1)]>=res[1])/n_eff
     result<-list(n_eff=n_eff, prop=prop, idx=idx, stat=stat, T2vec=res) # T2obj
     class(result)<-"T2obj"
     result
}

T2m<-function(rsobj,idx=NULL,stat="var", ...){

     T2m.Var.stat<-function(x){       # calculates statistic T2m
        var(rowSums(x[,idx, drop=FALSE]))
     }
     T2m.MAD1.stat<-function(x){       # calculates statistic T2m
        y<-rowSums(x[,idx, drop=FALSE])           # mean absolute deviation
        mean(abs(y-mean(y)))
     }
     T2m.MAD2.stat<-function(x){       # calculates statistic T2m
        mad(rowSums(x[,idx, drop=FALSE]),constant=1) # unscaled median absolute deviation
     }
     T2m.Range.stat<-function(x){     # calculates statistic T2m
        diff(range(rowSums(x[,idx, drop=FALSE])))
     }
     n<-rsobj$n
     n_eff<-rsobj$n_eff
     k<-rsobj$k                      # number of columns of matrices
     if(is.null(idx))
         stop("No item(s) for subscale specified (use idx!)")
     res<-switch(stat,
          "var"=rstats(rsobj,T2m.Var.stat),
          "mad1"=rstats(rsobj,T2m.MAD1.stat),
          "mad2"=rstats(rsobj,T2m.MAD2.stat),
          "range"=rstats(rsobj,T2m.Range.stat),
          stop("stat must be one of \"var\", \"mad1\", \"mad2\", \"range\"")
     )
     res<-unlist(res)
     prop<-sum(res[2:(n_eff+1)]<=res[1])/n_eff
     result<-list(n_eff=n_eff, prop=prop, idx=idx, stat=stat, T2mvec=res) # T2mobj
     class(result)<-"T2mobj"
     result
}


T4<-function(rsobj,idx=NULL,group=NULL,alternative="high", ...){

     T4.stat<-function(x){      # calculates statistic T4
        sign*sum(rowSums(x[gr,idx,drop=FALSE]))
     }
     n_eff<-rsobj$n_eff                         # number of simulated matrices
     n_tot<-rsobj$n_tot                         # number of all matrices
     k<-rsobj$k                                 # number of items
     if(is.null(idx))
         stop("No item(s) for subscale specified (use idx!)")
     if(length(idx)==k)  # rh 2011-03-03
         stop("Subscale containing all items gives meaningless results for T4.")
     if(is.null(group))
         stop("No group specified (use group!)")
     if(!is.logical(group))     # added rh 2011-03-03
         stop("group must be of type \"logical\" (e.g., group = (age==1) )")
     if(alternative=="high")
        sign <- 1
     else if(alternative=="low")
        sign <- -1
     else
        stop("alternative incorrectly specified! (use either \"high\" or \"low\")")

     gr<-as.logical(group)                      # group definition (logical)
     res<-rstats(rsobj,T4.stat)
     res<-unlist(res)
     prop<-sum(res[2:(n_tot)]>=res[1])/n_eff
     gr.nam <- deparse(substitute(group))
     gr.n <- sum(group)
     result<-list(n_eff=n_eff, prop=prop, idx=idx, gr.nam=gr.nam, gr.n=gr.n, T4vec=res, alternative=alternative)   # T4obj
     class(result)<-"T4obj"
     result
}
# removed in version 0.14-5
#T7<-function(rsobj,idx=NULL, ...){
#     T7.stat<-function(x){      # calculates statistic T7
#        calcT7<-function(i,j){  # calculates sum for all items in subscale
#          if(sitscor[i]>sitscor[j]){
#              sum(submat[,j]>submat[,i])   #
#              # t<-table(submat[,i],submat[,j])    # odds ratio gives the same result
#              # OR<-t[1]*t[4]/(t[2]*t[3])
#              # 1/OR
#          } else
#              NA
#        }
#        submat<-x[,idx]
#        submat<-submat[,order(itscor,decreasing=TRUE)]
#        RET<-unlist(lapply(1:(m-1), function(i) lapply((i+1):m, function(j) calcT7(i,j))))
#        RET
#     }
#
#     n_eff<-rsobj$n_eff                         # number of simulated matrices
#     n_tot<-rsobj$n_tot                         # number of all matrices
#     k<-rsobj$k                                 # number of items
#     if(is.null(idx))
#         stop("No items for subscale specified (use idx!)")
#     else if (length(idx)<2)
#         stop("At least 2 items have to be specified with idx!")
#     submat<-rsextrmat(rsobj,1)[,idx]
#     itscor<-colSums(submat)
#     names(itscor)<-colnames(submat)<-idx
#
#     submat<-submat[,order(itscor,decreasing=TRUE)]
#     sitscor<-sort(itscor,decreasing=TRUE)      # sorted itemscore
#     m<-length(itscor)
#
#     resList<-rstats(rsobj,T7.stat)
#     res<-sapply(resList,sum,na.rm=TRUE)
#     prop<-sum(res[2:(n_eff+1)]>=res[1])/n_eff
#     result<-list(n_eff=n_eff, prop=prop, itscor=itscor, T7vec=res)   # T7obj
#     class(result)<-"T7obj"
#     result
#}
#T7a<-function(rsobj,idx=NULL, ...){
#     T7a.stat<-function(x){      # calculates statistic T7a
#        calcT7a<-function(i,j){  # calculates sum for single Itempair
#          if(sitscor[i]>sitscor[j]){
#              sum(submat[,j]>submat[,i])   #
#              # t<-table(submat[,i],submat[,j])    # odds ratio gives the same result
#              # OR<-t[1]*t[4]/(t[2]*t[3])
#              # 1/OR
#          } else
#              NA
#        }
#        submat<-x[,idx]
#        submat<-submat[,order(itscor,decreasing=TRUE)]
#        RET<-unlist(lapply(1:(m-1), function(i) lapply((i+1):m, function(j) calcT7a(i,j))))
#        RET
#     }
#
#     n_eff<-rsobj$n_eff                         # number of simulated matrices
#     n_tot<-rsobj$n_tot                         # number of all matrices
#     k<-rsobj$k                                 # number of items
#     if(is.null(idx))
#         stop("No items for subscale specified (use idx!)")
#     else if (length(idx)<2)
#         stop("At least 2 items have to be specified with idx!")
#     submat<-rsextrmat(rsobj,1)[,idx]
#     itscor<-colSums(submat)
#     names(itscor)<-colnames(submat)<-idx
#     submat<-submat[,order(itscor,decreasing=TRUE)]
#     sitscor<-sort(itscor,decreasing=TRUE)      # sorted itemscore
#     m<-length(itscor)
#
#     res<-rstats(rsobj,T7a.stat)
#     res<-do.call(cbind, lapply(res,as.vector)) # converts result list to matrix
#     T7avec<-apply(res, 1, function(x) sum(x[2:(n_tot)]>=x[1])/n_eff)
#     T7anam<-NULL
#     for (i in 1:(m-1)) for(j in (i+1):m )
#          T7anam<-c(T7anam, paste("(",names(sitscor[i]),">",names(sitscor[j]),")",sep="",collapse=""))
#     names(T7avec)<-T7anam
#     result<-list(n_eff=n_eff, prop=T7avec,itscor=itscor)    # T7aobj
#     class(result)<-"T7aobj"
#     result
#}

T10<-function(rsobj, splitcr="median", ...){
      calc.groups<-function(x,splitcr){
        if (length(splitcr) > 1)  {        # numeric vectors converted to factors
            if (length(splitcr) != nrow(x)) {
                stop("Mismatch between length of split vector and number of persons!")
            }
            splitcr <- as.factor(splitcr)
            if (length(levels(splitcr))>2) {
                stop("Split vector defines more than 2 groups (only two allowed)!")
            }
            spl.lev <- levels(splitcr)
            #spl.gr <- paste(spl.nam, spl.lev, sep = " ")  # not necessary for the time being
            hi <- splitcr==spl.lev[1] # first level is high group
        } else if (!is.numeric(splitcr)) {
            spl.nam <- splitcr
            if (splitcr == "median") {
                spl.gr <- c("Raw Scores <= Median", "Raw Scores > Median")
                rv <- rowSums(x)
                rvsplit <- median(rv)
                hi <- rv > rvsplit
            }
            if (splitcr == "mean") {
                spl.gr <- c("Raw Scores < Mean", "Raw Scores >= Mean")
                rv <- rowSums(x)
                rvsplit <- mean(rv)
                hi <- rv > rvsplit
            }
        }
        list(hi=hi,spl.nam=spl.nam) # spl.nam is returned due to lex scoping even if not defined here
      }
      T10.stat<-function(x){      # calculates statistic T10 for one matrix
        nij.hi<-unlist(lapply(1:k,function(i) lapply(1:k, function(j) sum(x[hi,i]>x[hi,j]))))
        nij.low<-unlist(lapply(1:k,function(i) lapply(1:k, function(j) sum(x[!hi,i]>x[!hi,j]))))
        nji.hi<- unlist(lapply(1:k,function(i) lapply(1:k, function(j) sum(x[hi,i]<x[hi,j]))))
        nji.low<- unlist(lapply(1:k,function(i) lapply(1:k, function(j) sum(x[!hi,i]<x[!hi,j]))))
        RET<-sum(abs(nij.hi*nji.low-nij.low*nji.hi))
        RET
      }
      spl.nam <- deparse(substitute(splitcr))
      n_eff<-rsobj$n_eff                         # number of simulated matrices
      n_tot<-rsobj$n_tot                         # number of all matrices
      k<-rsobj$k                                 # number of columns of matrices
      obj<-rsextrobj(rsobj,1,1)                  # extract first matrix
      x<-matrix(obj$inpmat,obj$n,obj$k)
      ans <- calc.groups(x,splitcr)      # calculate grouping vector (logical)
      hi<-ans$hi
      hi.n<-sum(hi)
      low.n<-sum(!hi)

      res<-rstats(rsobj,T10.stat)                # for each matrix calculate T10

      res<-unlist(res)
      prop<-sum(res[2:(n_eff+1)]>=res[1])/n_eff
      result<-list(n_eff=n_eff, prop=prop,spl.nam=ans$spl.nam,hi.n=hi.n,low.n=low.n,T10vec=res)  # T10obj
      class(result)<-"T10obj"
      result
}


T11<-function(rsobj, ...){
      T11.stat<-function(x){
         as.vector(cor(x))
      }
      calc.T11<-function(x){      # calculates statistic T11 for one matrix
         sum(abs(x-rho))
      }
      n_eff<-rsobj$n_eff                         # number of simulated matrices
      n_tot<-rsobj$n_tot                         # number of all matrices
      k<-rsobj$k                                 # number of columns of matrices
      res<-rstats(rsobj,T11.stat)                # for each matrix calculate all r_ij's

      cormats <- matrix(unlist(res),nrow=k*k)    # k*k x n_tot matrix, each colum contains one corr matrix
      rho<-apply(cormats[,2:n_tot],1,mean)       # vector of estimated "real" rho_ij's
      T11obs<-calc.T11(cormats[,1])              # vector of observed r_ij's
      prop<-sum(apply(cormats[, 2:n_tot],2,calc.T11)>=T11obs)/n_eff
      result<-list(n_eff=n_eff, prop=prop, T11r=cormats[,1], T11rho=rho)   # T11obj
      class(result)<-"T11obj"
      result
}

print.MLobj<-function(x,...){
  print(x$MLres)
  cat("'exact' p-value =", x$prop, " (based on", x$n_eff, "sampled matrices)\n\n")
}

print.Tmdobj<-function(x,...){
  txt1<-"\nNonparametric RM model test: Tmd (Multidimensionality)"
  writeLines(strwrap(txt1, exdent=4))
  cat("    (correlation of subscale person scores)\n")
  cat("Number of sampled matrices:", x$n_eff,"\n")
  cat("Subscale 1 - Items:", x$idx1,"\n")
  cat("Subscale 2 - Items:", x$idx2,"\n")
  cat("Observed correlation:", x$Tmdvec[1],"\n")
  cat("one-sided p-value:",x$prop,"\n\n")
}

print.Tpbisobj<-function(x,...){
  txt1<-"\nNonparametric RM model test: Tpbis (discrimination)"
  writeLines(strwrap(txt1, exdent=4))
  cat("    (pointbiserial correlation of test item vs. subscale)\n")
  cat("Number of sampled matrices:", x$n_eff,"\n")
  cat("Test Item:", x$idxt,"\n")
  cat("Subscale  - Items:", x$idxs,"\n")
  cat("one-sided p-value (rpbis too low):",x$prop,"\n\n")
}

print.T1obj<-function(x,alpha=0.05,...){
  txt1<-"\nNonparametric RM model test: T1 (local dependence - increased inter-item correlations)\n"
  writeLines(strwrap(txt1, exdent=4))
  cat("    (counting cases with equal responses on both items)\n")
  cat("Number of sampled matrices:", x$n_eff,"\n")
  cat("Number of Item-Pairs tested:", length(x$prop),"\n")
  cat("Item-Pairs with one-sided p <", alpha,"\n")
  T1mat<-x$T1mat
  idx<-which(T1mat<alpha,arr.ind=TRUE)
  val<-T1mat[which(T1mat<alpha)]
  names(val)<-apply(idx,1,function(x) paste("(",x[2],",",x[1],")",sep="",collapse=""))
  if (length(val)>0)
     print(round(val,digits=3))
  else
     cat("none\n\n")
}

print.T1mobj<-function(x,alpha=0.05,...){
  txT1m<-"\nNonparametric RM model test: T1m (multidimensionality - reduced inter-item correlations)\n"
  writeLines(strwrap(txT1m, exdent=4))
  cat("    (counting cases with equal responses on both items)\n")
  cat("Number of sampled matrices:", x$n_eff,"\n")
  cat("Number of Item-Pairs tested:", length(x$prop),"\n")
  cat("Item-Pairs with one-sided p <", alpha,"\n")
  T1mmat<-x$T1mmat
  idx<-which(T1mmat<alpha,arr.ind=TRUE)
  val<-T1mmat[which(T1mmat<alpha)]
  names(val)<-apply(idx,1,function(x) paste("(",x[2],",",x[1],")",sep="",collapse=""))
  if (length(val)>0)
     print(round(val,digits=3))
  else
     cat("none\n\n")
}

print.T1lobj<-function(x,alpha=0.05,...){
  txt1<-"\nNonparametric RM model test: T1 (learning - based on item pairs)\n"
  writeLines(strwrap(txt1, exdent=4))
  cat("    (counting cases with reponsepattern (1,1) for item pair)\n")
  cat("Number of sampled matrices:", x$n_eff,"\n")
  cat("Number of sampled matrices:", x$n_eff,"\n")
  cat("Number of Item-Pairs tested:", length(x$prop),"\n")
  cat("Item-Pairs with one-sided p <", alpha,"\n")
  T1lmat<-x$T1lmat
  idx<-which(T1lmat<alpha,arr.ind=TRUE)
  val<-T1lmat[which(T1lmat<alpha)]
  names(val)<-apply(idx,1,function(x) paste("(",x[2],",",x[1],")",sep="",collapse=""))
  if (length(val)>0)
     print(round(val,digits=3))
  else
     cat("none\n\n")
}

print.T2obj<-function(x,...){
  prop<-x$prop
  idx<-x$idx
  stat<-x$stat
  statnam<-switch(stat,
     "var"="variance",
     "mad1"="mean absolute deviation",
     "mad2"="median absolute deviation",
     "range"="range"
  )
  txt<-"\nNonparametric RM model test: T2 (local dependence - model deviating subscales)\n"
  writeLines(strwrap(txt, exdent=4))
  cat("    (increased dispersion of subscale person rawscores)\n")
  cat("Number of sampled matrices:", x$n_eff,"\n")
  cat("Items in subscale:", idx,"\n")
  cat("Statistic:", statnam,"\n")
  cat("one-sided p-value:",prop,"\n\n")
#  cat("    (proportion of sampled",statnam," GE observed)\n\n")
}

print.T2mobj<-function(x,...){
  prop<-x$prop
  idx<-x$idx
  stat<-x$stat
  statnam<-switch(stat,
     "var"="variance",
     "mad1"="mean absolute deviation",
     "mad2"="median absolute deviation",
     "range"="range"
  )
  txt<-"\nNonparametric RM model test: T2m (multidimensionality - model deviating subscales)\n"
  writeLines(strwrap(txt, exdent=4))
  cat("    (decreased dispersion of subscale person rawscores)\n")
  cat("Number of sampled matrices:", x$n_eff,"\n")
  cat("Items in subscale:", idx,"\n")
  cat("Statistic:", statnam,"\n")
  cat("one-sided p-value:",prop,"\n\n")
#  cat("    (proportion of sampled",statnam," GE observed)\n\n")
}
print.T4obj<-function(x,...){
  prop<-x$prop
  idx<-x$idx
  gr.nam<-x$gr.nam
  gr.n<-x$gr.n
  alternative<-x$alternative
  cat("\nNonparametric RM model test: T4 (Group anomalies - DIF)\n")
  txt<-paste("    (counting", alternative, "raw scores on item(s) for specified group)\n", collapse="")
  writeLines(strwrap(txt, exdent=4))
  cat("Number of sampled matrices:", x$n_eff,"\n")
  cat("Items in Subscale:", idx,"\n")
  cat("Group:",gr.nam,"  n =",gr.n,"\n")
  cat("one-sided p-value:",prop,"\n\n")
#  cat("    (proportion of sampled raw scores GE observed)\n\n")
}

# removed in version 0.14-5
#print.T7obj<-function(x,...){
#  prop<-x$prop
#  cat("\nNonparametric RM model test: T7 (different discrimination - 2PL)\n")
#  txt<-"    (counting cases with response 1 on more difficult and 0 on easier item)\n"
#  writeLines(strwrap(txt, exdent=4))
#  cat("Number of sampled matrices:", x$n_eff,"\n")
#  cat("Item Scores:\n")
#  print(x$itscor)
#  cat("one-sided p-value:",prop,"\n\n")
#}
#print.T7aobj<-function(x,...){
#  prop<-x$prop
#  cat("\nNonparametric RM model test: T7a (different discrimination - 2PL)\n")
#  txt<-"    (counting cases with response 1 on more difficult and 0 on easier item)\n"
#  writeLines(strwrap(txt, exdent=4))
#  cat("Number of sampled matrices:", x$n_eff,"\n")
#  cat("Item Scores:\n")
#  print(x$itscor)
#  cat("\nItem-Pairs: (i>j ... i easier than j)\n\n")
#  print(round(prop,digits=3))
#}
print.T10obj<-function(x,...){
  spl.nam<-x$spl.nam
  prop<-x$prop
  hi.n<-x$hi.n
  low.n<-x$low.n
  txt<-"\nNonparametric RM model test: T10 (global test - subgroup-invariance)\n"
  writeLines(strwrap(txt, exdent=4))
  cat("Number of sampled matrices:", x$n_eff,"\n")
  cat("Split:",spl.nam,"\n")
  cat("Group 1: n = ",hi.n,"  Group 2: n =",low.n,"\n")
  cat("one-sided p-value:",prop,"\n\n")
#  cat("    (proportion of sampled statistics GE observed)\n\n")
}
print.T11obj<-function(x,...){
  prop<-x$prop
  txt<-"\nNonparametric RM model test: T11 (global test - local dependence)\n"
  writeLines(strwrap(txt, exdent=4))
  txt<-"    (sum of deviations between observed and expected inter-item correlations)\n"
  writeLines(strwrap(txt, exdent=4))
  cat("Number of sampled matrices:", x$n_eff,"\n")
  cat("one-sided p-value:",prop,"\n\n")
#  cat("    (proportion of sampled sums GE observed)\n\n")
}

