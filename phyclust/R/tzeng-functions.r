##-------------------------------------------------------------
## modified getI13.fun on Apr3, 2005; use it to replace I13.fun; in the new version, it can be
##                       used for both phased or unphased data
##-------------------------------------------------------------
## modified get.deriBigMatB.fun on Apr8, 2005; add dimensionality check
##   using chkdim.Multiply.fun
##         chkdim.Add.fun
##-------------------------------------------------------------


get.I.fun = function(y, Xt.mat, bigB.dPI, subPI.lst) {
  ## this function calculates I13 %*% I33i %*% t(I13), with input being:
  ## y:          data, length [nsim]
  ## Xt.mat:     data; dim [nsim, RR]
  ## bigB.dPI:   derivative info, calculated using construct.bigB.dPI,
  ##               and see generalform.tex for explaination;
  ## subPI.lst:  list of PIs with frequencies and names listed
  I13       = getI13.fun(y, Xt.mat, bigB.dPI)
  I13       = I13[-1,]  ## the first row is taken off, the base-line row. 
  I33i      = I33Inv.fun(subPI.lst, Xt.mat)
  I13t      = t(I13)
  Imat      = I13 %*% I33i %*% I13t
  ans       = list(I13 = I13, I33i = I33i, Imat = Imat)
  return(ans)
}


I33Inv.fun = function(subPI.lst, Xt.mat)
{
  ## this calculates I33-inverse matrix, and the names follow the same
  ## order as subPI.lst, and not Xt.mat. So Xt.mat's names have to agree
  ## with subPI.lst;
  
  ## subPI.lst: list of PIs with frequencies and names listed
  ## Xt.mat:    design matrix, of [nsim, RR]

  a = unlist(subPI.lst)
  a = a^2
  
  a.name = names(a)
  x.name = colnames(Xt.mat)
  b      = match(a.name, x.name)
  RR     = length(a.name)
  if(sum(b == (1:RR)) != RR)
    print("colnames(Xt.mat) and names(unlist(subPI.lst))) do not match - I33Inv.fun");
  Xt.mat = Xt.mat[,b]

  X = colSums(Xt.mat)
  ans = diag(a / (2 * X))
  rownames(ans) = a.name
  colnames(ans) = a.name

  return(ans)
}

       
#WCC getI13.fun = function(t1=y, Xt.mat, bigB.dPI) {
getI13.fun = function(t1, Xt.mat, bigB.dPI) {
  ## This function returns the I_{13} matrix
  ## Modified Apr3,2005; change y-ybar to t1, where t1= (y - mean(y))/a for phased data
  ##                                                t1= rep((y-mu)/a, nreps) * post for unphased data
  ## IMPORTANT!!!!!!!
  ## colnames(Xt.mat) should agree with rownames(bigB.dPI[[1]]), and
  ##  Xt.mat is sorted according to the rownames of bigB.dPI, so the
  ##  return matrix, I13's dimension names agree with bigB.dPI[[i]], and
  ##  not with colnames(Xt.mat)
  
  ## y is of length(n);
  ## Xt.mat is of dimension [n, RR]
  ## bigB.dPI is a list of length [RR], with each element being a
  ##   matrix of dimension [RR, RRstar]
  
  n      = length(t1)
  RR     = nrow(bigB.dPI[[1]])     ## original PI dimension;
  RRstar = ncol(bigB.dPI[[1]])     ## reduced PI dimension;
  I13    = matrix(0, RRstar, RR)   ## dimension of I13;
  o.name = rownames(bigB.dPI[[1]]) ## original PI names
  r.name = colnames(bigB.dPI[[1]]) ## reduced PI names
  
  ## check if o.name and Xt.mat names agree, and this sorts Xt.mat regardlessly
  x.name = colnames(Xt.mat)
  b      = match(o.name, x.name)
  if(sum(b == (1:RR)) != RR){    print("colnames(Xt.mat) and rownames(bigB.dPI[[1]]) do not match - getI13.fun")}
  Xt.mat = Xt.mat[,b]
  yX    = colSums( t1 * Xt.mat)  ## this should be the same as using ymat
  for(p in 1:RRstar)
    {
      for(k in 1:RR)
        {
          ## bigB.dPI[[k]][h,p] ==>
          ## k: pi that bigB.dPI is taking derivative w.r.t. to;
          ## h: pi that is in the original PI configuration;
          ## p: pi that is in the reduced PI configuration;
          aa = yX %*% bigB.dPI[[k]][,p];
          I13[p, k] = aa
        }
    }
  colnames(I13) = o.name;
  rownames(I13) = r.name;
  return(I13)
}


 ## I13.fun = function(y, Xt.mat, bigB.dPI) {
 ##   ## This function returns the I_{13} matrix
 ## 
 ##   ## IMPORTANT!!!!!!!
 ##   ## colnames(Xt.mat) should agree with rownames(bigB.dPI[[1]]), and
 ##   ##  Xt.mat is sorted according to the rownames of bigB.dPI, so the
 ##   ##  return matrix, I13's dimension names agree with bigB.dPI[[i]], and
 ##   ##  not with colnames(Xt.mat)
 ## 
 ##   ## y is of length(n);
 ##   ## Xt.mat is of dimension [n, RR]
 ##   ## bigB.dPI is a list of length [RR], with each element being a
 ##   ##   matrix of dimension [RR, RRstar]
 ##   
 ##   n      = length(y)
 ##   RR     = nrow(bigB.dPI[[1]])     ## original PI dimension;
 ##   RRstar = ncol(bigB.dPI[[1]])     ## reduced PI dimension;
 ##   I13    = matrix(0, RRstar, RR)   ## dimension of I13;
 ##   ybar   = mean(y)
 ##   o.name = rownames(bigB.dPI[[1]]) ## original PI names
 ##   r.name = colnames(bigB.dPI[[1]]) ## reduced PI names
 ## 
 ##   ## check if o.name and Xt.mat names agree, and this sorts Xt.mat
 ##   ## regardlessly
 ##   x.name = colnames(Xt.mat)
 ##   b      = match(o.name, x.name)
 ##   if(sum(b == (1:RR)) != RR)
 ##     print("colnames(Xt.mat) and rownames(bigB.dPI[[1]]) do not match - I13.fun");
 ##   Xt.mat = Xt.mat[,b]
 ##   
 ##   ## ymat = matrix( rep( y - ybar, RR), ncol = RR)
 ##   ## yX   = colSums(ymat * Xt.mat)
 ##   yX    = colSums( (y - ybar) * Xt.mat)  ## this should be the same as using ymat
 ## 
 ##   for(p in 1:RRstar)
 ##     {
 ##       for(k in 1:RR)
 ##         {
 ##           ## bigB.dPI[[k]][h,p] ==>
 ##           ## k: pi that bigB.dPI is taking derivative w.r.t. to;
 ##           ## h: pi that is in the original PI configuration;
 ##           ## p: pi that is in the reduced PI configuration;
 ##           a = yX %*% bigB.dPI[[k]][,p];
 ##           I13[p, k] = a
 ##         }
 ##     }
 ##   colnames(I13) = o.name;
 ##   rownames(I13) = r.name;
 ##   return(I13)
 ## }



zero.mat.fun = function(mat) {
  ## making the input matrix zero in every cell;
  nr = nrow(mat);
  nc = ncol(mat);
  rn = rownames(mat);
  cn = colnames(mat);
  mat = matrix(0, nr, nc);
  rownames(mat) = rn;
  colnames(mat) = cn;
  return(mat)
}

index.k.fun<-function(len.subPI = 6) {
  ## this function creats a bit strange index that keeps track of
  ## B matrices in the BB-list, with each B indicates how it transforms
  ## layers of haplotypes to their appropriate ancestors;

  ##             __[[1]]  
  ##            |   
  ##  pi(0)     1   __[[2]]
  ##  pi(1)     2  |   __[[3]]
  ##  pi(2)     4  3  |   __[[4]]
  ##  pi(3)     7  6  5  |   __[[5]] 
  ##  pi(4)    11 10  9  8  |   __[[6]]
  ##  pi(5)    16 15 14 13 12  |   __[[7]]
  ##  pi(6)    22 21 20 19 18 17  |   __[[8]]
  ##  pi(7)    29 28 27 26 25 24 23  |
  ##  pi(8)    37 36 35 34 33 32 31 30

  ## the numbers indicate which member of the BB-list belongs to which layer of pi;
  ## [[num]] is the return list in which the vertical numbers are sorted;
  ## NOTE: because pi(num) starts at 0, so len.subPI = 6 ==> pi(5)!!!!
  ## NOTE: the last [[num]] is always set to be NA!!!
  
  sseq <- list()
  KK <- len.subPI - 1
  if(KK == 0)
    {
      sseq[[1]] <- 1
    } else {
      sseq[[1]] <- 1 + cumsum(1:KK)
      if(KK > 1)
        {
          for(k in 2:KK) {
            sseq[[k]] <- sseq[[k-1]][-1]-1
          }
        }
      sseq[[1]] = c(1, sseq[[1]])
      sseq[[(KK+1)]] = NA
    }
  return(sseq)
}


##--------
      chkdim.Multiply.fun<-function(mat){
        chk1<-length(mat)>0
        chk2<-all(dim(mat)>0)
        chk3<-sum(is.na(mat))==0
        return(list("chk1"=chk1,"chk2"=chk2, "chk3"=chk3))
      }

#WCC      chkdim.Add.fun<-function(mat, ncolShouldBe=RRstar){
      chkdim.Add.fun<-function(mat, ncolShouldBe){
        if(length(mat)==0)
          {
            chk1=chk2=0
          }else{
            chk1<-all(dim(mat)>0)
            chk2<-(ncol(mat) == ncolShouldBe)
          }
        return(chk1*chk2)
      }


#WCC get.deriBigMatB.fun<-function(BB.deri, subPI.lst, subPi.num=k)
get.deriBigMatB.fun<-function(BB.deri, subPI.lst, subPi.num)
  { 
    ## modified on Apr 8 2005, add the compatability check of dimensionalith when getting
    ## preBigBB.deri and BigMatB, see "*******".
    ##
    ## this calculate one bigB.dPI (w.r.t. to one pi). For reference,
    ## see generalform.tex;

    ## NOTE: this function is not easy to use, because it does not say,
    ##  specifically which pi that the function is calculating w.r.t.;
    ##  Instead, one has to know how to construct the correct BB.deri,
    ##  which should have the appropriate cells filled out (see
    ##  construct.bigB.dPI) by derivative information. 

    ## BB.deri: It's a BB-list, but for this function, it should have
    ##         appropriate cells filled out with derivative calculated
    ##         information.
    ## subPI.lst:  list of PIs with frequencies and names listed
    ## subPi.num:  the number of "sub-group of PIs" that the derivative
    ##              pi belongs. The number should not exceed length(subPI.lst); 
    
    mmat         <-get.allpossROUTE.fun()
    r.ele        <-c( 1,1,2,4,8,16,32,64,128,256)
    cumsum.r     <-cumsum(r.ele)
    lenPI        <-length(subPI.lst)  ## 4
    lenBB        <-length(BB.deri)     ## 7=1+1+2+3
    lenROUTE     <-cumsum.r[lenPI]    ## 8=1+1+2+4
    RRstar       <-length(subPI.lst[[1]])
    ### changed on 03/11/05 (remove "if" condition)
    mat.route<-as.matrix(mmat[(1:lenROUTE),(1:lenBB)]);
    
    mat.route[,1]<-0 ### this adjustment (made on 03/12/05 is for calculating
                     ### deri.B only, because BB.deri[[1]] is alwyas a 0 matrix.
                     ### Hence has to be exclude from the product calculation,
                     ### o.w. making the final product=0
    pi.ele  <-rep((1:length(r.ele)), r.ele)
    j.index <-pi.ele[1:cumsum.r[lenPI]] ###  indexing a mat B is for the level of H^{(j)}

    preBigBB.deri<-list()
    preBigBB.deri[[1]]<-BB.deri[[1]]
    if(lenROUTE>1)
      {
        for(kk in 2:lenROUTE)
          {
            bb<-BB.deri[mat.route[kk,]]
            BB<-bb[[length(bb)]]
            if(length(bb)>1)
              {
                ##***** check the compability of dimenationality  before multiplying
                chk<-sapply(bb, chkdim.Multiply.fun)
                if(all(chk==T))
                  {
                    for(jj in max(1,(length(bb)-1)):1){BB<-BB%*%bb[[jj]]}### range of jj is modified on 03/12/05 as well
                  }
                else
                  {
                    BB<-matrix(0, nrow=length(subPI.lst[[ j.index[[kk]]  ]]), ncol=RRstar)
                    colnames(BB) = names(subPI.lst[[1]])
                    rownames(BB) = names(subPI.lst[[ j.index[[kk]]  ]])
                  }
            } 
            preBigBB.deri[[kk]]<-BB
          }
      }

    ## If subPi.num is not in the first (reduced group),
    ##  then some of the preBigBB.deri matrices have to be zeroed.     
    index = index.k.fun(lenPI)
    if(subPi.num > 1)
      {
        for(kk in 1:(subPi.num-1))
          {
            index.k = index[[kk]]
            for(jj in index.k)
              {
                preBigBB.deri[[jj]] = zero.mat.fun(preBigBB.deri[[jj]])
              }
          }
      }
    
    indx<-split((1:length(preBigBB.deri)), j.index)
    BigMatB<-NULL
    for(jj in 1:length(indx))
      {
        tmpB<-preBigBB.deri[ indx[[jj]] ]
        if(length(tmpB)==1)
          {
            BigMatB<-rbind(BigMatB, tmpB[[1]])
          }else {
            ##***** check the compability of dimenationality  before adding
            key<-sapply(tmpB,chkdim.Add.fun, ncolShouldBe=RRstar); pos<-(1:length(tmpB))[key==1]
            if(all(pos==0))
              {
                sumB = NULL
              } else {
                Add.tmpB<-tmpB[pos]
                sumB    <-Add.tmpB[[1]]
                for(ll in 2:length(Add.tmpB)){sumB<-sumB+Add.tmpB[[ll]]}
              }
            BigMatB<-rbind(BigMatB, sumB)
          }
      }
    return(BigMatB)
  }
            
## old version 
## get.deriBigMatB.fun<-function(BB.lst, subPI.lst, subPi.num=k)
##   { 
##     ## this calculate one bigB.dPI (w.r.t. to one pi). For reference,
##     ## see generalform.tex;
## 
##     ## NOTE: this function is not easy to use, because it does not say,
##     ##  specifically which pi that the function is calculating w.r.t.;
##     ##  Instead, one has to know how to construct the correct BB.lst,
##     ##  which should have the appropriate cells filled out (see
##     ##  construct.bigB.dPI) by derivative information. 
## 
##     ## BB.lst: It's a BB-list, but for this function, it should have
##     ##         appropriate cells filled out with derivative calculated
##     ##         information.
##     ## subPI.lst:  list of PIs with frequencies and names listed
##     ## subPi.num:  the number of "sub-group of PIs" that the derivative
##     ##              pi belongs. The number should not exceed length(subPI.lst); 
##     
##     mmat         <-get.allpossROUTE.fun()
##     r.ele        <-c( 1,1,2,4,8,16,32,64,128,256)
##     cumsum.r     <-cumsum(r.ele)
##     lenPI        <-length(subPI.lst)  ## 4
##     lenBB        <-length(BB.lst)     ## 7=1+1+2+3
##     lenROUTE     <-cumsum.r[lenPI]    ## 8=1+1+2+4
## 
##     ### changed on 03/11/05 (remove "if" condition)
##     mat.route<-as.matrix(mmat[(1:lenROUTE),(1:lenBB)]);
## 
##     
##     mat.route[,1]<-0 ### this adjustment (made on 03/12/05 is for calculating
##                      ### deri.B only, because BB.lst[[1]] is alwyas a 0 matrix.
##                      ### Hence has to be exclude from the product calculation,
##                      ### o.w. making the final product=0
##     pi.ele<-rep((1:10), r.ele)
##     j.index <-pi.ele[1:cumsum.r[lenPI]] ###  indexing a mat B is for the level of H^{(j)}
## 
##     preBigBB.lst<-list()
##     preBigBB.lst[[1]]<-BB.lst[[1]]
##     if(lenROUTE>1)
##       {
##         for(kk in 2:lenROUTE)
##           {
##             bb<-BB.lst[mat.route[kk,]]
##             BB<-bb[[length(bb)]]
##             if(length(bb)>1){ for(jj in max(1,(length(bb)-1)):1){BB<-BB%*%bb[[jj]]} } ### range of jj is modified on 03/12/05 as well
##             preBigBB.lst[[kk]]<-BB
##           }
##       }
## 
##     ## If subPi.num is not in the first (reduced group),
##     ##  then some of the preBigBB.lst matrices have to be zeroed.     
##     index = index.k.fun(lenPI)
##     if(subPi.num > 1)
##       {
##        for(kk in 1:(subPi.num-1))
##          {
##            index.k = index[[kk]]
##            for(jj in index.k)
##              {
##                preBigBB.lst[[jj]] = zero.mat.fun(preBigBB.lst[[jj]])
##              }
##          }
##       }
## 
## 
##     indx<-split((1:length(preBigBB.lst)), j.index)
##     BigMatB<-NULL
##     for(jj in 1:length(indx))
##       {
##         tmpB<-preBigBB.lst[ indx[[jj]] ]
##         if(length(tmpB)==1)
##           {
##             BigMatB<-rbind(BigMatB, tmpB[[1]])
##           }else{
##             sumB<-0
##             for(ll in 1:length(tmpB)){sumB<-sumB+tmpB[[ll]]}
##             BigMatB<-rbind(BigMatB, sumB)
##           }
##       }
##     return(BigMatB)
##   }


  construct.dB.list.fun = function(CC.mat, subPI.lst) {
    ## modified on Apr 10, 2005, see **** parts to take care of "NA" and "numeric(0)" in CC.mat
    ##
    ## this constructs the dB.list, see reference in generalform.tex for
    ## every CC.mat, calculate, w.r.t. every PI in that CC.mat,
    ## partial(B)/partial(pi). See pt1 -- pt6 calculations in
    ## generalform.tex. The answer it returns is a list of RR (length of
    ## subPI.lst).
    
    ## Note: Because the input CC.mat is of different dimension, the
    ##        output is a list of RR matrices, each with the same dim
    ##        as the input CC.mat matrix.
    ## Note: The derivative is taken with respect to ALL pi's, and if
    ##        any pi is not part of the CC.mat, then the answers are
    ##        automatically zero.
    ## Note: If the CC.mat is all zero, the outcome is automatically
    ##        all ZERO. 
    
    ## CC.mat is ONE CC matrix;
    ## subPI.lst: list of PIs with frequencies and names listed

    ##********* begin modificiation on 04.10.2005 ****************
    if(sum(is.na(CC.mat)) > 0  | length(CC.mat)==0 )
      {
        ## when CC.mat is NA or numeric(0), make dB.list=NULL
          full.PI  = unlist(lapply(subPI.lst, names))
          l        = length(full.PI)
          dB.list  = list()
          for(i in 1:l)
            {
              dB.list[[i]] = numeric(0)
              attr(dB.list[[i]], "name") = full.PI[i]
            }
          ##********* begin modificiation on 04.10.2005 ****************
      } else if(all(CC.mat == 0))
        {
          full.PI  = unlist(lapply(subPI.lst, names))
          l        = length(full.PI)
          dB.list  = list()
          for(i in 1:l)
            {
              dB.list[[i]] = CC.mat * 0
              attr(dB.list[[i]], "name") = full.PI[i]
            }
        } else {
          target.PI = colnames(CC.mat)  #pi^C in the doc
          full.PI = unlist(lapply(subPI.lst, names))
          
          for(i in 1:length(subPI.lst))
            {
              if(!is.na(sum(match(target.PI, names(subPI.lst[[i]])))))
                {
                  k = i;
                  break;
                }
            }
          subPI = subPI.lst[[k]]         ##freq of target.PI
          
          ##pt1 = pt2 = pt3 = pt4 = pt5 = pt6 = CC.mat;
          
          nr  = nrow(CC.mat)
          nc  = ncol(CC.mat)
          
          pt5 = matrix(rep(subPI, nr), nrow = nr, byrow=T)
          
          deno = rowSums(CC.mat * pt5)
          a    = 1/deno;       a[deno==0] = 0
          pt1 = matrix(rep(a, nc), ncol = nc)
          
          pt3 = pt1^2
          pt4 = CC.mat
          
          l = length(full.PI)
          dB.list = list()
          
          ## dB.list is a list of length l, with each element being a
          ## derivative taken with respect to the name attributes. Therefore,
          ## if the name is not an element of target.PI, that dB.list is left
          ## to be zero. 
          
          for(i in 1:l)
            {
              dB.list[[i]] = matrix(0, nr, nc)
              dB.list[[i]] = CC.mat * 0
              attr(dB.list[[i]], "name") = full.PI[i]
            }
          
          no.zero.list = match(full.PI, target.PI)
          
          for(j in 1:l)
            {
              if(!is.na(no.zero.list[j]))
                {
                  i = no.zero.list[j]
                  pt2 = matrix(0, nr, nc)
                  pt2[,i] = CC.mat[,i]
                  
                  pt6 = matrix(rep(pt2[,i], nc), ncol = nc)
                  dB.list[[j]] = pt1 * pt2 - pt3 * pt4 * pt5 * pt6
                  attr(dB.list[[j]], "name") = full.PI[j]
                }
            }
        }
    
    return(dB.list)
  }


construct.bigB.dPI = function(dB.all, subPI.lst, BB.lst) 
{
  ## this function constructs the (bigB.dPI) expression in
  ## generalform.tex. In reality, it calculates derivative of B
  ## w.r.t. to ALL PI's. So for RR-pi, the function returns a list of
  ## length(RR), each is a matrix of the same dimension as BB.

  ## Note: this function calls get.deriBigMatB.fun(), the usage of
  ##        which is not very intuitive. Basically it refills a structure
  ##        that is the same as BB.lst, and use a similar way to construct
  ##        the derivative matrices for every pi. This function can be
  ##        better written later. 
  
  #dB.all: is returned as functions of CC, such as 
  ##    dB.all = list()
  ##    for(i in 1:length(CC)) {
  ##      dB.all[[i]] = construct.dB.list.fun(CC[[i]], subPI.lst)
  ##    }

  #subPI.lst: list of PIs with frequencies and names listed, the
  ##            names have to agree with CC!!!!! 
  #BB.lst: BB.list 

  full.PI = unlist(lapply(subPI.lst, names))
  reduced.PI = names(subPI.lst[[1]])
  l = length(full.PI)
  m = length(reduced.PI)
  bigB.dPI = list()
  for(i in 1:l) {
    bigB.dPI[[i]] = matrix(0, l, m)
    rownames(bigB.dPI[[i]]) = full.PI
    colnames(bigB.dPI[[i]]) = reduced.PI
    attr(bigB.dPI[[i]], "name") = full.PI[i]
  }

  len.subPI = length(subPI.lst)
  index = index.k.fun(len.subPI)  

  for(i in 1:l)
    {
      ##let's do full.PI[i] one at a time;
      current.PI = full.PI[i];
      ##find out which subPI that current.PI belongs to;
      for(j in 1:len.subPI)
        {
          if(!is.na(match(current.PI, names(subPI.lst[[j]]))))
            {
              k = j;
              break;
            }
        }
      ## current.PI belongs to subPI.lst[[k]];
      
      ## make a copy of BB, and change components according to build
      ##   bigB.dPI[[i]]
      BB.deri = BB.lst
      
      ## index associated with k: 
      index.k = index[[k]]
      if(!is.na(sum(index.k)))
        {
          for (jj in index.k)
            {
              BB.deri[[jj]] = dB.all[[jj]][[i]]
            }
          ## Since BB.deri has the same structure as the BB.lst, we are
          ##  using a modified old function to find the derivative matrix. 
          bigB.dPI[[i]] <-get.deriBigMatB.fun(BB.deri, subPI.lst, k)
          attr(bigB.dPI[[i]], "name") = full.PI[i]
        }
    }
  return(bigB.dPI)
}




get.allpossROUTE.fun<-function()
{
  ## when this file is modified, function "getPIstar.fun" need to be changed as well!!!!!!
  m1 <-list(1)
  m2 <-  lapply(m1, 2, FUN = c)
  m3 <-c(lapply(m2, 3, FUN = c),
         lapply(m1, 4, FUN = c))
  m4 <-c(lapply(m3, 5, FUN = c),
         lapply(m2, 6, FUN = c),
         lapply(m1, 7, FUN = c))
  m5 <-c(lapply(m4, 8, FUN = c),
         lapply(m3, 9, FUN = c),
         lapply(m2, 10, FUN = c),
         lapply(m1, 11, FUN = c))
  m6 <-c(lapply(m5, 12, FUN = c),
         lapply(m4, 13, FUN = c),
         lapply(m3, 14, FUN = c),
         lapply(m2, 15, FUN = c),
         lapply(m1, 16, FUN = c))
  m7 <-c(lapply(m6, 17, FUN = c),
         lapply(m5, 18, FUN = c),
         lapply(m4, 19, FUN = c),
         lapply(m3, 20, FUN = c),
         lapply(m2, 21, FUN = c),
         lapply(m1, 22, FUN = c))
  m8 <-c(lapply(m7, 23, FUN = c),
         lapply(m6, 24, FUN = c),
         lapply(m5, 25, FUN = c),
         lapply(m4, 26, FUN = c),
         lapply(m3, 27, FUN = c),
         lapply(m2, 28, FUN = c),
         lapply(m1, 29, FUN = c))
  m9 <-c(lapply(m8, 30, FUN = c),
         lapply(m7, 31, FUN = c),
         lapply(m6, 32, FUN = c),
         lapply(m5, 33, FUN = c),
         lapply(m4, 34, FUN = c),
         lapply(m3, 35, FUN = c),
         lapply(m2, 36, FUN = c),
         lapply(m1, 37, FUN = c))
  m10<-c(lapply(m9, 38, FUN = c),
         lapply(m8, 39, FUN = c),
         lapply(m7, 30, FUN = c),
         lapply(m6, 41, FUN = c),
         lapply(m5, 42, FUN = c),
         lapply(m4, 43, FUN = c),
         lapply(m3, 44, FUN = c),
         lapply(m2, 45, FUN = c),
         lapply(m1, 46, FUN = c))
  
  m = c(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10)
  d = 10
  r = c(1, rep(2, d-1))^(c(0,0:(d-2)))
  di = d * (d-1) / 2 + 1
  
  l = rep(1:d, r)
  mat = matrix(0, sum(r), di)
  for(i in 1:sum(r)) {
    mat[i,m[[i]]] = m[[i]]
  }
  row.names(mat)<-paste("Step ",
                        c(1,2,rep(3,2), rep(4,4), rep(5,8), rep(6,16), rep(7,32), rep(8,64), rep(9,128), rep(10, 256)),
                        ".",
                        c(1,1,   (1:2),    (1:4),    (1:8),    (1:16),    (1:32),    (1:64),    (1:128), (1:256)),   sep="")
  dimnames(mat)[[2]]<-paste("mat",c(1,2,rep(3,2), rep(4,3), rep(5,4), rep(6,5),rep(7,6)),".(",(1:di),")",sep="")
  return(mat)
}


    
