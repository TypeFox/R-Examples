

#getcut.fun<-function(pp.org,nn=2*ngeno, plot=0){
#      pp     <-rev(sort(pp.org))
#      ct     <-round(pp*nn)
#      dimen  <-log((1:length(pp)), base=2 )/ct
#      info   <-cumsum(pp*log(1/pp))
#      netinfo<-info-dimen
#      cutpos <-netinfo==max(netinfo)
#      if(plot==1){
#        plot( netinfo, pp,type="b",cex=0.5)
#        barplot(pp.org)
#        abline(h=pp[cutpos])
#      }
#      return( (1:length(pp))[cutpos])
#    }


    final.BigMatB.fun<-function(mmle.cscn, hapreserv, loci, digit){
      ini.cscn     <-get.ini.possHap.hapReserv.fun(mmle.cscn, loci, digit, hapreserv)
      mut.mat.cscn <-ini.cscn$prepsi
      subPI.cscn   <-get.subPI.fun(mmle.cscn,  mut.mat=mut.mat.cscn)
      ## get BB matrix
      jnk.cscn     <-get.preBB.skip.fun(subPI.cscn, loci, digit, mmle.cscn)
      preBB.cscn   <-jnk.cscn$preBB.lst
      indexBB.cscn <-jnk.cscn$indexBB
      medBB.cscn   <-get.medBB.skip.fun(subPI.lst=subPI.cscn, preBB.lst=preBB.cscn, indexBB.lst=indexBB.cscn,loci, digit, mmle.cscn)
      fullBB.cscn  <-get.fullBB.fun(medBB.lst=medBB.cscn,  subPI.lst=subPI.cscn)
      BB.cscn      <-get.BB.fun(fullBB.cscn)
      BigMatB.cscn <-getBigMatB.fun(BB.cscn, subPI.cscn)
      return(list(fullBB.cscn=fullBB.cscn,
                  ##BB.cscn=BB.cscn,
                  BigMatB.cscn =BigMatB.cscn ,
                  subPI.cscn =  subPI.cscn
                  ))
  }



##--------------------------------------------------------------------------------------
## modified in "get.fullBB.fun":    if(aa+bb==2 && is.null(medBB.lst[[rr]][[mm]])==0)
## 
## on 03/10/05: modified in
##   get.preBB.skip.fun: use mmle.cscn to construct matB
##   get.medBB.skip.fun: use mmle.cscn to construct matB
##   getPIstar.fun:      remove the "if" statement
## on 03/11/05:  add a new function:
##   getBigMatB.fun: BB.lst=BB.cs, subPI.lst=subPI.cs
##
## on 04/06/05: add dimensionality check in function
##   getBigMatB.fun, by using chkdim.Multiply.fun, and chkdim.Add.fun
##   modified part is annotated by ******
##
##--------------------------------------------------------------------------------------


##----
get.step.mut.fun<-function(loci, digit,  Ncol, Nrow,NAMEcol, NAMErow){
  ## Ncol  & NAMEcol = col variable
  ## Nrow  & NAMErow = row variable
  UU<-matrix(0, ncol=loci, nrow=Ncol)
  for(i in 1:loci){    UU[,i]<-as.numeric(substring(NAMEcol,   (i-1)*digit+1,i*digit))}
  VV<-matrix(0, ncol=loci, nrow=Nrow)
  for(i in 1:loci){  VV[,i]<-as.numeric(substring(NAMErow,(i-1)*digit+1,i*digit))}
  step.mut<-matrix(0, nrow=Nrow, ncol=Ncol)
  for(v in 1:Nrow){
    tmp<-rep(0,Ncol); for(i in 1:loci){tmp<-tmp+abs(VV[v,i]-UU[,i])}
    step.mut[v,]<-tmp}
  row.names(step.mut)<-NAMErow
  dimnames(step.mut)[[2]]<-NAMEcol
  return(step.mut)}



get.ini.possHap.hapReserv.fun<-function(em, loci, digit, prehapReserv){
  possHap  <-names(em)
  RR       <-length(possHap)
  tmp      <-prehapReserv[  match(possHap, prehapReserv)]
  hapReserv<-tmp[is.na(tmp)==0]
  Rstar     <-length(hapReserv)
  ## cat("RR=", RR, " Rstar=", Rstar,"\n")
  prepsi    <-get.step.mut.fun(loci, digit,  RR, Rstar, possHap, hapReserv )
  return(list(possHap=possHap,
              hapReserv=hapReserv,
              RR=RR,
              Rstar=Rstar,
              prepsi=prepsi))
}

#################################################
## get matrix BB
#################################################
get.subPI.fun<-function(em, mut.mat){
  subPI.lst<-list()
  allhap <-names(em)
  if(nrow(mut.mat)==ncol(mut.mat)){    subPI.lst[[1]]<-em  }else{
    mmin   <-apply(mut.mat, 2, min)
    tab    <-table(mmin)
    subPI.lst<-list()
    subPI.lst[[1]]<-tmp<-em[mmin==0]
    for(rr in 2:(1+max(mmin))){
      del  <-match(names(tmp), names(em))
      subPI.lst[[rr]]<-tmp1<- (em[-del])[apply(matrix(mut.mat[,-del], ncol=(length(em)-length(del)))==(rr-1), 2, sum)>0]
      tmp<-c(tmp1, tmp)
    }
  }
  return(subPI.lst)
}
#WCC get.preBB.skip.fun<-function(subPI.lst=subPI.cs, loci, digit, mmle.cscn){
get.preBB.skip.fun<-function(subPI.lst, loci, digit, mmle.cscn){
  preBB.lst<-list(); indexBB.lst<-list()
  len      <-length(subPI.lst)
  preBB.lst[[1]]  <-diag(rep(1, length(subPI.lst[[1]])));  dimnames(preBB.lst[[1]])<-list(names(subPI.lst[[1]]),names(subPI.lst[[1]]))
  indexBB.lst[[1]]<-0
  if(len>1)
    {
      sum.subPI.lst<-sapply(subPI.lst, sum)  
      for(rr in 2:len){
        if(sum.subPI.lst[rr]==0){## no "rr"-step descendants
          preBB.lst[[rr]]<-NA
        } else
        if(sum.subPI.lst[rr]> 0)
          {
            index.all<-0; ss<-0
            while(index.all==0)
              {## in order to find the right "ss"-step ancestors to obtain mut.step
                ss<-ss+1
                if(sum.subPI.lst[rr-ss]>0)
                  { 
                    mut.step<-as.data.frame(get.step.mut.fun(loci, digit, length(subPI.lst[[rr-ss]]), length(subPI.lst[[rr]]),
                                                             names(subPI.lst[[rr-ss]]), names(subPI.lst[[rr]])))
                    mmin      <-apply(mut.step, 1, min)
                    index.all<-max(mmin==ss)
                  } else {index.all<-0}
              }
            names(mmin)<-names(subPI.lst[[rr]])
            indexBB.lst[[rr]]<-mmin
            theta<-0
            ## changed on 03/10/2005, should use mmle.cscn instead of pcs or pcn to construct matB
            poolPI<-mmle.cscn[match(names(subPI.lst[[rr-ss]]), names(mmle.cscn))]
            if(length(poolPI)>1){
                tmp <-(theta^(mut.step-ss))%*% diag(poolPI) }else{
                tmp <-(theta^(mut.step-ss)) *       poolPI  }
            preBB.lst[[rr]]<- tmp/    apply(tmp, 1,sum)
            dimnames(preBB.lst[[rr]])<-dimnames(mut.step)
          }
      }## end of rr-loop
    }## end of if(len>1)
  return(list(preBB.lst=preBB.lst,                indexBB.lst=indexBB.lst))
}


##  subPI.lst<-subPI.cs[[j]][[ii]];preBB.lst<-preBB.cs[[j]][[ii]];indexBB.lst<-indexBB.cs[[j]][[ii]]
##  rm(medBB.lst, subPI.lst, preBB.lst, indexBB.lst, len, rr, ss,tt,tmp, tmp1, tmp0, rowsum.tmp,pos,sstop, chk, prechk)

#WCC get.medBB.skip.fun<-function(subPI.lst=subPI.cn, preBB.lst=preBB.cn, indexBB.lst=indexBB.cn,loci, digit, mmle.cscn){
get.medBB.skip.fun<-function(subPI.lst, preBB.lst, indexBB.lst, loci, digit, mmle.cscn){
  ## if some hap couldn't find a direct 1-step ancestor, then need to use medBB to construct matB
  medBB.lst<-list()
  medBB.lst[[1]]<-list()
  len      <-length(subPI.lst)
  medBB.lst[[1]][[1]]<-preBB.lst[[1]]
  if(len>1)
    {
      for(rr in 2:len)
        {
          medBB.lst[[rr]]<-list()
          if(prod(indexBB.lst[[rr]]==1))## i.e., each hap has its 1-step mutant ancestor
            {
              medBB.lst[[rr]][[1]]<-preBB.lst[[rr]]
            } else
          if(sum(preBB.lst[[rr]]=='NaN')==0)## i.e., from preBB we've found its closet mutant ancestor
            {
              medBB.lst[[rr]][[1]]<-preBB.lst[[rr]]
            }else{## i.e., some hap dont have 1-step ancestors, need to prepare "BB-1step", "BB-2step", ...etc
              ##--- BB-1step ---
              tmp0<-preBB.lst[[rr]];tmp0[tmp0=='NaN']<-0
              medBB.lst[[rr]][[1]]<-tmp0
              chk<-1; ss<-(min(indexBB.lst[[rr]])+1)
              ##for(ss in (min(indexBB.lst[[rr]])+1):(rr-1)){ ##--- BB-"ss"step ---
              while(chk>0 & ss<rr)
                {
                  pos<-(indexBB.lst[[rr]]>=ss )
                  if(sum(pos)==0) ## i.e., all hap have no "ss" step ancetors
                    {
                      medBB.lst[[rr]][[ss]]<-NULL
                    }else{## i.e., some hap have "ss" step ancestors
                      ## need to get the BBmatrix for rr to rr-ss
                      ## note rr=# of mutations needed for current group of hap to become the reserved hap
                      sstop<-0; tt<-ss
                      while(sstop==0)
                        {## in order to find the right "tt"-step ancestors to obtain mut.step
                          if(sum(subPI.lst[[rr-tt]])==0)
                            {
                              tt<-tt+1
                            }else{
                              mut.step<-as.data.frame(get.step.mut.fun(loci, digit, Nrow=length(pos),NAMErow=names(subPI.lst[[rr]]),
                                                                       Ncol=length(subPI.lst[[rr-tt]]),NAMEcol = names(subPI.lst[[rr-tt]])))
                              mmin     <-apply(mut.step, 1, min)
                              sstop   <-prod(mmin[pos]==tt); if(sstop==0){tt<-tt+1}
                              ## here I assume all hap will find their ancestor at the same tt-step, this need to be modified later
                            }
                        }## end of while
                      theta<-0
                      ## changed on 03/10/2005, should use mmle.cscn instead of pcs or pcn to construct matB
                      poolPI<-mmle.cscn[match(names(subPI.lst[[rr-tt]]), names(mmle.cscn))]
                      if(length(poolPI)>1)
                        {
                          tmp<-(theta^(mut.step-tt))%*% diag(poolPI)
                        }else{
                          tmp<-(theta^(mut.step-tt))  *      poolPI
                        }
                      rowsum.tmp<-apply(tmp, 1,sum); rowsum.tmp[rowsum.tmp==0]<-1
                      tmp1<- tmp/rowsum.tmp
                      medBB.lst[[rr]][[ss]]<-tmp1*0
                      medBB.lst[[rr]][[ss]][pos,]<-tmp1[pos,]
                      dimnames(medBB.lst[[rr]][[ss]])<-dimnames(mut.step)
                    }
                  ## check if each hap finds its ancestor
                  if(length(medBB.lst[[rr]])>1){
                    tmp<-lapply(medBB.lst[[rr]],function(mat){if(1-is.null(mat)){apply(mat, 1, sum)}else{NA}})
                    prechk<-sapply(tmp, sum)
                  }else{  prechk<-apply(medBB.lst[[rr]][[1]],1, sum)}
                  chk<-sum(prechk==0, na.rm=T)   ## chk=0 means all category can find its closest ancestor; chk>0 means some havent found yet.
                  if(chk>0){ss<-ss+1}
                }## end of while
            }
        }## end of rr-loop
    }
  return(medBB.lst)
}




#WCC get.fullBB.fun<-function(medBB.lst=medBB.cs,subPI.lst=subPI.cs)
get.fullBB.fun<-function(medBB.lst,subPI.lst)
{
  name.lst <-lapply(subPI.lst, names)
  dim.vec  <-sapply(subPI.lst, length)
  len      <-length(dim.vec)

  fullBB.lst<-list()
  fullBB.lst[[1]]<-medBB.lst[[1]]
  if(len>1)
    {
      fullBB.lst[[2]]<-medBB.lst[[2]]
      if(len>2)
        {
          for(rr in 3:length(medBB.lst))
            {
              fullBB.lst[[rr]]<-list()
              for(kk in (1:(rr-1)))
                {
                  mat<-matrix(0, nrow=dim.vec[rr], ncol=dim.vec[rr-kk])
                  row.names(mat)    <-name.lst[[rr]]
                  dimnames(mat)[[2]]<-name.lst[[rr-kk]]
                  fullBB.lst[[rr]][[kk]]<-mat
                }
              
              for(kk in (1:(rr-1)))
                {
                  for(mm in 1:length(medBB.lst[[rr]]))
                    {
                      aa<-sum(is.na(match(dimnames(fullBB.lst[[rr]][[kk]])[[2]], dimnames( medBB.lst[[rr]][[mm]])[[2]] )))==0
                      bb<-sum(is.na(match(dimnames( medBB.lst[[rr]][[mm]])[[2]], dimnames(fullBB.lst[[rr]][[kk]])[[2]] )))==0
                      if(aa+bb==2 && is.null(medBB.lst[[rr]][[mm]])==0)
	              ##   if(aa+bb==2)
                        {
                          fullBB.lst[[rr]][[kk]]<-medBB.lst[[rr]][[mm]]
                        }
                    }## end of mm-loop
                }## end of kk-loop

              
            }## end of rr-loop
        }## end of if(len>2)
    }## end of if (len>1)
  return(fullBB.lst)

}


get.BB.fun<-function(fullBB.lst){
  BB.lst  <-list()
  len     <-c(1,(1:100))
  cumlen  <-cumsum(len)
  BB.lst[[1]]<-fullBB.lst[[1]][[1]]
  if(length(fullBB.lst)>1)
    {
      BB.lst[[ 2 ]]<-fullBB.lst[[2]][[1]]
      if(length(fullBB.lst)>2)
        {
          for(rr in 3:length(fullBB.lst))
            {
              for(kk in 1:length(fullBB.lst[[rr]]))
                {
                  BB.lst[[ cumlen[rr-1]+kk ]]<-fullBB.lst[[rr]][[kk]]
                }
            }
        }
    }
  return(BB.lst)
}


  
  

  
get.allpossROUTE.fun<-function()
{
  ### when this file is modified, function "getPIstar.fun" need to be changed as well!!!!!!
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

#WCC getBigMatB.fun<-function(BB.lst=BB.cs, subPI.lst=subPI.cs){
getBigMatB.fun<-function(BB.lst, subPI.lst){
   ## modified on Apr 6, 2005; now check for the compatibility of dimensionality before adding or mulitplying	
   ## see ****** part. functions used:  chkdim.Multiply.fun, and  chkdim.Add.fun
    mmat         <-get.allpossROUTE.fun()
    r.ele        <-c( 1,1,2,4,8,16,32,64,128,256)
    cumsum.r     <-cumsum(r.ele)
    lenPI        <-length(subPI.lst)  # 4
    lenBB        <-length(BB.lst)     # 7=1+1+2+3
    lenROUTE     <-cumsum.r[lenPI]    # 8=1+1+2+4
    RRstar       <-length(subPI.lst[[1]]) ##*********
    
    ## changed on 03/11/05 (remove "if" condition)
    mat.route<-as.matrix(mmat[(1:lenROUTE),(1:lenBB)])
    pi.ele   <-rep((1:10), r.ele)
    j.index  <-pi.ele[1:cumsum.r[lenPI]] ## indexing a mat B is for the level of H^{(j)}

    preBigBB.lst<-list()
    for(kk in 1:lenROUTE)## it's fine to have kk start from "1"; when lenROUTE=1, BB=bb[[1]]=Identity Matrix
      {
        bb<-BB.lst[mat.route[kk,]]
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
        preBigBB.lst[[kk]]<-BB
      }
    
    indx<-split((1:length(preBigBB.lst)), j.index)
    BigMatB<-NULL
    for(jj in 1:length(indx))
      {
        tmpB<-preBigBB.lst[ indx[[jj]] ]
        if(length(tmpB)==1)
          {
            BigMatB<-rbind(BigMatB, tmpB[[1]])
          }else{
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


##   getBigMatB.fun<-function(BB.lst=BB.cs, subPI.lst=subPI.cs)
##     {
##       mmat         <-get.allpossROUTE.fun()
##       r.ele        <-c( 1,1,2,4,8,16,32,64,128,256)
##       cumsum.r     <-cumsum(r.ele)
##       lenPI        <-length(subPI.lst)  # 4
##       lenBB        <-length(BB.lst)     # 7=1+1+2+3
##       lenROUTE     <-cumsum.r[lenPI]    # 8=1+1+2+4
##   
##       ## changed on 03/11/05 (remove "if" condition)
##       mat.route<-as.matrix(mmat[(1:lenROUTE),(1:lenBB)])
##       pi.ele   <-rep((1:10), r.ele)
##       j.index  <-pi.ele[1:cumsum.r[lenPI]] ## indexing a mat B is for the level of H^{(j)}
##   
##       preBigBB.lst<-list()
##       for(kk in 1:lenROUTE)
##         {
##           bb<-BB.lst[mat.route[kk,]]
##   	## in "getPIstar.fun" there is such if, but it seems to be useless so i removed
##   	## on 03/10/05
##           ## if(sum(is.na(bb))>0){
##           ##   ##          prePIstar.lst[[kk]]<-rep(0, length(subPI.lst[[1]]))
##           ##   cat("dd=",dd); print(" is.na(bb)>0")
##           ## }else{ below...}
##           BB<-bb[[length(bb)]]
##           for(jj in max(1,(length(bb)-1)):1){BB<-BB%*%bb[[jj]]}
##           preBigBB.lst[[kk]]<-BB
##         }
##       indx<-split((1:length(preBigBB.lst)), j.index)
##       BigMatB<-NULL
##       for(jj in 1:length(indx))
##         {
##           tmpB<-preBigBB.lst[ indx[[jj]] ]
##           if(length(tmpB)==1)
##             {
##               BigMatB<-rbind(BigMatB, tmpB[[1]])
##             }else{
##               sumB<-0
##               for(ll in 1:length(tmpB)){sumB<-sumB+tmpB[[ll]]}
##               BigMatB<-rbind(BigMatB, sumB)
##             }
##         }
##       return(BigMatB)
##     }
##   
##   
##   



##rm(BB.lst, subPI.lst, prePIstar.lst, mmat, r.ele, cumsum.r, lenPI, lenBB, lenROUTE, mat.route, pi.ele, PI,bb, BB, kk, jj,PIstar)
#WCC getPIstar.fun<-function(BB.lst=BB.cs, subPI.lst=subPI.cs)
getPIstar.fun<-function(BB.lst, subPI.lst)
  {
    mmat         <-get.allpossROUTE.fun()
    r.ele        <-c( 1,1,2,4,8,16,32,64,128,256)
    cumsum.r     <-cumsum(r.ele)
    lenPI        <-length(subPI.lst)  # 4
    lenBB        <-length(BB.lst)     # 7=1+1+2+3
    lenROUTE     <-cumsum.r[lenPI] # 8=1+1+2+4

    ## changed on 03/11/05 (remove "if" condition)
    mat.route<-as.matrix(mmat[(1:lenROUTE),(1:lenBB)])
    pi.ele   <-rep((1:10), r.ele)
    PI       <-subPI.lst[pi.ele[1:cumsum.r[lenPI]]]

    prePIstar.lst<-list()
    for(kk in 1:lenROUTE)
      {
        bb<-BB.lst[mat.route[kk,]]
        if(sum(is.na(bb))>0){
          prePIstar.lst[[kk]]<-rep(0, length(subPI.lst[[1]]))
        }else{
          BB<-bb[[length(bb)]]
          for(jj in max(1,(length(bb)-1)):1){BB<-BB%*%bb[[jj]]}
          prePIstar.lst[[kk]]<-PI[[kk]] %*% BB
        }
      }

    PIstar<-colSums(matrix(unlist(prePIstar.lst), nrow=length(prePIstar.lst), byrow=T))
    names(PIstar)<-names(subPI.lst[[1]])
    return(PIstar)
  }

