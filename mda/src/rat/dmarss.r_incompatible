subroutine marss(nx,n,p,nclass,y,x,w,tagx,maxorder,mmax,penalty,thresh,forwstep,interms,prune,bx,fullin,lenb, bestgcv, bestin, flag, cut,dir, res,alpha,beta,
scrat,iscrat,trace)
#
# multivariate mars
#
# input 
#
# integer nx- number of rows of x, bx, and tagx
# integer n- number of observations
# integer p- number of variables
# integer nclass - number of classes (resposne variables)
# double precision y(n,nclass)- matrix of response values
# double precision x(nx,p) -matrix of predictors
# double precision w(n)- weights (ge 0)
# integer tagx(nx,p) - tag array for sorting x
# integer maxorder- maxumum order of interaction (mi in Friedman's notation)
#
# integer mmax- maximum number of basis functions (nk in  Friedman's notation)
#     if it is desired to add say 3 terms to an existing model, one can set mmax=interms+3
#     and set interms appropriately (see below)
#
# double precision penalty- cost per knot selection (df in Friedman's notation)
#
# double precision thresh- forward stepwise threshold - forw step stops   if (rss/rssnull).le.thresh)
#    to prevent numerical problems with overfit solutions
#
# logical forwstep- if true, forward stepwise is carried out. if false, procedure
#   starts with basis matrix bx and carried out backward stepwise. This saves
#   a great deal of time in choosing the penalty parameter 
#
# integer interms- number of terms in existing model. for a fresh start interms=1.
#
# logical prune- true if pruning (backward stepwise deletion) is desired
#
#
#
# 
# 
# double precision bx(nx,mmax)- basis matrix- used on input only if interms>1 or forwstep=false
#
# integer fullin(mmax)- basis column indicator. columns marked with a 1 indicate active
#    columns-  used on input only if interms>1 or forwstep=false
#    
# integer lenb- number of active columns of bx -
#      used on input only if interms>1 or forwstep=false
#
# double precision res(n,nclass) - residuals - used as input only if interms>1, in which case they
# should be the residuals from the current fitted model
#
# double precision scrat(1)- scratch array of length at least 
# 1+n+ 2*n*mmax +4*mmax*mmax +3*mmax + 3*mmax*nclass +3*nclass
#                                               +28n+51 (for linpack)
#
# integer iscrat(1)- scratch array of length at least 4*mmax
#
#logical trace - should there be a printing of progress
#
#
# output
#
# double precision bx(nx,mmax)- basis matrix. only first lenb columns are active. first column is 1
# integer fullin(mmax)- basis column indicator. columns marked with a 1 indicate
#    the linearly indep set 
# integer lenb- number of active columns of bx
# double precision bestgcv- gcv for best model
# integer bestin(mmax) basis column indicator for best model. 1 indicates term is in model
# integer flag(mmax,p)- ijth element is 1 if predictor j is in model term i 
#    (zero otherwise). note for a zero column of bx (fullin=0), the row of flag
#    will be all 1, for a technical reason, it should be ignored of course.
#
# double precision cut(mmax,p) - ijth element cutpoint for variable j in model term i
# double precision dir(mmax,p) - ijth element is 1 if factor involving predictor j in term i is
#  of the form (x-t)+ ; it is -1 if it has the form (t-x)+
#
# double precision res(n,nclass) - residuals from final fit
# double precision alpha(nclass)- intercepts from final fit
# double precision beta(mmax,nclass)- slopes from final fit



implicit double precision (a-h,o-z)

integer nx, n,p,nclass,tagx(nx,p),maxorder,mmax,bestin(mmax),flag(mmax,p),
fullin(mmax)
double precision y(n,nclass),x(nx,p),w(n),bx(nx,mmax),bestgcv,cut(mmax,p),dir(mmax,p),res(nx,nclass),
alpha(nclass),beta(mmax,nclass)

double precision scrat(1)
integer iscrat(1)
logical forwstep, prune, trace, tracec
common tracec
tracec=trace



# iscrat must be of len at least 3*mmax
# scrat must be of length at least 
#  1+n+ 2*n*mmax +4*mmax*mmax +3*mmax + 3*mmax*nclass +3*nclass

#1. bxorth
#2. bxorthm
#3. cov
#4. covsy
#5. ybar
#6. scr1
#7. scr5
#8. scr6
#9.  temp

#10. bxsc
#11. r
#12. betasc,
#13
#14
#15  work

len1=n*mmax
len2=mmax
len3=mmax*mmax
len4=mmax*nclass
len5=nclass
len6=mmax
len7=mmax
len8=nclass
len9=n


len10=n*mmax
len11=mmax*mmax
len12=mmax*nclass
len13=mmax*mmax
len14=mmax*mmax

n1=1
n2=n1+len1
n3=n2+len2
n4=n3+len3
n5=n4+len4
n6=n5+len5
n7=n6+len6
n8=n7+len7
n9=n8+len8
n10=n9+len9
n11=n10+len10
n12=n11+len11
n13=n12+len12
n14=n13+len13
n15=n14+len14




call marsnew1(nx,n,p,nclass,y,x,w,tagx,maxorder,mmax,bx,bestgcv, bestin,fullin,lenb,flag, cut,dir,res,
alpha,beta,penalty,thresh,forwstep,interms,prune,
scrat,scrat(n2),scrat(n3),scrat(n4),scrat(n5),
scrat(n6),scrat(n7),scrat(n8),scrat(n9),
scrat(n10),scrat(n11),scrat(n12),scrat(n13),scrat(n14),scrat(n15),
iscrat,iscrat(1+mmax),iscrat(1+2*mmax),iscrat(1+3*mmax))

return
end

subroutine marsnew1(nx,n,p,nclass,y,x,w,tagx,maxorder,mmax,bx,bestgcv, bestin,fullin,lenb,flag, cut,dir,res,alpha,beta,penalty,thresh,forwstep,interms,prune,
bxorth,bxorthm,cov,covsy,ybar,
scr1,scr5,scr6, temp, 
  bxsc, r, betasc,varsc,var,work,
termlen,in, tempin, qpivot)
# input n,p,nclass,x,tagx,y,w,mmax,maxorder

# nclass is # of response variables

# output- bestin, bestgcv, flag, cut,dir, bx,res

implicit double precision (a-h,o-z)

integer n,nterms2,p,mmax,flag(mmax,p),tagx(nx,p),termlen(mmax), nclass,fullin(mmax)

double precision cov(mmax,mmax),covsy(mmax,nclass),
critmax,x(nx,p),bx(nx,mmax),bxorth(n,mmax),bxorthm(mmax),y(n,nclass),ybar(nclass),
scr1(mmax),scr5(mmax),scr6(nclass)

double precision temp(n),w(n), cut(mmax,p),dir(mmax,p),alpha(nclass),beta(mmax,nclass),  bxsc(n,mmax), r(mmax,mmax), dofit,
  res(nx,nclass),betasc(mmax,nclass), varsc(mmax,mmax), var(mmax,mmax), stopfac,work(1)
  
integer   tempin(mmax), bestin(mmax),qrank, qpivot(mmax)

logical forwstep,go, prune, newform, cvar, trace
common trace

#
# forw stepwise stops if gcv/gcvnull > stopfac, to avoid adding lots of
# noisy terms
#
double precision rtemp(4)
integer itemp(4)

tolbx=.01
stopfac=10.0
#stopfac=10e9
prevcrit=10e9
#write(6,*) "in marsnew mmax=",mmax, " forwstep=",forwstep, "interms=",interms

#call intpr("M1",2,n,1)
#call dblepr("penalty",7,penalty,1)
#call intpr("maxo",4,maxorder,1)

if(interms.eq.1) {dofit=0}
 else{dofit=0;do j=2,lenb {dofit=dofit+fullin(j)}
     nterms=interms
}


if(forwstep){
fullin(1)=1
do i=2,mmax{
 fullin(i)=0
}

do i=1,n {
  w(i)=1
}

do i=1, mmax{
  termlen(i)=0
  do j=1, p{
   flag(i,j)=0
   cut(i,j)=0
}}



nterms=1;  nterms2=2;
do i=1,n {
   bx(i,1)=1
  bxorth(i,1)=1.0/dsqrt(dfloat(n))
}

bxorthm(1)=1/dsqrt(dfloat(n))

do i=1,n {
 do j=1, mmax{
   bx(i,j)=0.0
}}

do i=1,n {
   bx(i,1)=1
}

do k=1, nclass{
 ybar(k)=0.0
 do i=1,n {
   ybar(k)=ybar(k)+y(i,k)/n
}}
if(interms.eq.1){
rssnull=0.0
 do k=1, nclass{ 
 do i=1,n {
   rssnull=rssnull+(y(i,k)-ybar(k))**2
 }}
}
else{
 rssnull=0.0
 do k=1, nclass{ 
 do i=1,n {
   rssnull=rssnull+res(i,k)**2
 }}
 
}
 


rss=rssnull

cmm= (1+dofit) + penalty*(.5*dofit)
gcvnull=(rssnull/n)/(1.0-cmm/n)**2

#write(6,*) "initial rss=", rssnull, " initial gcv=", gcvnull
if(trace)call dblepr("initial rss=",11,rssnull,1)
if(trace)call dblepr("initial gcv=",11,gcvnull,1)

lenb=1
ii=interms-1
go=.true.
while( (ii.lt.(mmax-1)).and.((rss/rssnull).gt.thresh).and.go){
  ii=ii+2

#  write(6,*) "bef addtrm",ii, nterms

   do i1=1, nterms{
     do i2=1, nterms{
     cov(i1,i2)=0
  }}

  do j=1, nterms{
    cov(j,j)=0.0
    do i=1,n {
      cov(j,j)=cov(j,j)+(bxorth(i,j)-bxorthm(j))*(bxorth(i,j)-bxorthm(j))
    }
  }
  


  
 do k=1,nclass{
   
   do j=1, nterms{
     covsy(j,k)=0.0
     do i=1,n {
      covsy(j,k)=covsy(j,k)+(y(i,k)-ybar(k))*bxorth(i,j)
      }
    }
  }




  do ik=1,mmax { tempin(ik)=fullin(ik)}

  call addtrm(nx,bx,tempin,bxorth,bxorthm,p,n,nclass,rss,prevcrit,
     cov,covsy,y,ybar,x,tagx,w,termlen,mmax,tolbx,
  nterms,flag,maxorder,scr1,scr5,scr6,imax,jmax,kmax,critmax, newform,bxsc, r, betasc, temp)
 
# check whether to accept term
doftemp=dofit
 doftemp=doftemp+1
if((imax>1).and.(newform)) {doftemp=doftemp+1}




temprss=rss-critmax
cmm= (1+doftemp) + penalty*(.5*doftemp)
gcv=(temprss/n)/(1.0-cmm/n)**2

go=.false.
if(((critmax/rss).gt.thresh).and.((gcv/gcvnull).lt.stopfac)){
   go=.true.
   dofit=doftemp

    rss=rss-critmax
    kk=tagx(imax,jmax)
#write(6,256) jmax, imax ,kmax, critmax,x(kk,jmax), rss,gcv, dofit
256 format(" ","adding term"," jmax=",i3, "  imax=",i3 ,"  kmax=",i3, "  critmax= ",f8.2,
"  cutp=", f9.5," rss=",f8.2, " gcv=",f8.2, " dofit=",f9.3)
#    write(6,*) 
itemp(1)=jmax;itemp(2)=imax;itemp(3)=kmax;
rtemp(1)=critmax;rtemp(2)=x(kk,jmax);rtemp(3)=rss;rtemp(4)=gcv

if(trace) call intpr("adding term ",12,ii,1)
if(trace)call intpr("var, sp index, parent",21,itemp,3)
if(trace)call dblepr("critmax cut rss gcv",19,rtemp,4)

prevcrit=critmax
  do j=1,p {
    flag(ii,j)=flag(kmax,j)
   flag(ii+1,j)=flag(kmax,j)
   cut(ii,j)=cut(kmax,j)
   cut(ii+1,j)=cut(kmax,j)
   dir(ii,j)=dir(kmax,j)
   dir(ii+1,j)=dir(kmax,j)
   

  }
  termlen(ii)=termlen(kmax)+1
  termlen(ii+1)=termlen(kmax)+1


 


 do i=1,n { temp(i)=x(tagx(i,jmax),jmax)}
 temp1=temp(imax)

  fullin(ii)=1;

  if((imax.gt.1).and.(newform)) {fullin(ii+1)=1 }

        flag(ii,jmax)=1;
       flag(ii+1,jmax)=1;
       cut(ii,jmax)=temp1
       cut(ii+1,jmax)=temp1
        dir(ii,jmax)=1
        dir(ii+1,jmax)=-1
#
#this is to prevent trying to split a 0 column later

 if(fullin(ii+1).eq.0)
  {
#do j=1,p {flag(ii+1,j)=1}
#   
       termlen(ii+1)=maxorder+1
  }


  
  do i=1,n { if( (x(i,jmax)-temp1).gt.0) bx(i,ii)=bx(i,kmax)*(x(i,jmax)-temp1)
              if((temp1-x(i,jmax)).ge.0) bx(i,ii+1)=bx(i,kmax)*(temp1-x(i,jmax))
           }
 


 
    if(nterms.eq.1){
      temp1=0.0;  do i=1,n { temp1=temp1+bx(i,2)/n; } 
            do i=1,n {bxorth(i,2)=bx(i,2)-temp1;}
    }
    else{
##     write(6,*) "bef call to orthreg"


 call orthreg(n,n,nterms,bxorth,fullin, bx(1,ii),bxorth(1,nterms2))
          }

 if(fullin(ii+1).eq.1){


 call orthreg(n,n,nterms+1,bxorth,fullin, bx(1,ii+1),bxorth(1,nterms2+1))
#
  }
  else {do i=1,n {bxorth(i,nterms2+1)=0}}

  bxorthm(nterms2)=0.0 ; 
  bxorthm(nterms2+1)=0.0

   do i=1,n {
     bxorthm(nterms2)=bxorthm(nterms2)+bxorth(i,nterms2)/n
     bxorthm(nterms2+1)=bxorthm(nterms2+1)+bxorth(i,nterms2+1)/n
   }

 

   temp1=0.0;temp2=0.0;do i=1,n {temp1=temp1+bxorth(i,nterms2)**2
                                 temp2=temp2+bxorth(i,nterms2+1)**2 }
##write(6,*) "norm of new columns",temp1,temp2
   if (temp1.gt.0.0) {do i=1,n { bxorth(i,nterms2)  =bxorth(i,nterms2)/dsqrt(temp1) }}
   if (temp2.gt.0.0) {do i=1,n { bxorth(i,nterms2+1)=bxorth(i,nterms2+1)/dsqrt(temp2) }}
   lenb=lenb+2
  nterms=nterms+2; nterms2=nterms2+2;

}}

#write(6,*) "stopping forw stepwise"
#if((rss/rssnull).le.thresh) write(6,*) "rss ratio=",rss/rssnull
#if((critmax/rss).le.thresh) write(6,*) "crit ratio=",critmax/rss
#if((gcv/gcvnull).gt.stopfac) write(6,*) "gcv ratio=",gcv/gcvnull

rtemp(1)=rss/rssnull;rtemp(2)=critmax/rss;rtemp(3)=gcv/gcvnull
if(trace)call dblepr("stopping forw step; rss crit and gcv ratios",43,rtemp,3)
if(trace){
 if((rss/rssnull).le.thresh) call dblepr("rss ratio=",10,rss/rssnull,1)
 if((critmax/rss).le.thresh) call dblepr ("crit ratio=",11,critmax/rss,1)
call dblepr("critmax",7,critmax,1)
call dblepr("rss",3,rss,1)
 if((gcv/gcvnull).gt.stopfac) call dblepr("gcv ratio=",10,gcv/gcvnull,1)
}

# write(6,*) "after forward step  ","go=",go," rss=",rss, "nterms=", nterms


}


#call intpr("M2",2,fullin,10)
#call intpr("nterms",6,nterms,1)
dofit= -1
do i=1,nterms{
 bestin(i)=fullin(i)
 dofit=dofit+fullin(i)
}
#call intpr("M21",3,n,1)


if(trace)call intpr("aft forw step",13,nterms,1)

##added by Trevor 5/24/2015
cvar=.false.
call  qrreg(nx,n,mmax,lenb,nclass,bx,bxsc,bestin,y,qpivot,qrank,beta,res,rss,cvar,var,varsc,scr1, work) 


#call dblepr("rss",3,rss,1)
#call dblepr("beta",4,beta,4)
#call intpr("M3",2,n,1)

nt=dofit+1
if(qrank< nt) {
#write(6,*) "singular matrix"
 do i=qrank+1,nt{ 
  bestin(qpivot(i))=0
  fullin(qpivot(i))=0
#write(6,*) "removing col ",qpivot(i)
  dofit=dofit-1;
} }


##write(6,*) "bef mlin", nterms,dofit,(bestin(j),j=1,nterms),y(1,1),n,mmax,nclass,w(1)

cvar=.true.
#call intpr("M23",3,fullin,11)
#call intpr("bestin",6,bestin,11)
#
#call intpr("lenb",4,lenb,1)
#call intpr("mmax",4,mmax,1)
#call dblepr("res",3,res,10)
#call dblepr("bx",2,bx(1,1),5)
#call dblepr("bx",2,bx(1,2),5)
#call dblepr("bx",2,bx(n,11),1)
#call dblepr("y",1,y,5)
#call intpr("nx",2,nx,1)
#call intpr("n",1,n,1)




#write(6,*) "rank=",qrank, "rss=",rss


 





rssfull=rss

# dofit is the number of nonconstant, indep basis fns in the model



 cmm= (1+dofit) + penalty*(.5*dofit)
bestgcv=(rss/n)/(1.0-cmm/n)**2

#write(6,*) "full model,gcv=",bestgcv, " rss=",rssfull, "dofit=",dofit
rtemp(1)=bestgcv; rtemp(2)=rssfull;rtemp(3)=dofit
if(trace)call dblepr("full model: gcv rss dofit",25,rtemp,3)
if(trace)call intpr("terms",5,fullin,lenb)

# after forward stepwise, nterms and lenb = # of columns of bx
# dofit is number of nonzero columns 



if(prune){
 
do i=1,mmax{ tempin(i)=bestin(i)}
#write(6,*) "initial tempin", (tempin(i),i=1,nterms)

#call intpr("M5",2,n,1)
while(dofit>0 ) {
    
     
    jo=1;rsstemp=10e99;minterm=0
     do ii=2, lenb {

       if(tempin(ii).eq.1) {
        jo=jo+1
          temp7=0.0
          do kc=1,nclass{
           temp7=temp7+beta(jo,kc)**2/var(jo,jo)
          }

        if(temp7 < rsstemp) {minterm=ii; rsstemp=temp7}
      


      }}
         rss=rss+rsstemp
        dofit=dofit-1; 
        cmm= (1.0+dofit) + penalty*(.5*dofit)
         gcv=(rss/n)/(1.0-cmm/n)**2
        
         tempin(minterm)=0
      

#         write(6,100) minterm,gcv, rss,dofit,(tempin(ik),ik=1,lenb)
100 format(" ","pruning, minterm= ",i4, " gcv=",f9.3,2x, " rss=",f9.3,2x," dof=",f9.3,
" model= ",60(i1,1x))
   
 
    if(gcv< bestgcv) {bestgcv=gcv;do i=1,mmax {bestin(i)=tempin(i);}}


  if(dofit > 0) {  
cvar=.true.
call  qrreg(nx,n,mmax,lenb,nclass,bx,bxsc,tempin,y,qpivot,qrank,beta,res,rss,cvar,var,varsc,scr1,work)


          }
}



call  qrreg(nx,n,mmax,lenb,nclass,bx,bxsc,bestin,y,qpivot,qrank,beta,res,rss,cvar,var,varsc,
scr1, work)
#write(6,101) bestgcv,rss,(bestin(ik),ik=1,lenb)
101 format(" ","best model gcv=",f9.3," rss=",f9.3,2x,"model= ",60(i1,1x))
if(trace)call intpr("best model",10,bestin,lenb)
if(trace)call dblepr(" gcv=",4,bestgcv,1)


}





return
end
subroutine addtrm(nx,bx,tempin,bxorth,bxorthm,p,n,nclass,
rss,prevcrit,cov,covsy,y,ybar,x,tagx,w,termlen,mmax,tolbx,
nterms,flag,maxorder,scr1,scr5,scr6,imax,jmax,kmax,critmax, newform,bxsc,r,betasc, scrat)


implicit double precision (a-h,o-z)

integer n,nterms,nterms2,p,mmax,flag(mmax,p),v,tagx(nx,p),termlen(mmax), nclass, tempin(mmax), minspan, iendspan

double precision cov(mmax,mmax),covsy(mmax,nclass),
critmax,x(nx,p),bx(nx,mmax),bxorth(n,mmax),bxorthm(mmax),y(n,nclass),ybar(nclass),
scr1(mmax),scr5(mmax),scr6(nclass), bxsc(n,mmax), r(mmax,mmax),
 betasc(mmax,nclass), scrat(n),w(n)

double precision temp1, temp2, scr2,sumb, sumbx, su, st, tem
logical newform, tnewform, trace
common trace
critmax=0;jmax=0;imax=0;kmax=0;


# finds best term to add to current mars model







##write(6,*) "in addtrm n,p, nclass, nterms,mmax, y(1,1), x(1,1), tagx(1,1),flag(1,1),bx(1,1#),bxorth(1,1)", n,p, nclass, nterms,mmax, y(1,1), x(1,1), tagx(1,1),flag(1,1),bx(1,1),bxor#th(1,1)





do m=1,nterms {

 nm=0
 do jjj=1,n {
  if(bx(jjj,m).gt.0) {nm=nm+1}
 }

# tem=-(1d0/(n*nm))*dlog(1d0 - 5d-2)  n should be p thanks to  Gints Jekabsons/Stephen Milborrow
 tem=-(1d0/(p*nm))*dlog(1d0 - 5d-2)

 minspan= -1d0*(dlog(tem)/dlog(2d0))/2.5

# tem=(5d-2)/n   n should be p thanks to  Gints Jekabsons/Stephen Milborrow
 tem=(5d-2)/p
 iendspan=3d0-dlog(tem)/dlog(2d0)

# # write(6,*) "minspan,iendspan",minspan,iendspan

if(termlen(m)< maxorder){

 do v=1,p {

  if(flag(m,v).eq.0){
# check if model term type is already in model, to avoid a linear dependence


  tnewform=.true.
  mm=1
  while((mm.le.nterms).and.tnewform) {
     mm=mm+1
    if(tempin(mm).eq.1){
     tnewform=.false.
     if(flag(mm,v).ne.1){tnewform=.true.;go to 9911}
     do j=1,p {if(j.ne.v) { if(flag(mm,j).ne.flag(m,j)) {tnewform=.true.;go to 9911 }}}
    }
9911 continue
}

# if new form of term, fit term bm*x
# should it be bxorth below?

if(tnewform) {

nterms2=nterms+1

     do i=1,n {
        scrat(i)=x(i,v)*bx(i,m)
      }
#
#

if(nterms>1){

 call orthreg(n,n,nterms,bxorth,tempin, scrat,bxorth(1,nterms2))
   }
else{
   tem=0;do i=1,n {tem=tem+scrat(i)/n;}
    do i=1,n{ bxorth(i,2)=scrat(i)-tem}
   }


    bxorthm(nterms2)=0.0 ; 
   do i=1,n {
     bxorthm(nterms2)=bxorthm(nterms2)+bxorth(i,nterms2)/n 
   }
   temp1=0.0;do i=1,n {temp1=temp1+bxorth(i,nterms2)**2}                               
   if (temp1.gt.tolbx) {do i=1,n { bxorth(i,nterms2)=bxorth(i,nterms2)/dsqrt(temp1) }}
     else {do i=1,n {bxorth(i,nterms2)=0};tnewform=.false.;}

    do i1=1, nterms2{
         cov(i1,nterms2)=0.0
         cov(nterms2, i1)=0.0
    }

 
      cov(nterms2,nterms2)=1


  do kc=1,nclass{
     covsy(nterms2,kc)=0.0
     do i=1,n {
      covsy(nterms2,kc)=covsy(nterms2,kc)+(y(i,kc)-ybar(kc))*bxorth(i,nterms2)

      }
   }


    critnew=0.0
    do kc=1,nclass {
      temp1=0
      do i=1,n { temp1=temp1+y(i,kc)*bxorth(i,nterms2)}
      critnew=critnew+temp1**2
     }
## write(6,*) "geoff",nterms,v,critnew
    if(critnew.gt.critmax) {
       jmax=v
       critmax=critnew
       imax=1
       kmax=m
      }

}


if(tnewform) {nterms2=nterms+1; nterms21=nterms+2} else{nterms2=nterms; nterms21=nterms+1;critnew=0.0}

#if((nterms.eq.7).and.(v.eq.2)) {
# call dblepr("bxorth",6,bxorth(1,nterms2),10)
# call dblepr("temp1",5,temp1,1)
# call dblepr("critmax",7,critmax,1)
 
#}





# try other knot locations





     do kc=1, nclass{
      covsy(nterms21,kc)=0
     }
     do ii=1,nterms21 {
        cov(ii,nterms21)=0
        cov(nterms21,ii)=0
     }
   do kc=1,nclass {
     scr6(kc)=0
   }
    do ii=1,nterms21{
        scr1(ii)=0;
    }
scr2=0;su=0;st=0;sumbx2=0;sumb=0.0;sumbx=0.0;
   for (k=n-1;k>0;k=k-1) {
 
   do i=1,nterms2 {
       kk=tagx(k,v)
       kk1=tagx(k+1,v)
       scr1(i)=scr1(i)+(bxorth(kk1,i)-bxorthm(i))*bx(kk1,m)
       cov(i,nterms21)=cov(i,nterms21)+ (x(kk1,v)-x(kk,v))*scr1(i)
       cov(nterms21,i)=cov(i,nterms21)   
      }
     scr2=scr2+(bx(kk1,m)**2)*x(kk1,v)

     sumbx2=sumbx2+bx(kk1,m)**2
    sumb=sumb+bx(kk1,m)
    sumbx=sumbx+bx(kk1,m)*x(kk1,v)
     su=st
    st=sumbx-sumb*x(kk,v)
     cov(nterms21,nterms21)= cov(nterms21,nterms21)+ 
      (x(kk1,v)-x(kk,v))*(2*scr2-sumbx2*(x(kk,v)+x(kk1,v)))+ 
          ( (su*su)-(st*st) )/n


   crittemp=critnew
     do kc=1, nclass{

     scr6(kc)=scr6(kc)+(y(kk1,kc)-ybar(kc))*bx(kk1,m)
     covsy(nterms21,kc)=covsy(nterms21,kc )+(x(kk1,v)-x(kk,v))*scr6(kc)

    temp1=covsy(nterms21,kc)
    temp2=cov(nterms21,nterms21)

     do jk=1,nterms2 {
        temp1=temp1-covsy(jk,kc)*cov(jk,nterms21)
         temp2=temp2-cov(jk,nterms21)*cov(jk,nterms21)
       }
#
# this has to be fixed!!
#



 if(cov(nterms21,nterms21)>0){
        if((temp2/cov(nterms21,nterms21)) > tolbx) {critadd=(temp1*temp1)/temp2} else {critadd=0.0}}
 else{critadd=0}

 crittemp=crittemp+critadd

#if((nterms.eq.7).and.(v.eq.2).and.(critadd.gt.0)) {

#call intpr("k",1,k,1)
 #call dblepr("temp2",5,temp2,1)
#call dblepr("cov",3,cov(nterms21,nterms21),1)
# call dblepr("critadd",7,critadd,1)
# call dblepr("crittemp",8,crittemp,1)

#}



 
 if(crittemp.gt.(1.01*rss)) {crittemp=0.0}
 if(crittemp.gt.(2*prevcrit)) {crittemp=0.0}

 }
         
   


if(k.gt.1){
  k0=tagx(k-1,v)
}
#
# decide whether to accept term
#

#if((k.eq.70).and.(v.eq.10)){
##write(6,*) "HI",k,v,critmax,minspan,mod(k,minspan),iendspan,n-iendspan,x(kk,v),x(kk0,v)
#}
#call intpr("k",1,k,1)
#call intpr("v",1,v,1)
#call dblepr("crittemp",8,crittemp,1)
#call intpr("minspan",7,minspan,1)
#call intpr("iendspan",8,iendspan,1)
     if((crittemp.gt.critmax).and.(mod(k,minspan).eq.0).and.
(k.ge.iendspan).and.(k.le.(n-iendspan)).and.(bx(kk1,m).gt.0).and.
(.not.((k.gt.1).and.(x(kk,v).eq.x(k0,v))) )) {
       jmax=v
       critmax=crittemp
       imax=k
       kmax=m
       newform=tnewform

       
      }
   


   } }
9999 continue
 }}


}



return
end




