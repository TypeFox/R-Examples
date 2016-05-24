# This is a special version of K-dO rq modified to compute the LR process
# Note that the sol array is now p+3 by J and the third row contains the
# value of the objective function at the tau specified in row one.
#
#     This software is in the public domain  and may be freely
#     used and redistributed for non-commercial purposes.  No guarantees
#     are offered or implied.  comments, bug reports, etc are welcome
#     and should be sent to roger@ysidro.econ.uiuc.edu or to
#
#        roger koenker
#        department of economics
#        university of illinois
#        champaign, illinois, 61820
#
#
subroutine rqbr(m,nn,m5,n3,n4,a,b,t,toler,ift,x,e,s,wa,wb,nsol,ndsol,sol,dsol,lsol,h,qn,cutoff,ci,tnmat,big,lci1)
#
#     m = number of observations
#     n = number of parameters
#     m5 = m+5
#     n3 = n+3
#     n4 = n+4
#     a is the x matrix
#     b is the y vector
#     t, the desired quantile
#     toler, smallest detectable |x-y|/x machine precision to the 2/3
#     ift exit code:
#     		0-ok
#    		else dimensions inconsistent see below
#     x the parameter estimate betahat
#     e is the residual vector
#     s is an integer work array
#     wa is a real work array
#     wb is another real work array
#     nsol is an estimated (row) dimension of the primal solution array say 3*m; minimum value = 2
#     ndsol is an estimated (row) dimension of the dual solution array say 3*m; minimum value = 2
#     sol is the primal solution array
#     dsol is the dual solution arry
#     lsol is the actual dimension of the solution arrays
#     h is the matrix of basic observations indices
#     qn is the vector of residuals variances from the projection of each column of the x matrix on the remaining columns
#     cutoff, the critical point for N(0,1)
#     ci is the matrix of confidence intervals
#     tnmat is the matrix of the JGPK rank test statistics
#     big, large positive finite floating-point number
#     utilization:  if you just want a solution at a single quantile you
#     neednt bother with sol, nsol, etc, if you want all the solutions
#     then set theta to something <0 and sol and dsol will return all the
#     estimated quantile solutions.
#     the algorithm  is a slightly modified version of algorithm as 229
#     described in koenker and dorey, computing regression quantiles,
#     applied statistics, pp. 383-393.
#
integer i,j,k,kl,kount,kr,l,lsol,m,m1,m2,m3,m4,m5,ift
integer n,n1,n2,n3,nsol,ndsol,out,s(m),h(nn,nsol)
integer nn,n4,idxcf
logical stage,test,init,iend,lup
logical lci1,lci2,skip
double precision a1,aux,b1,big,d,dif,pivot,smax,t,t0,t1,tnt
double precision min,max,toler,zero,half,one,two
double precision b(m),sol(n3,nsol),a(m,nn),x(nn),wa(m5,n4),wb(m)
double precision sum,e(m),dsol(m,ndsol)
double precision qn(nn),cutoff,ci(4,nn),tnmat(4,nn),tnew,told,tn
parameter( zero = 0.d0)
parameter( one = 1.d0)
parameter( two = 2.d0)
#
#  check dimension parameters
#
n=nn
ift = 0
wa(m+2,nn+1) = one
if (m5!=m+5)
  ift = 3
if (n3!=n+3)
  ift = 4
if (n4!=n+4)
  ift = 5
if (m<=zero||n<=zero)
  ift = 6
if (ift<=two) {
#
#  initialization
#
  half = one/two
  iend = .true.
  lci2 = .false.
  lup = .true.
  skip = .false.
  idxcf = 0
  tnew = zero
  tn = zero
  m1 = m+1
  n1 = n+1
  n2 = n+2
  m2 = m+2
  m3 = m+3
  m4 = m+4
  do j = 1,n{
    x(j) = zero
    }
  do i = 1,m
    e(i) = zero
  if (t<zero||t>one) {
    t0 = one/(two*float(m))
    t1 = one-t0
    t = t0
    iend = .false.
    lci1 = .false.
    }
  repeat { #0
    do i = 1,m {
      k = 1
      do j = 1,nn
        if (k<=nn){
          if(j==idxcf)
            skip = .true.
          else
            skip = .false.
          if(!skip){
            wa(i,k) = a(i,j)
            k = k+1
            }
          }  
      wa(i,n4) = n+i
      wa(i,n2) = b(i)
      if (idxcf != 0)
        wa(i,n3) = tnew*a(i,idxcf)
      else
        wa(i,n3) = zero
      wa(i,n1) = wa(i,n2)-wa(i,n3)
      if (wa(i,n1)<zero)
        do j = 1,n4
          wa(i,j) = -wa(i,j)
      }
    do j = 1,n {
      wa(m4,j) = j
      wa(m2,j) = zero
      wa(m3,j) = zero
      do i = 1,m {
        aux = sign(one,wa(m4,j))*wa(i,j)
        wa(m2,j) = wa(m2,j)+aux*(one-sign(one,wa(i,n4)))
        wa(m3,j) = wa(m3,j)+aux*sign(one,wa(i,n4))
        }
      wa(m3,j) = two*wa(m3,j)
      }
    dif = zero
    init = .false.
    if(!lci2){
#compute the columns means
      do k = 1,n {
        wa(m5,k) = zero
        do i = 1,m{
          wa(m5,k) = wa(m5,k)+a(i,k)
          }
        wa(m5,k) = wa(m5,k)/float(m)
        }
      }
    lsol = 1
    kount = 0
    repeat { #1
#
# compute new marginal costs
#
      do j = 1,n{
        wa(m1,j) = wa(m2,j)+wa(m3,j)*t
        }
      if (!init) {
#
# stage 1
#
# determine the vector to enter the basis
#
        stage = .true.
        kr = 1
        kl = 1
        go to 30
        }
      repeat { #2
#
# stage 2
#
        stage = .false.
        repeat { #3
#
# determine the vector to enter the basis
#
          max = -big
          do j = kr,n {
            d = wa(m1,j)
            if (d<zero) {
              if (d>(-two))
                next 1
              d = -d-two
              }
            if (d>max) {
              max = d
              in = j
              }
            }
          if (max<=toler)
            break 2
          if (wa(m1,in)<=zero) {
            do i = 1,m4
              wa(i,in) = -wa(i,in)
            wa(m1,in) = wa(m1,in)-two
            wa(m2,in) = wa(m2,in)-two
            }
          repeat { #4
#
# determine the vector to leave the basis
#
            k = 0
            do i = kl,m {
              d = wa(i,in)
              if (d>toler) {
                k = k+1
                wb(k) = wa(i,n1)/d
                s(k) = i
                test = .true.
                }
              }
            repeat { #5
              if (k<=0)
                test = .false.
              else {
                min = big
                do i = 1,k
                  if (wb(i)<min) {
                    j = i
                    min = wb(i)
                    out = s(i)
                    }
                wb(j) = wb(k)
                s(j) = s(k)
                k = k-1
                }
#
# check for linear dependence in stage 1
#
              if (!test&&stage)
                break 1
              if (!test)
                break 5
              pivot = wa(out,in)
              if (wa(m1,in)-pivot-pivot<=toler){
                go to 10
                }
              do j = kr,n3 {
                d = wa(out,j)
                wa(m1,j) = wa(m1,j)-d-d
                wa(m2,j) = wa(m2,j)-d-d
                wa(out,j) = -d
                }
              wa(out,n4) = -wa(out,n4)
              } #5
            do i = 1,m4 {
              d = wa(i,kr)
              wa(i,kr) = wa(i,in)
              wa(i,in) = d
              }
            kr = kr+1
            go to 20
#
# pivot on wa(out,in)
#
            10  do j = kr,n3
              if (j!=in)
                wa(out,j) = wa(out,j)/pivot
            do i = 1,m3
              if (i!=out) {
                d = wa(i,in)
                do j = kr,n3
                  if (j!=in)
                    wa(i,j) = wa(i,j)-d*wa(out,j)
                }
            do i = 1,m3
              if (i!=out)
                wa(i,in) = -wa(i,in)/pivot
            wa(out,in) = one/pivot
            d = wa(out,n4)
            wa(out,n4) = wa(m4,in)
            wa(m4,in) = d
            kount = kount+1
            if (!stage)
              break 1
#
# interchange rows in stage 1
#
            kl = kl+1
            do j = kr,n4 {
              d = wa(out,j)
              wa(out,j) = wa(kount,j)
              wa(kount,j) = d
              }
            20  if (kount+kr==n1)
              break 2
            30  max = -one
            do j = kr,n
              if (abs(wa(m4,j))<=n) {
                d = abs(wa(m1,j))
                if (d>max) {
                  max = d
                  in = j
                  }
                }
            if (wa(m1,in)<zero)
              do i = 1,m4
                wa(i,in) = -wa(i,in)
            } #4
          } #3
        } #2
      if (kr==1) {
        do j = 1,n {
          d = abs(wa(m1,j))
          if (d<=toler||two-d<=toler){
            ift = 1
            wa(m2,nn+1) = zero
            go to 80
            }
          }
        }
      80 kount = 0
      sum = zero
      if (!lci2){
        do i = 1,kl-1 {
          k = wa(i,n4)*sign(one,wa(i,n4))
          x(k) = wa(i,n1)*sign(one,wa(i,n4))
          }
        }
      do i = 1,n {
        kd = abs(wa(m4,i))-n
        dsol(kd,lsol) = one+wa(m1,i)/two
        if (wa(m4,i)<zero)
          dsol(kd,lsol) = one-dsol(kd,lsol)
        if (!lci2){
          sum = sum + x(i)*wa(m5,i)
          sol(i+3,lsol) = x(i)
          h(i,lsol) = kd
          }
        }
      do i = kl,m {
        kd = abs(wa(i,n4))-n
        if (wa(i,n4)<zero)
          dsol(kd,lsol) = zero
        if (wa(i,n4)>zero)
          dsol(kd,lsol) = one
        }
      if (!lci2){
        sol(1,lsol) = smax
        sol(2,lsol) = sum
        sum = zero
        do j=kl,m{
          d = wa(j,n1)*sign(one,wa(j,n4))
          sum = sum + d*(smax + half*(sign(one,d) - one))
	 }
        sol(3,lsol) = sum
        do i=1,m
          dsol(i,lsol+1) = dsol(i,lsol)
        }
      if (lci2){
# compute next theta
        a1 = zero
        do i = 1,m {
          a1 = a1+a(i,idxcf)*(dsol(i,lsol)+t-one)
          }
        tn = a1/sqrt(qn(idxcf)*t*(one-t))
        if (abs(tn)<cutoff){
          if (lup)
            smax = big
          else
            smax = -big
          do i =1,kl-1 {
            k = wa(i,n4)*sign(one,wa(i,n4))
            sol(k,1) = wa(i,n2)*sign(one,wa(i,n4))
            sol(k,2) = wa(i,n3)*sign(one,wa(i,n4))/tnew
            }
          do i = kl,m {
            a1 = zero
            b1 = zero
            k = wa(i,n4)*sign(one,wa(i,n4))-n
            l = 1
            do j = 1,n{
              if (j==idxcf)
                l = l+1
              a1 = a1 + a(k,l)*sol(j,1)
              b1 = b1 + a(k,l)*sol(j,2)
              l = l+1
              }
            tnt = (b(k)-a1)/(a(k,idxcf)-b1)
            if (lup){
              if (tnt>tnew)
                if (tnt<smax){
                  smax = tnt
                  out = i
                  }
              }
            else {
              if (tnt<tnew)
                if (tnt>smax){
                  smax = tnt
                  out = i
                  }
              }
            }
          if (lup){
            told = tnew 
            tnew = smax + toler
            ci(3,idxcf) = told - toler
            tnmat(3,idxcf) = tn
            if (!(tnew < big-toler)){
              ci(3,idxcf) = big
              ci(4,idxcf) = big
              tnmat(3,idxcf) = tn
              tnmat(4,idxcf) = tn
              lup = .false.
              go to 70
              }
            }
          else{
            told = tnew 
            tnew = smax - toler
            ci(2,idxcf) = told + toler
            tnmat(2,idxcf) = tn
            if (!(tnew > -big+toler)){
              ci(2,idxcf) = -big
              ci(1,idxcf) = -big
              tnmat(2,idxcf) = tn
              tnmat(1,idxcf) = tn
              lup = .true.
              go to 60
              }
            }
#update the new marginal cost
          do i = 1,m{
            wa(i,n3) = wa(i,n3)/told*tnew
            wa(i,n1) = wa(i,n2) - wa(i,n3)
            }
          do j = kr,n3{
            d = wa(out,j)
            wa(m1,j) = wa(m1,j) -d -d
            wa(m2,j) = wa(m2,j) -d -d
            wa(out,j) = -d
            }
          wa(out,n4) = -wa(out,n4)
          init = .true.
          }
        else{
          if (lup){
            ci(4,idxcf) = tnew - toler
            tnmat(4,idxcf) = tn
            lup = .false.
            go to 70
            }
          else{
            ci(1,idxcf) = tnew + toler
            tnmat(1,idxcf) = tn
            lup = .true.
            go to 60
            }
          }
        }
      if ((iend)&&(!lci2))
        go to 40
      if (!lci2){
        init = .true.
        lsol = lsol+1
        do i = 1,m
          s(i) = zero
        do j = 1,n
          x(j) = zero
#
#  compute next t
#
        smax = two
        do j = 1,n {
          b1 = wa(m3,j)
          a1 = (-two-wa(m2,j))/b1
          b1 = -wa(m2,j)/b1
          if (a1>=t)
            if (a1<smax) {
              smax = a1
              dif = (b1-a1)/two
              }
          if (b1>t)
            if (b1<smax) {
              smax = b1
              dif = (b1-a1)/two
              }
          }
        tnt = smax+toler*(one+abs(dif))
        if (tnt>=t1+toler)
          iend = .true.
        t = tnt
        if (iend)
          t = t1
        }
      } #1
    wa(m2,nn+1) = two
    ift = 2
    go to 50
    40  if (lsol>2) {
      sol(1,1) = zero
      #sol(2,1) = zero
      sol(3,1) = zero
      sol(1,lsol) = one
      #sol(2,lsol) = zero
      sol(3,lsol) = zero
      do i = 1,m {
        dsol(i,1) = one
        dsol(i,lsol) = zero
        dsol(i,lsol+1) = zero
        }
      }
    l = kl-1
    do i = 1,l
      if (wa(i,n1)<zero)
        do j = kr,n4
          wa(i,j) = -wa(i,j)
    50  sum = zero
    if(!lci2){
      do i = kl,m {
        k = wa(i,n4)*sign(one,wa(i,n4))
        d = wa(i,n1)*sign(one,wa(i,n4))
        sum = sum+d*sign(one,d)*(half+sign(one,d)*(t-half))
        k = k-n
        e(k) = d
        }
      wa(m2,n2) = kount
      wa(m1,n2) = n1-kr
      wa(m1,n1) = sum
      }
    if (wa(m2,nn+1)==two)
      break 1
    if (!lci1)
      break 1
    if (!lci2){
      lci2 = .true.
      n = nn-1
      n1 = n+1
      n2 = n+2
      n3 = n+3
      n4 = n+4
    60 idxcf = idxcf+1
      if (idxcf>nn){
        break 1
	}
    70  if (lup){
        tnew = x(idxcf)+toler
        told = tnew
        ci(3,idxcf) = x(idxcf)
        tnmat(3,idxcf) = zero
        }
      else{
        tnew = x(idxcf)-toler
        told = tnew
        ci(2,idxcf) = x(idxcf)
        tnmat(2,idxcf) = zero
        }
      }
    } #0
  # restore the original value of dual when ci is true
  do i=1,m
    dsol(i,lsol) = dsol(i,lsol+1)
  }
return
end
