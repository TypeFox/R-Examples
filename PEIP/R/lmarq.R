lmarq <-
function(afun,ajac,p0,tol,maxiter)
{ 
  Vnorm <-function(X)
    {
      return( sqrt(sum(X^2)) )
    }


###   Initialize p and oldp.

 func =  get(afun)
  jacob =  get(ajac)

 
p=p0
fp=Vnorm(func(p))^2
oldp=p0*2
oldfp=fp*2

###  Initialize lambda.
lambda=0.0001


###  The main loop.  While the current solution isn't good enough, keep
###  trying...  Stop after maxiter iterations in the worst case.
### 
iter=0
while(iter <= maxiter)
  {
  ###  Compute the Jacobian.
  J=jacob(p)

  ###  Compute rhs=-J'*f
  rhs=  - t(J) %*% func(p)

  ###  Check the termination criteria.
  ###  all criteria are relative measires
  ### 
  ###  the current rhs must be small enough we won't move much AND
  ###  the change in the norm of f must not have changed to much last step AND
  ###  the point can't have moved too far last step
  if ((Vnorm(rhs)< sqrt(tol)*(1+abs(fp))) &(abs(oldfp-fp)<tol*(1+abs(fp))) & (Vnorm(oldp-p)<sqrt(tol)*(1+Vnorm(p))))
    {
    pstar=p
    return(list(pstar=pstar,iter=iter))
     } 

  ###  We use a clever trick here.  The least squares problem
  ### 
  ###   min || [ J              ] s - [ -F ] ||
  ###       || [ sqrt(lambda)*I ]     [ 0  ] ||
  ### 
  ###  Has the normal equations solution
  ### 
  ###   s=-inv(J'*J+lambda*I)*J'*F
  ### 
  ###  which is precisely the LM step.  We can solve this least squares problem
  ###  more accurately using the QR factorization then by computing 
  ###  inv(J'*J+lambda*I) explicitly.
  ###


  myrhs=c( -func(p),  rep(0,length(p))  )


  s= Ainv( rbind(J, sqrt(lambda)*diag(1,length(p)) )  , myrhs )



  ###  compute the new chisq value
  f1 = func(p+s)
  fpnew=Vnorm(f1)^2

  ###  If the chisq improves then f is improved, make the step, and decrease lambda
  if(fpnew < fp)
  {
    oldp=p
    oldfp=fp
    p=p+s
    fp=fpnew
    lambda=lambda/2
    if (lambda <10^(-12))
      {
      lambda=1.0e-12
       } #
  }
  else
    {
    ###  Didn't improve f, increase lambda, and try again.
    lambda=lambda*2.5
    if (lambda >10^16)
      {
      lambda=10^16
       } #
     } #

  ### update the iteration count
  iter=iter+1
   } #

###   Return, max iters exceeded.
pstar=p


 return(list(pstar=pstar,iter=iter)) 

}
