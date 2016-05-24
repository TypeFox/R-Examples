## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work
## Author: Andrew Hooker

line_search <- function(f_name, f_options, l, u, x, f, g, d, options=list()){
  #determine maximum step size
  
  if(!isfield(options, 'ftol')){
    options$ftol=1e-3
  }
  if(!isfield(options, 'gtol')){
    options$gtol=0.9
  }
  if(!isfield(options, 'xtol')){
    options$xtol=0.1
  }
  
  xtol=options$xtol
  gtol=options$gtol
  ftol=options$ftol
  
  fixed=(x<=l | x>=u)
  stpmx=Inf
  temp1=Inf
  for(i in find_matlab(!fixed)){
    dk=d[i]
    if(dk<0){
      temp2=l[i]-x[i]
      if(temp2>=0){
        temp1=0
      } else {
        if(dk*stpmx < temp2) temp1=temp2/dk
      } 
    } else {
      temp2=u[i]-x[i]
      if(temp2<=0){
        temp1=0
      } else if(!is.na(dk*stpmx)) {
        if(dk*stpmx>temp2) temp1=temp2/dk
      }
    }
    stpmx = min(temp1,stpmx)
  }
  stp=min(1/norm(d,type="2"),stpmx)
  stpmin=0
  ##calc directional derivative
  gp=t(g)%*%d
  
  xtrapl=1.1
  xtrapu=4
  bracket=FALSE
  stage=1
  finit=f
  ginit=gp
  gtest=ftol*ginit
  width=stpmx-stpmin
  width1=width/0.5
  stx=0
  fx=finit
  gx=ginit
  sty=0
  fy=finit
  gy=ginit
  stmin=0
  stmax=stp + xtrapu*stp
  while(TRUE){
    f_options[[1]] <- x+stp*d
    returnArgs <- do.call(f_name,f_options)
    f <- returnArgs[[1]]
    g <- returnArgs[[2]]
    gp=t(g)%*%d
    ftest = finit+stp*gtest
    if(stage==1 && f<=ftest && gp>=0){
      stage=2
    }
    if(f<=ftest && abs(gp) <=gtol*(-ginit)){
      break
    }
    
    if((bracket && (stp <= stmin || stp>=stmax))){
      break
    }
    if((bracket && stmax-stmin < xtol*stmax)){
      break
    }
    if((stp == stpmx && f <= ftest && gp <= gtest) ){
      break
    }
    if((stp == stpmin && (f > ftest || gp >= gtest))){
      break
    }
    if(stage==1 && f<=fx && f>ftest ){
      fm=f-stp*gtest
      fxm=fx-stx*gtest
      fym=fy-sty*gtest
      gm=gp-gtest
      gxm=gx-gtest
      gym=gy-gtest
      returnArgs <- compute_step(stx,fxm,gxm,sty,fym,gym,stp,fm,gm,bracket,stmin,stmax) 
      stx <- returnArgs[[1]]
      fx <- returnArgs[[2]]
      dx <- returnArgs[[3]]
      sty <- returnArgs[[4]]
      fy <- returnArgs[[5]]
      dy <- returnArgs[[6]]
      stp <- returnArgs[[7]]
      bracket <- returnArgs[[8]]
      fx=fxm+stx*gtest
      fy=fym+sty*gtest
      gx=gxm+gtest
      gy=gym+gtest
    } else {
      returnArgs <- compute_step(stx,fx,gx,sty,fy,gy,stp,f,gp,bracket,stmin,stmax) 
      stx <- returnArgs[[1]]
      fx <- returnArgs[[2]]
      dx <- returnArgs[[3]]
      sty <- returnArgs[[4]]
      fy <- returnArgs[[5]]
      dy <- returnArgs[[6]]
      stp <- returnArgs[[7]]
      bracket <- returnArgs[[8]]
    }
    #Decide if a bisection step is needed.
    if((bracket) ){
      if(abs(sty-stx)>0.66*width1){
        stp=stx+0.5*(sty-stx)
      }
      width1=width
      width=abs(sty-stx)
    }
    if((bracket) ){
      stmin=min(stx,sty)
      stmax=max(stx,sty)
    } else {
      stmin=stp+xtrapl*(stp-stx)
      stmax=stp+xtrapu*(stp-stx)
    }  
    stp=max(stp,stpmin)
    stp=min(stp,stpmx)
    if(bracket && (stp<=stmin || stp>=stmax) || (bracket && stmax-stmin < xtol*stmax)){
      stp=stx
    }
    
  }
  x=x+stp*d
  return(list( f= f,g=g,x=x)) 
}

compute_step <- function(stx,fx,dx,sty,fy,dy,stp,fp,dp,bracket,stpmin,stpmax){
  sgnd=dp*(dx/abs(dx))
  #First case: A higher function value. The minimum is bracketed. 
  #If the cubic step is closer to stx than the quadratic step, the 
  #cubic step is taken, otherwise the average of the cubic and 
  #quadratic steps is taken.
  if(fp>fx ){
    theta=3*(fx-fp)/(stp-stx)+dx+dp
    s=max(matrix(c(abs(theta),abs(dx),abs(dp)),nrow=1,byrow=T))
    gamma=s*sqrt((theta/s)^2-(dx/s)*(dp/s))
    if(stp<stx){
      gamma=-gamma
    }
    p=(gamma-dx)+theta
    q=((gamma-dx)+gamma)+dp
    r=p/q
    stpc=stx+r*(stp-stx)
    stpq=stx + ((dx/((fx-fp)/(stp-stx)+dx))/2)*(stp-stx)
    if(abs(stpc-stx)<abs(stpq-stx)){
      stpf=stpc
    } else {
      stpf=stpc+(stpq-stpc)/2
    }
    bracket=TRUE
    #Second case: A lower function value and derivatives of opposite 
    #sign. The minimum is bracketed. If the cubic step is farther from
    #stp than the secant step, the cubic step is taken, otherwise the
    #secant step is taken.
  } else if(sgnd<0) {
    
    theta = 3*(fx - fp)/(stp - stx) + dx + dp
    s = max(matrix(c(abs(theta),abs(dx),abs(dp)),nrow=1,byrow=T))
    gamma = s*sqrt((theta/s)^2 - (dx/s)*(dp/s))
    if((stp > stx) ){
      gamma = -gamma
    }
    p = (gamma - dp) + theta
    q = ((gamma - dp) + gamma) + dx
    r = p/q
    stpc = stp + r*(stx - stp)
    stpq = stp + (dp/(dp - dx))*(stx - stp)
    if(abs(stpc-stp) > abs(stpq-stp)){
      stpf = stpc
    } else {
      stpf = stpq
    }
    bracket=TRUE
    #Third case: A lower function value, derivatives of the same sign,
    # and the magnitude of the derivative decreases.
  } else if (abs(dp) < abs(dx)){ 
    #The cubic step is computed only if the cubic tends to infinity
    #in the direction of the step or if the minimum of the cubic
    #is beyond stp. Otherwise the cubic step is defined to be the
    #secant step.
    theta = 3*(fx - fp)/(stp - stx) + dx + dp
    s = max(matrix(c(abs(theta),abs(dx),abs(dp)),nrow=1,byrow=T))
    
    #The case gamma = 0 only arises if the cubic does not tend
    #to infinity in the direction of the step.
    
    gamma = s*sqrt(max(0,(theta/s)^2-(dx/s)*(dp/s)))
    if((stp>stx) ){
      gamma = -gamma
    }
    p = (gamma - dp) + theta
    q = (gamma + (dx - dp)) + gamma
    r = p/q
    if((r <0 && gamma!=0)){
      stpc = stp + r*(stx - stp)
    } else if(stp > stx){
      stpc = stpmax
    } else {
      stpc = stpmin
    }
    stpq = stp + (dp/(dp - dx))*(stx - stp)
    
    if(bracket){
      
      #A minimizer has been bracketed. If the cubic step is
      #closer to stp than the secant step, the cubic step is
      #taken, otherwise the secant step is taken.
      
      if(abs(stpc-stp)<abs(stpq-stp)){
        stpf = stpc
      } else {
        stpf = stpq
      }
      if((stp > stx)){
        stpf = min(stp+0.66*(sty-stp),stpf)
      } else {
        stpf = max(stp+0.66*(sty-stp),stpf)
      }
    } else {
      #A minimizer has not been bracketed. If the cubic step is
      #farther from stp than the secant step, the cubic step is
      #taken, otherwise the secant step is taken.
      
      if(abs(stpc-stp) > abs(stpq-stp)){
        stpf = stpc
      } else {
        stpf = stpq
      }
      stpf = min(stpmax,stpf)
      stpf = max(stpmin,stpf)
    }
    #Fourth case: A lower function value, derivatives of the same sign, 
    #and the magnitude of the derivative does not decrease. If the
    #minimum is not bracketed, the step is either stpmin or stpmax,
    #otherwise the cubic step is taken.
  } else {
    if(bracket){
      theta = 3*(fp - fy)/(sty - stp) + dy + dp
      s = max(matrix(c(abs(theta),abs(dy),abs(dp)),nrow=1,byrow=T))
      gamma = s*sqrt((theta/s)^2 - (dy/s)*(dp/s))
      if(stp > sty){
        gamma = -gamma
      }
      p = (gamma - dp) + theta
      q = ((gamma - dp) + gamma) + dy
      r = p/q
      stpc = stp + r*(sty - stp)
      stpf = stpc
    } else if (stp > stx){
      stpf = stpmax
    } else {
      stpf = stpmin
    }
  }
  if((fp >fx)){
    sty=stp
    fy=fp
    dy=dp
  } else {
    if(sgnd <0){
      sty = stx
      fy = fx
      dy = dx
    }
    stx = stp
    fx = fp
    dx = dp
  }
  stp=stpf
  return(list( stx= stx,fx=fx,dx=dx,sty=sty,fy=fy,dy=dy,stp=stp,bracket=bracket)) 
}
