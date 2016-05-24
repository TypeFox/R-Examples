#Mills lower bound on the integral of a standard normal over an interval
#Need to use something else near zero, where mills is bad.
mills.lb = function(a,b){
  t=5#threshold for switching to approximation
  if(b<a){return(0)}
  if(a>=t){
    return(a*exp(-a^2/2)/(1+a^2)-exp(-b^2/2)/b) 
  }
  if(b <= -t){
    return(mills.lb(-b,-a))
  }
  #Note, like in the rest of this program, I'm not dividing by sqrt(2*pi), so we need to rescale here
  sqrt(2*pi)*pnorm(min(b,t))-sqrt(2*pi)*pnorm(max(a,-t))+mills.lb(a,-t)+mills.lb(t,b)
}


#Truncate the interval.  We chop off the ends near infinity, being careful 
#so that the chopped tail is guaranteed to be close enough to its mills approximation
#a,b are the left and right endpoints, z is the mid point
#delta is the multiplicative error limit of the truncation on the final fraction (roughly)
truncate.interval = function(a,b,z,delta=1e-16){
  #Initialize some stuff
  L.extra = 0#Extra probability to be added for the truncation on the left
  R.extra = 0#Extra probability to be added for the truncation on the right
  a.new = a#truncated interval bounds
  b.new = b#truncated interval bounds
  
  #We need bounds on the integrals 
  RL.lb = mills.lb(a,b)
  R.lb = mills.lb(z,b)
  
  #Now we bound the error we can tolerate in the tail approximations
  eps.R = min(delta*R.lb,delta*RL.lb/2)
  eps.L = delta*RL.lb/2
  
  #For now, only truncate infinite end points
  #Might want to change this one day, if we have trouble with super wide but finite intervals
  if (b==Inf){
    f = function(x){x^2+log(1+x^2)+log(eps.R)}#encodes error of mills approximation
    b.new = uniroot(f,c(1.1,1000))$root
    b.new = max(b.new,z+1)#Don't truncate past z
    R.extra = exp(-b.new^2/2)/b.new
  }
  if (a==-Inf){
    f = function(x){x^2+log(1+x^2)+log(eps.L)}#encodes error of mills approximation
    a.new = -uniroot(f,c(1.1,1000))$root
    a.new = min(a.new,z-1)#Don't truncate past z
    L.extra = exp(-a.new^2/2)/a.new
  }
  
  list(a=a.new,b=b.new,L.extra=L.extra,R.extra=R.extra,z=z)
}


#Approximates integral_a^b e^{-x^2/2+offset^2/2} dx
# offset is used to make ratios slightly more stable
# defaults to offset=0
# Note that I've left out 1/sqrt(2*pi), you can add it in if you like
approx.int = function(a,b,n=1000,offset=0){
  delta = (b-a)/n #Step size, may want to vary in the future
  x = seq(from=a,to=b,by=delta)
  y = -x^2/2 + offset^2/2 # On the log scale
  m = diff(y)/diff(x) # Line segment slopes
  de = diff(exp(y)) # Difference on original scale
  sum(de/m) #Sum of integrals of line segments (closed form)
}

#Uses approx.int to evaluate int_x^b phi(z)dz / int_a^b phi(z)dz
#Right now offsets everything for a little more stability
#Uses truncation to handle infinite endpoints
max.approx.frac = function(a,b,x,mu=0,n=1000){
  returns = numeric(length(mu))
  for(i in 1:length(returns)){
    truncation = truncate.interval(a-mu[i],b-mu[i],x-mu[i])
    #Our offset will use the smaller of a and b in absolute value
    offset = min(abs(truncation$a),abs(truncation$b))
    #The truncation also shifts by the mean, so we don't need to do it again for the end points
    #but we do need to use the center z returned by truncation, rather than x, to match
    left = approx.int(truncation$a,truncation$z,n,offset)+truncation$L.extra
    right = approx.int(truncation$z,truncation$b,n,offset)+truncation$R.extra
    returns[i] = right/(left+right)
  }
  returns
}
