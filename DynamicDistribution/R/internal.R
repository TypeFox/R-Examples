#'@keywords internal
#'@export 
#'Normal
#'Chi_Square
#'Exponential 
#'Rayleigh
#'F_dis
#'Hypergeometric
#'Con_Uniform
#'Student_t
#'Gamma_dis
#'Beta_dis
#'Laplace
#'Logistic
#'Lognormal
#'Pareto
#'Cauchy
#'Inverse_Gaussian
#'Geometric
#'Dis_Uniform
#'Bernoulli
#'Binomial
#'Poisson
#'Negative_Binomial
#'Logarithmic_Series

Chi_Square<-function(x=NULL,para)
{
  n<-para
  n<-seq(round(min(n)),round(max(n)),1)
  n=as.matrix(n)
  if(all(n<=0))
    stop("Error! The parameter n must be a positive integer!")
  
  if(!is.numeric(x))
  { 
    if(all(n<=2))
      stop("In order to calculate Mode, the maximum value of the parameter n must be larger than 2!")
    
    n_Mode<-subset(n,n>2)
    Mode_seq=n_Mode-2    
    y_max<-max(exp(-Mode_seq/2)*(Mode_seq)^(n_Mode/2-1)/(2^(n_Mode/2))/gamma(n_Mode/2))
    return(c(0.0001,2.5*max(Mode_seq),y_max))
  }
  if(is.numeric(x))
  {
    density<-ifelse(x>=0,x^((n-2)/2)*exp(-x/2)/(2^(n/2)*gamma(n/2)),0)
    fun<-function(x){x^((n-2)/2)*exp(-x/2)/(2^(n/2)*gamma(n/2))}
    cdf<-c()
    for(i in x)
    {
      cdf<-c(cdf,as.numeric(ifelse(i>0,integrate(fun,0.0001,i),0))) 
    }
    Mean<-n
    Variance<-2*n
    if(n>=2)
      Mode=n-2
    else
      Mode=Inf
    Skewness<-(8/n)^(1/2)
    Kurtosis<-3+12/n
    return(list(density=density,
                cdf=cdf,
                Mean=Mean,
                Variance=Variance,
                Mode=Mode,
                Skewness=Skewness,
                Kurtosis=Kurtosis,
                n=n,
                'n',
                "Chi-Square Distribution"))
  }
}

Exponential<-function(x=NULL,para)
{
  b<-para
  if(any(b<=0))
    stop("Error! The parameter must be positive integer!")
  if(!is.numeric(x))
  {
    b=subset(b,select=(b>0))
    b_mode<-b
    y_max<-max(exp(-b_mode/b)/b)
    return(c(0,2*max(b_mode),y_max))
  } 
  
  if(is.numeric(x))
  {
    density<-ifelse(x>=0,(exp(-x/b))/b,0)
    fun<-function(x){(exp(-x/b))/b}
    cdf<-c()
    for(i in x)
    {
      cdf<-c(cdf,as.numeric(ifelse(i>=0,integrate(fun,0,i),0)))
    }
    
    Mean<-b
    Variance<-b^2
    Mode<-0
    Skewness<-2
    Kurtosis<-9
    
    return(list(density=density,
                cdf=cdf,
                Mean=Mean,
                Variance=Variance,
                Mode=Mode,
                Skewness=Skewness,
                Kurtosis=Kurtosis,
                b=b,
                'b',
                "Exponential Distribution"))
  }
}
Rayleigh<-function(x=NULL,para)
{
  b<-para
  if(any(b<=0))
    stop("Error! The parameter must be positive integer!")
  
  if(!is.numeric(x))
  {
    b<-subset(b,select=(b>0))
    Mode_seq=b
    y_max<-max(Mode_seq*exp(-Mode_seq^2/(2*b^2))/b^2)
    return(c(0,2.5*max(Mode_seq),y_max))
  }
  if(is.numeric(x))
  {
    density<-ifelse(x>=0,x*exp(-x^2/(2*b^2))/b^2,0)
    fun<-function(x){x*exp(-x^2/(2*b^2))/b^2}
    cdf<-c()
    for(i in x)
    {
      cdf<-c(cdf,as.numeric(ifelse(i>=0,integrate(fun,0,i),0))) 
    }
    Mean<-b*sqrt(pi/2)
    Variance<-(2-pi/2)*b^2
    Mode<-b
    Skewness<-2*(pi-3)*sqrt(pi)/(4-pi)^(3/2)
    Kurtosis<-(32-3*pi^2)/(4-pi)^2
    
    return(list(density=density,
                cdf=cdf,
                Mean=Mean,
                Variance=Variance,
                Mode=Mode,
                Skewness=Skewness,
                Kurtosis=Kurtosis,
                b=b,
                'b',
                "Rayleigh Distribution"))
  }
}

Student_t<-function(x=NULL,para)
{
  n<-para
  n<-seq(round(min(n)),round(max(n)),1)
  n=as.matrix(n)
  if(any(n<=0))
    stop("Error! The parameter must be positive integer!")
  
  if(!is.numeric(x))
  {
    if(all(n<=2))
      stop("In order to calculate Variance, the maximum value of the parameter n must larger than 2!")
    
    n_Variance<-subset(n,n>2)
    StaVariance_seq=sqrt(as.vector(n_Variance/(n_Variance-2)))
    y_max<-max(gamma((n+1)/2)/gamma(n/2)/sqrt(n*pi))
    return(c(-3*max(StaVariance_seq),3*max(StaVariance_seq),y_max))
  }
  if(is.numeric(x))
  {
    density<-(gamma((n+1)/2))/((gamma(n/2))*(sqrt(n*pi))*((1+x^2/n)^((n+1)/2)))
    fun<-function(x){(gamma((n+1)/2))/((gamma(n/2))*(sqrt(n*pi))*((1+x^2/n)^((n+1)/2)))}
    cdf<-c()
    for(i in x)
    {
      cdf<-c(cdf,integrate(fun,-Inf,i)$value) 
    }
    if(n>1)
      Mean<-0
    else
      Mean<-Inf
    if(n>2)
      Variance<-n/(n-2)
    else
      Variance<-Inf
    Mode<-0
    Skewness<-0
    if(n>4)
      Kurtosis<-(3*(n-2))/(n-4)
    else
      Kurtosis<-Inf
    
    return(list(density=density,
                cdf=cdf,
                Mean=Mean,
                Variance=Variance,
                Mode=Mode,
                Skewness=Skewness,
                Kurtosis=Kurtosis,
                n=n,
                'n',
                "Student's t Distribution"))
  }
}

Beta_dis<-function(x=NULL,para,const_par)
{
  a<-para[1,];b<-para[2,]
  const_par=const_par
  if(any(a<=0)||any(b<=0))
    stop("Error! The parameter must be positive integer!")
  if(const_par[1]<=0||const_par[2]<=0)
    stop("Error! The constant parameter must be positive integer!")
  
  if(!is.numeric(x))
  {
    if(all(a<=1)||all(b<=1))
      stop("Error! In order to calculate the mode, the largest number of a,b must larger than 2!")
    if(const_par[1]<=1||const_par[2]<=1)
      stop("Error! In order to calculate the mode, constant a,b must larger than 2!")
    
    a_Mode=subset(a,a>1)
    b_Mode=subset(b,b>1)
    Mode_seq1<-(a_Mode-1)/(a_Mode+const_par[2]-2)
    y_max1<-max((Mode_seq1^(a_Mode-1))*((1-Mode_seq1)^(const_par[2]-1))/beta(a_Mode,const_par[2]))
    Mode_seq2<-(const_par[1]-1)/(const_par[1]+b_Mode-2)
    y_max2<-max((Mode_seq2^(const_par[1]-1))*((1-Mode_seq2)^(b_Mode-1))/beta(const_par[1],b_Mode))
    return(c(0,1,y_max1,0,1,y_max2))
  }
  
  if(is.numeric(x))
  {
    
    density<-(x^(a-1))*((1-x)^(b-1))/beta(a,b)
    fun<-function(x){(x^(a-1))*((1-x)^(b-1))/beta(a,b)}
    cdf<-c()
    for(i in x)
    {
      cdf<-c(cdf,as.numeric(ifelse(i>=0,integrate(fun,0,i),0)))
    }
    
    Mean<-a/(a+b)
    Variance<-a*b/(((a+b)^2)*(a+b+1))
    Mode<-ifelse(a>1&b>1,(a-1)/(a+b-2),Inf)
    Skewness<-2*(b-a)*((a+b+1)^0.5)/((a+b+2)*((a*b)^0.5))
    Kurtosis<-3*(a+b+1)*(2*((a+b)^2)+a*b*(a+b-6))/(a*b*(a+b+2)*(a+b+3))
    return(list(density=density,cdf=cdf,Mean=Mean,Variance=Variance,Mode=Mode,Skewness=Skewness,Kurtosis=Kurtosis,a=a,'a',b=b,'b',"Beta Distribution"))
  }
}

Cauchy<-function(x=NULL,para,const_par)
{
  a<-para[1,];b<-para[2,]
  const_par=const_par
  if(any(b<=0))
    stop("Error! The parameter b must be positive number!")
  if(const_par[2]<=0)
    stop("Error! The constant parameter b must be positive number!")
  
  if(!is.numeric(x))
  {
    a=a
    b=subset(b,b>0)
    Mode_seq1<-a
    y_max1<-max(1/(pi*const_par[2]*(1+((Mode_seq1-a)/const_par[2])^2)))
    Mode_seq2<-const_par[1]
    y_max2<-max(1/(pi*b*(1+((Mode_seq2-const_par[1])/b)^2)))
    return(c(-4*const_par[2]+min(a),4*const_par[2]+max(a),y_max1,-4*max(b)+const_par[1],4*max(b)+const_par[1],y_max2))
  }
  
  if(is.numeric(x))
  {
    
    density<-1/(pi*b*(1+((x-a)/b)^2))
    fun<-function(x){1/(pi*b*(1+((x-a)/b)^2))}
    cdf<-c()
    for(i in x)
    {
      cdf<-c(cdf,integrate(fun,-Inf,i)$value)
    }
    
    Mean<-Inf
    Variance<-Inf
    Mode<-a
    Skewness<-Inf
    Kurtosis<-Inf
    return(list(density=density,cdf=cdf,Mean=Mean,Variance=Variance,Mode=Mode,Skewness=Skewness,Kurtosis=Kurtosis,a=a,'a',b=b,'b',"Cauchy Distribution"))
  }
}

Con_Uniform<-function(x=NULL,para,const_par)
{
  a<-para[1,];b<-para[2,]
  const_par=const_par
  if(max(a)>=const_par[2]||const_par>=min(b))
    stop("Error! Parameter a should be less than b")
  
  if(!is.numeric(x))
  {
    y_max1<-1.2*max(1/(const_par[2]-a))
    y_max2<-1.2*max(1/(b-const_par[1]))
    return(c(min(a),const_par[2],y_max1,const_par[1],max(b),y_max2))
  }
  
  if(is.numeric(x))
  {
    density<-ifelse(x<a,0,ifelse(x>b,0,1/(b-a)))
    fun<-function(x){1/(b-a)}
    cdf<-c()
    for(i in x)
    {
      cdf<-c(cdf,as.numeric((i-a)/(b-a)))
    }
    
    Mean<-(b+a)/2
    Variance<-((b-a)^2)/12
    Mode<-Inf
    Skewness<-0
    Kurtosis<-Inf
    
    return(list(density=density,cdf=cdf,Mean=Mean,Variance=Variance,Mode=Mode,Skewness=Skewness,Kurtosis=Kurtosis,a=a,'a',b=b,'b',"Continuous Uniform Distribution"))
  }
}

F_dis<-function(x=NULL,para,const_par)
{
  m<-para[1,];n<-para[2,]
  const_par=const_par
  if(any(m<=0)||any(n<=0))
    stop("Error! The parameter must be positive integer!")
  if(const_par[1]<=0||const_par[2]<=0)
    stop("Error! The constant parameter must be positive integer!")
  
  if(!is.numeric(x))
  {
    if(all(m<=2))
      stop("Error! In order to calculate the mode, the largest number of m must larger than 2!")
    if(const_par[1]<=2)
      stop("Error! In order to calculate the mode, constant m must larger than 2!")
    
    m_Mode=subset(m,m>2)
    n=subset(n,n>0)
    Mode_seq1<-const_par[2]*(m_Mode-2)/m_Mode/(const_par[2]+2)
    y_max1<-max(gamma((const_par[2]+m_Mode)/2)*((m_Mode/2)^(m_Mode/2))*(Mode_seq1^(m_Mode/2-1))/gamma(m_Mode/2)/gamma(const_par[2]/2)/((const_par[2]/2)^(m_Mode/2))/((1+m_Mode*Mode_seq1/const_par[2])^((m_Mode+const_par[2])/2)))
    Mode_seq2<-n*(const_par[1]-2)/const_par[1]/(n+2)
    y_max2<-max(gamma((const_par[1]+n)/2)*((const_par[1]/2)^(const_par[1]/2))*(Mode_seq2^(const_par[1]/2-1))/gamma(const_par[1]/2)/gamma(n/2)/((n/2)^(const_par[1]/2))/((1+const_par[1]*Mode_seq2/n)^((const_par[1]+n)/2)))
    return(c(0.0001,4*max(Mode_seq1),y_max1,0.0001,4*max(Mode_seq2),y_max2))
  }
  
  if(is.numeric(x))
  {
    
    density<-ifelse(x>=0,gamma((m+n)/2)/gamma(m/2)/gamma(n/2)*((m/2)^(m/2))*(x^(m/2-1))/((n/2)^(m/2))/((1+m/n*x)^(m/2+n/2)),0)
    fun<-function(x){gamma((m+n)/2)/gamma(m/2)/gamma(n/2)*((m/2)^(m/2))*(x^(m/2-1))/((n/2)^(m/2))/((1+m/n*x)^(m/2+n/2))}
    cdf<-c()
    for(i in x)
    {
      cdf<-c(cdf,as.numeric(ifelse(i>=0.0001,integrate(fun,0.0001,i),0)))
    }
    
    Mean<-ifelse(n>2,n/(n-2),Inf)
    Variance<-ifelse(n>4,2*n^2*(n+m-2)/m/(n-2)^2/(n-4),Inf)
    Mode<-ifelse(m>2,n*(m-2)/(n+2)/m,Inf)
    Skewness<-ifelse(n>6,(2*m+n-2)*(8*(n-4))^0.5/(n-6)/(m*(m+n-2))^0.5,Inf) 
    Kurtosis<-ifelse(n>8,3+12*((n-2)*2*(n-4)+m*(m+n-2)*(5*n-22))/(m*(n-6)*(n-8)*(m+n-2)),Inf)
    return(list(density=density,cdf=cdf,Mean=Mean,Variance=Variance,Mode=Mode,Skewness=Skewness,Kurtosis=Kurtosis,m=m,'m',n=n,'n',"F Distribution"))
  }
}


Gamma_dis<-function(x=NULL,para,const_par)
{
  a<-para[1,];b<-para[2,]
  const_par=const_par
  if(any(a<=0)||any(b<=0))
    stop("Error! The parameter must be positive integer!")
  if(const_par[1]<=0||const_par[2]<=0)
    stop("Error! The constant parameter must be positive integer!")
  
  if(!is.numeric(x))
  {
    if(all(a<=1))
      stop("Error! In order to calculate the mode, the largest number of a must larger than 2!")
    if(const_par[1]<=1)
      stop("Error! In order to calculate the mode, constant a must larger than 2!")
    
    if(any(a<=1))
    {
      a_Mode=a[-which(a<=1)]
      b_Mode=b[-which(a<=1)]
    }
    else
    {
      a_Mode=a
      b_Mode=b
    }
    Mode_seq1<-const_par[2]*(a_Mode-1)
    sq1<-(a_Mode*(const_par[2]^2))^0.5
    y_max1<-max((exp(-Mode_seq1/const_par[2]))*(Mode_seq1^(a_Mode-1))/(gamma(a_Mode)*(const_par[2]^a_Mode)))
    Mode_seq2<-b_Mode*(const_par[1]-1)
    sq2<-(const_par[1]*(b_Mode^2))^0.5
    temp1<-unique(which.max(sq1))
    temp2<-unique(which.max(sq2))
    y_max2<-max((exp(-Mode_seq2/b_Mode))*(Mode_seq2^(const_par[1]-1))/(gamma(const_par[1])*(b_Mode^const_par[1])))
    return(c(max(Mode_seq1[temp1]-3*sq1[temp1],0.0001),Mode_seq1[temp1]+3*sq1[temp1],y_max1,max(Mode_seq2[temp2]-3*sq2[temp2],0.0001),Mode_seq2[temp2]+3*sq2[temp2],y_max2))
  }
  
  if(is.numeric(x))
  {
    
    density<-ifelse(x>=0,(exp(-x/b))*(x^(a-1))/(gamma(a)*(b^a)),0)
    cdf<-c()
    for(i in x)
    {
      cdf<-c(cdf,pgamma(i,shape=a,scale=b))
    }
    
    Mean<-a*b
    Variance<-a*(b^2)
    Mode<-ifelse(a>1,b*(a-1),Inf)
    Skewness<-2/(a^0.5)
    Kurtosis<-3+(6/a)
    return(list(density=density,cdf=cdf,Mean=Mean,Variance=Variance,Mode=Mode,Skewness=Skewness,Kurtosis=Kurtosis,a=a,'a=',b=b,'b=',"Gamma Distribution"))
  }
}

Inverse_Gaussian<-function(x=NULL,para,const_par)
{
  lamda<-para[1,];mu<-para[2,]
  const_par=const_par
  if(any(lamda<=0)||any(mu<=0))
    stop("Error! The parameter lamda and mu must be positive numbers!")
  if(const_par[1]<=0||const_par[2]<=0)
    stop("Error! The constant parameter lamda and mu must be positive numbers!")
  
  if(!is.numeric(x))
  {
    lamda=subset(lamda,lamda>0)
    mu=subset(mu,mu>0)
    Mode_seq1<-const_par[2]*(((1+(9*(const_par[2]^2)/(4*(lamda^2))))^0.5)-(3*const_par[2])/(2*lamda))
    y_max1<-max(((lamda/(2*pi*(Mode_seq1^3)))^0.5)*(exp((-lamda*((Mode_seq1-const_par[2])^2))/(2*(const_par[2]^2)*Mode_seq1))))
    Mode_seq2<-mu*(((1+(9*(mu^2)/(4*(const_par[1]^2))))^0.5)-(3*mu)/(2*const_par[1]))
    y_max2<-max(((const_par[1]/(2*pi*(Mode_seq2^3)))^0.5)*(exp((-const_par[1]*((Mode_seq2-mu)^2))/(2*(mu^2)*Mode_seq2))))
    return(c(0.0001,10,y_max1,0.0001,10,y_max2))
  }
  
  if(is.numeric(x))
  {
    
    density<-((lamda/(2*pi*(x^3)))^0.5)*(exp((-lamda*((x-mu)^2))/(2*(mu^2)*x)))
    fun<-function(x){((lamda/(2*pi*(x^3)))^0.5)*(exp((-lamda*((x-mu)^2))/(2*(mu^2)*x)))}
    cdf<-c()
    for(i in x)
    {
      cdf<-c(cdf,as.numeric(integrate(fun,0.0001,i)$value))
    }
    
    Mean<-mu
    Variance<-(mu^3)/lamda
    Mode<-mu*(((1+9*(mu^2)/(4*(lamda^2)))^0.5)-(3*mu/(2*lamda)))
    Skewness<-3*((mu/lamda)^0.5)
    Kurtosis<-3+15*mu/lamda
    return(list(density=density,cdf=cdf,Mean=Mean,Variance=Variance,Mode=Mode,Skewness=Skewness,Kurtosis=Kurtosis,lamda=lamda,'lamda',mu=mu,'mu',"Inverse Gaussian Distribution"))
  }
}

Laplace<-function(x=NULL,para,const_par)
{
  a<-para[1,];b<-para[2,]
  const_par=const_par
  if(any(b<=0))
    stop("Error! The parameter must be positive integer!")
  if(const_par[2]<=0)
    stop("Error! The constant parameter must be positive integer!")
  
  if(!is.numeric(x))
  {
    a=a
    b=subset(b,b>0)
    Sd_seq1<-(2*(const_par[2]^2))^0.5
    y_max1<-max(1/(2*const_par[2]))
    Sd_seq2<-(2*(b^2))^0.5
    y_max2<-max(1/(2*b))
    return(c(-4*Sd_seq1+min(a),4*Sd_seq1+max(a),y_max1,-4*max(Sd_seq2)+const_par[1],4*max(Sd_seq2)+const_par[1],y_max2))
  }
  
  if(is.numeric(x))
  {
    
    density<-exp(-(abs(x-a))/b)/(2*b)
    fun<-function(x){exp(-(abs(x-a))/b)/(2*b)}
    cdf<-c()
    for(i in x)
    {
      cdf<-c(cdf,as.numeric(integrate(fun,-Inf,i)$value))
    }
    
    Mean<-a
    Variance<-2*(b^2)
    Mode<-a
    Skewness<-0
    Kurtosis<-6
    return(list(density=density,cdf=cdf,Mean=Mean,Variance=Variance,Mode=Mode,Skewness=Skewness,Kurtosis=Kurtosis,a=a,'a',b=b,'b',"Laplace Distribution"))
  }
}

Logistic<-function(x=NULL,para,const_par)
{
  a<-para[1,];b<-para[2,]
  # const_par=const_par
  if(any(b<=0))
    stop("Error! The parameter must be positive integer!")
  if(const_par[2]<=0)
    stop("Error! The constant parameter must be positive integer!")
  
  if(!is.numeric(x))
  {
    #  a=a
    #   b=subset(b,b>0)
    Sd_seq1<-((const_par[2]^2)*(pi^2)/3)^0.5
    y_max1<-max(1/(4*const_par))
    Sd_seq2<-((b^2)*(pi^2)/3)^0.5
    y_max2<-max(1/(4*b))
    return(c(-5*Sd_seq1+min(a),5*Sd_seq1+max(a),y_max1,-5*max(Sd_seq2)+const_par[1],5*max(Sd_seq2)+const_par[1],y_max2))
  }
  
  if(is.numeric(x))
  {
    
    density<-(exp(-(x-a)/b))/(b*((1+exp(-(x-a)/b))^2))
    #  fun<-function(x){(exp(-(x-a)/b))/(b*((1+exp(-(x-a)/b))^2))}
    cdf<-c()
    for(i in x)
    {
      cdf<-c(cdf,(1+exp(-(i-a)/b))^(-1))
    }
    
    Mean<-a
    Variance<-(b^2)*(pi^2)/3
    Mode<-a
    Skewness<-0
    Kurtosis<-4.2
    return(list(density=density,cdf=cdf,Mean=Mean,Variance=Variance,Mode=Mode,Skewness=Skewness,Kurtosis=Kurtosis,a=a,'a',b=b,'b',"Logistic Distribution"))
  }
}

Lognormal<-function(x=NULL,para,const_par)
{
  mu<-para[1,];sigma<-para[2,]
  const_par=const_par
  if(any(sigma<=0))
    stop("Error! The parameter must be positive integer!")
  if(const_par[2]<=0)
    stop("Error! The constant parameter must be positive integer!")
  
  if(!is.numeric(x))
  {
    mu=mu
    sigma=subset(sigma,sigma>0)
    Mode_seq1<-exp(mu-(const_par[2]^2))
    y_max1<-max(exp((-((log(Mode_seq1)-mu)^2))/(2*(const_par[2]^2)))/(((2*pi)^0.5)*Mode_seq1*pi))
    Mode_seq2<-exp(const_par[1]-(sigma^2))
    y_max2<-max(exp((-((log(Mode_seq2)-const_par[1])^2))/(2*(sigma^2)))/(((2*pi)^0.5)*Mode_seq2*pi))
    return(c(0.0001,3*max(Mode_seq1),y_max1,0.0001,3*max(Mode_seq2),y_max2))
  }
  
  if(is.numeric(x))
  {
    
    density<-exp((-((log(x)-mu)^2))/(2*(sigma^2)))/(((2*pi)^0.5)*x*pi)
    fun<-function(x){exp((-((log(x)-mu)^2))/(2*(sigma^2)))/(((2*pi)^0.5)*x*pi)}
    cdf<-c()
    for(i in x)
    {
      cdf<-c(cdf,as.numeric(integrate(fun,0.0001,i)$value))
    }
    
    Mean<-exp(mu+(sigma^2)/2)
    Variance<-(exp(sigma^2))*((exp(sigma^2))-1)*(exp(2*mu))
    Mode<-exp(mu-(sigma^2))
    Skewness<-(exp(sigma^2)+2)*(((exp(sigma^2))-1)^0.5)
    Kurtosis<-exp(4*(sigma^2))+2*exp(3*(sigma^2))+3*exp(2*(sigma^2))-3
    return(list(density=density,cdf=cdf,Mean=Mean,Variance=Variance,Mode=Mode,Skewness=Skewness,Kurtosis=Kurtosis,mu=mu,'mu',sigma=sigma,'sigma',"Lognormal Distribution"))
  }
}

Normal<-function(x=NULL,para,const_par)
{
  mu<-para[1,];sigma<-para[2,]
  if(any(sigma<=0))
    stop("Error! The parameter sigma must be positive!")
  if(!is.numeric(x))
  {
    # mu=mu
    #  sigma=subset(sigma,sigma>0)
    y_max1<-max(1/(const_par[2]*((2*pi)^0.5)))
    y_max2<-max(1/(sigma*((2*pi)^0.5)))
    return(c(-5*const_par[2]+min(mu),5*const_par[2]+max(mu),y_max1,-5*max(sigma)+const_par[1],5*max(sigma)+const_par[1],y_max2))
  }
  if(is.numeric(x))
  {
    density<-(exp(-((x-mu)^2)/(2*(sigma^2))))/(sigma*((2*pi)^0.5))
    # fun<-function(x){(exp(-((x-mu)^2)/(2*(sigma^2))))/(sigma*((2*pi)^0.5))}
    cdf<-c()
    for(i in x)
    {
      #cdf<-c(cdf,integrate(fun,-Inf,i)$value)
      cdf<-c(cdf,pnorm(i,mu,sigma))
    }
    
    Mean<-mu
    Variance<-sigma^2
    Mode<-mu
    Skewness<-0 
    Kurtosis<-3
    return(list(density=density,cdf=cdf,Mean=Mean,Variance=Variance,Mode=Mode,Skewness=Skewness,Kurtosis=Kurtosis,mu=mu,'mu',sigma=sigma,'sigma',"Normal Distribution"))
  }
}

Pareto<-function(x=NULL,para,const_par)
{
  a<-para[1,];b<-para[2,]
  const_par=const_par
  if(any(a<=0)||any(b<=0))
    stop("Error! The parameter must be positive integer!")
  if(const_par[1]<=0||const_par[2]<=0)
    stop("Error! The constant parameter must be positive integer!")
  
  if(!is.numeric(x))
  {
    a=subset(a,a>0)
    b=subset(b,b>0)
    Mode_seq1<-a
    y_max1<-max(const_par[2]*(a^const_par[2])/(Mode_seq1^(const_par[2]+1)))
    Mode_seq2<-const_par[1]
    y_max2<-max(b*(const_par[1]^b)/(Mode_seq2^(b+1)))
    return(c(min(a),6,y_max1,min(const_par[1]),6,y_max2))
  }
  
  if(is.numeric(x))
  {
    
    density<-ifelse(x<a,0,b*(a^b)/(x^(b+1)))
    fun<-function(x){b*(a^b)/(x^(b+1))}
    cdf<-c()
    for(i in x)
    {
      cdf<-c(cdf,as.numeric(integrate(fun,a,i)$value))
    }
    
    Mean<-ifelse(b>1,a*b/(b-1),Inf)
    Variance<-ifelse(b>2,b*(a^2)/(((b-1)^2)*(b-2)),Inf)
    Mode<-a
    Skewness<-ifelse(b>3,2*(b+1)*(((b-2)/b)^0.5)/(b-3),Inf)
    Kurtosis<-ifelse(b>4,3*(b-2)*(3*(b^2)+b+2)/(b*(b-3)*(b-4)),Inf)
    return(list(density=density,cdf=cdf,Mean=Mean,Variance=Variance,Mode=Mode,Skewness=Skewness,Kurtosis=Kurtosis,a=a,'a',b=b,'b',"Pareto Distribution"))
  }
}
Bernoulli<-function(k=NULL,para)
{
  p<-para
  if(any(p<=0||p>=1))
    stop("Error! The parameter must be positive and less than 1!")
  
  if(!is.numeric(k))
  {
    y_max<-max(p,1-p)
    return(c(0,1,y_max))
  }
  
  if(is.numeric(k))
  {
    density<-ifelse(k==1,p,1-p)
    fun<-function(k){ifelse(k==1,p,1-p)}
    cdf<-c()
    for(i in k)
    {
      cdf<-c(cdf,sum(fun(0:i)))
    }
    Mean<-p
    Variance<-p*(1-p)
    Mode<-ifelse(p<0.5,0,ifelse(p==0.5,0.5,1))
    Skewness<-(1-2*p)/((p*(1-p))^0.5)
    Kurtosis<-(1-6*p*(1-p))/(p*(1-p))
    
    return(list(density=density,
                cdf=cdf,
                Mean=Mean,
                Variance=Variance,
                Mode=Mode,
                Skewness=Skewness,
                Kurtosis=Kurtosis,
                p=p,
                'p',
                "Bernoulli Distribution"))
  }
}

Dis_Uniform<-function(k=NULL,para)
{
  N<-para
  N<-seq(round(min(N)),round(max(N)),1)
  if(any(N<=0))
    stop("Error! The parameter must be positive integer!")
  
  if(!is.numeric(k))
  {
    y_max<-max(1/N)
    return(c(1,max(N),y_max))
  }
  
  if(is.numeric(k))
  {
    density<-ifelse(k>N,0,(1/N))
    fun<-function(k){ifelse(k>N,0,(1/N))}
    cdf<-c()
    for(i in k)
    {
      cdf<-c(cdf,sum(fun(1:i)))
    }
    Mean<-(N+1)/2
    Variance<-(N-1)*(N+1)/12
    Mode<-Inf
    Skewness<-0
    Kurtosis<-3-6*(N^2+1)/(5*(N-1)*(N+1))
    
    return(list(density=density,
                cdf=cdf,
                Mean=Mean,
                Variance=Variance,
                Mode=Mode,
                Skewness=Skewness,
                Kurtosis=Kurtosis,
                N=N,
                'N',
                "Discrete Uniform Distribution"))
  }
}

Geometric<-function(k=NULL,para)
{
  p<-para
  if(any(p<=0||p>=1))
    stop("Error! The parameter must be positive and less than 1!")
  
  if(!is.numeric(k))
  {
    Mode_seq<-0
    y_max<-max(((1-p)^Mode_seq)*p)
    return(c(0,20,y_max))
  }
  
  if(is.numeric(k))
  {
    density<-((1-p)^k)*p
    fun<-function(k){((1-p)^k)*p}
    cdf<-c()
    for(i in k)
      cdf<-c(cdf,sum(fun(0:i)))
    Mean<-(1-p)/p
    Variance<-(1-p)/(p^2)
    Mode<-0
    Skewness<-(2-p)/((1-p)^0.5)
    Kurtosis<-9+(p^2)/(1-p)
    
    return(list(density=density,
                cdf=cdf,
                Mean=Mean,
                Variance=Variance,
                Mode=Mode,
                Skewness=Skewness,
                Kurtosis=Kurtosis,
                p=p,
                'p',
                "Geometric Distribution"))
  }
}

Logarithmic_Series<-function(k=NULL,para)
{
  theta<-para
  if(any(theta<=0||theta>=1))
    stop("Error! The parameter must be positive and less than 1!")
  
  if(!is.numeric(k))
  {
    return(c(1,20,1))
  }
  
  if(is.numeric(k))
  {
    density<--(theta^k)/((log(1-theta))*k)
    fun<-function(k){-(theta^k)/((log(1-theta))*k)}
    cdf<-c()
    for(i in k)
      cdf<-c(cdf,sum(fun(1:i)))
    a<--1/(log(1-theta))
    Mean<-a*theta/(1-theta)
    Variance<-a*theta*(1-a*theta)/((1-theta)^2)
    Mode<-Inf
    Skewness<-a*theta*(1+theta-3*a*theta+2*(a^2)*(theta^2))/((a*theta*(1-a*theta))^(1.5))
    Kurtosis<-(1+4*theta+(theta^2)-4*a*theta*(1+theta)+6*(a^2)*(theta^2)-3*(a^3)*(theta^3))/(a*theta*((1-a*theta)^2))
    
    return(list(density=density,
                cdf=cdf,
                Mean=Mean,
                Variance=Variance,
                Mode=Mode,
                Skewness=Skewness,
                Kurtosis=Kurtosis,
                theta=theta,
                'theta',
                "Logarithmic Series Distribution"))
  }
}

Poisson<-function(k=NULL,para)
{
  lamda<-para
  if(any(lamda<=0))
    stop("Error! The parameter must be positive!")
  
  if(!is.numeric(k))
  {
    Mode_seq<-floor(lamda)
    y_max<-max((exp(-lamda))*(lamda^Mode_seq)/(factorial(Mode_seq)))
    return(c(0,2*max(Mode_seq),y_max))
  }
  
  if(is.numeric(k))
  {
    density<-(exp(-lamda))*(lamda^k)/(factorial(k))
    fun<-function(k){(exp(-lamda))*(lamda^k)/(factorial(k))}
    cdf<-c()
    for(i in k)
      cdf<-c(cdf,sum(fun(0:i)))
    Mean<-lamda
    Variance<-lamda
    Mode<-floor(lamda)
    Skewness<-1/(lamda^0.5)
    Kurtosis<-1/(lamda^0.5)
    
    return(list(density=density,
                cdf=cdf,
                Mean=Mean,
                Variance=Variance,
                Mode=Mode,
                Skewness=Skewness,
                Kurtosis=Kurtosis,
                lamda=lamda,
                'lamda',
                "Poisson Distribution"))
  }
}

Binomial<-function(k=NULL,para,const_par)
{
  N<-para[1,];p<-para[2,]
  # const_par=const_par
  if(any(N<=0)||any(p<=0||p>=1))
    stop("Error! The parameter N must be positive integer and p must be positive and less than 1!")
  if(const_par[1]<=0||(const_par[2]<=0||const_par[2]>=1))
    stop("Error! The constant parameter N must be positive integer and p must be positive and less than 1!")
  
  if(!is.numeric(k))
  {
    Mode_seq1<-floor((N+1)*const_par[2])
    Mode_seq2<-floor((const_par[1]+1)*p)
    y_max1<-max((choose(N,Mode_seq1))*(const_par[2]^Mode_seq1)*((1-const_par[2])^(N-Mode_seq1)))
    y_max2<-max((choose(const_par[1],Mode_seq2))*(p^Mode_seq2)*((1-p)^(const_par[1]-Mode_seq2)))
    return(c(0,max(N),y_max1,0,const_par[1],y_max2))
  }
  
  if(is.numeric(k))
  {
    density<-ifelse(k>N,0,(choose(N,k))*(p^k)*((1-p)^(N-k)))
    fun<-function(k){ifelse(k>N,0,(choose(N,k))*(p^k)*((1-p)^(N-k)))}
    cdf<-c()
    for(i in k)
    {
      cdf<-c(cdf,sum(fun(0:i)))
    }
    
    Mean<-N*p
    Variance<-N*p*(1-p)
    Mode<-floor((N+1)*p)
    Skewness<-(1-2*p)/((N*p*(1-p))^0.5)
    Kurtosis<-3-6/N+1/(N*p*(1-p))
    return(list(density=density,cdf=cdf,Mean=Mean,Variance=Variance,Mode=Mode,Skewness=Skewness,Kurtosis=Kurtosis,N=N,'N',p=p,'p',"Binomial Distribution"))
  }
}

Negative_Binomial<-function(k=NULL,para,const_par)
{
  r<-para[1,];p<-para[2,]
  const_par=const_par
  if(any(r<=0)||any(p<=0||p>=1))
    stop("Error! The parameter r must be positive integer and p must be positive and less than 1!")
  if(const_par[1]<=0||(const_par[2]<=0||const_par[2]>=1))
    stop("Error! The constant parameter r must be positive integer and p must be positive and less than 1!")
  
  if(!is.numeric(k))
  {
    Mode_seq1<-floor((r-1)*(1-const_par[2])/const_par[2])
    Mode_seq2<-floor((const_par[1]-1)*(1-p)/p)
    y_max1<-max((choose(r+Mode_seq1-1,Mode_seq1))*(const_par[2]^r)*((1-const_par[2])^Mode_seq1))
    y_max2<-max((choose(const_par[1]+Mode_seq2-1,Mode_seq2))*(p^const_par[1])*((1-p)^Mode_seq2))
    return(c(0,2*max(Mode_seq1),y_max1,0,2*max(Mode_seq2),y_max2))
  }
  
  if(is.numeric(k))
  {
    density<-(choose(r+k-1,k))*(p^r)*((1-p)^k)
    fun<-function(k){(choose(r+k-1,k))*(p^r)*((1-p)^k)}
    cdf<-c()
    for(i in k)
    {
      cdf<-c(cdf,sum(fun(0:i)))
    }
    
    Mean<-r*(1-p)/p
    Variance<-r*(1-p)/(p^2)
    Mode<-floor((r-1)*(1-p)/p)
    Skewness<-(2-p)/((r*(1-p))^0.5)
    Kurtosis<-3+6/r+(p^2)/(r*(1-p))
    return(list(density=density,cdf=cdf,Mean=Mean,Variance=Variance,Mode=Mode,Skewness=Skewness,Kurtosis=Kurtosis,r=r,'r',p=p,'p',"Negative Binomial Distribution"))
  }
}

Hypergeometric<-function(k=NULL,para,const_par)
{
  n<-para[1,];M<-para[2,];N<-para[3,]
  #  const_par=const_par
  if(any(n<=0)||any(M<=0)||any(N<=0))
    stop("Error! The parameter n,M,N must be positive integer!")
  if(const_par[1]<=0||const_par[2]<=0||const_par[3]<=0)
    stop("Error! The constant parameter n,M,N must be positive integer!")
  if(any(n>const_par[3])||any(N<const_par[1]))
    stop("Error! The parameter n must be not larger than N!")
  
  if(!is.numeric(k))
  {
    Mode_seq1<-floor(((n+1)*(const_par[2]+1))/(const_par[3]+2))
    Mode_seq2<-floor(((const_par[1]+1)*(M+1))/(const_par[3]+2))
    Mode_seq3<-floor(((const_par[1]+1)*(const_par[2]+1))/(N+2))
    y_max1<-max(choose(const_par[2],Mode_seq1)*choose(const_par[3]-const_par[2],n-Mode_seq1)/choose(const_par[3],n))
    y_max2<-max(choose(M,Mode_seq2)*choose(const_par[3]-M,const_par[1]-Mode_seq2)/choose(const_par[3],const_par[1]))
    y_max3<-max(choose(const_par[2],Mode_seq3)*choose(N-const_par[2],const_par[1]-Mode_seq3)/choose(N,const_par[1]))
    return(c(min(max(0,const_par[2]-const_par[3]+n)),max(min(n,const_par[2])),y_max1,min(max(0,M-const_par[3]+const_par[1])),max(min(const_par[1],M)),y_max2,min(max(0,const_par[2]-N+n)),max(min(const_par[1],const_par[2])),y_max3))
  }
  
  if(is.numeric(k))
  {
    L<-max(0,M-N+n)
    U<-min(n,M)
    density<-ifelse(k<=U&k>=L,choose(M,k)*choose(N-M,n-k)/choose(N,n),0)
    fun<-function(k){ifelse(k<=U&k>=L,choose(M,k)*choose(N-M,n-k)/choose(N,n),0)}
    cdf<-c()
    for(i in k)
    {
      cdf<-c(cdf,sum(fun(0:i)))
    }
    
    Mean<-n*M/N
    Variance<-n*M*(1-M/N)*((N-n)/(N-1))/N
    Mode<-floor(((n+1)*(M+1))/(N+2))
    Skewness<-((N-2*M)*(N-2*n)*((N-1)^0.5))/((N-2)*((n*M*(N-M))^0.5)*(N-n))
    Kurtosis<-(((N^2)*(N-1))/(n*M*(N-M)*(N-2)*(N-3)*(N-n)))*((3*n*M*(N-M)*(6-n)/N)+(N*(N+1-6*n))+(6*(n^2))+(3*M*(N-M)*(n-2))-(18*(n^2)*M*(N-M)/(N^2)))
    return(list(density=density,cdf=cdf,Mean=Mean,Variance=Variance,Mode=Mode,Skewness=Skewness,Kurtosis=Kurtosis,n=n,'n',M=M,'M',N=N,'N',"Hypergeometric Distribution"))
  }
}
