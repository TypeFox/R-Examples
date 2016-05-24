ipd.sep <- function(
               effect,
               event0=NULL,
               event1=NULL,
               mean0=NULL,
               mean1=NULL,
               var0=NULL,
               var1=NULL,
               x0=NULL,
               x1=NULL,
               s20=NULL,
               s21=NULL,
               n0=NULL,
               n1=NULL,
               data,
               alpha=.05)
{

#REPEATED CALLS TO IPD.META.POWER.BASE.FUNCTION

    types <- c("point","lower","upper")

    if(!missing(data)){ #CHECK IF OBJECTS ARE IN GIVEN DATA FRAME
      
      CALL <- match.call()
      arg.names <- lapply(CALL,deparse)
      
      if(!is.null(arg.names$event0)&&arg.names$event0%in%names(data)){
        event0 <- data[,arg.names$event0]
      }

      if(!is.null(arg.names$event1)&&arg.names$event1%in%names(data)){
        event1 <- data[,arg.names$event1]
      }

      if(!is.null(arg.names$mean0)&&arg.names$mean0%in%names(data)){
        mean0 <- data[,arg.names$mean0]
      }

      if(!is.null(arg.names$mean1)&&arg.names$mean1%in%names(data)){
        mean1 <- data[,arg.names$mean1]
      }
      
      if(!is.null(arg.names$var0)&&arg.names$var0%in%names(data)){
        var0 <- data[,arg.names$var0]
      }

      if(!is.null(arg.names$var1)&&arg.names$var1%in%names(data)){
        var1 <- data[,arg.names$var1]
      }

      if(!is.null(arg.names$n0)&&arg.names$n0%in%names(data)){
        n0 <- data[,arg.names$n0]
      }

      if(!is.null(arg.names$n1)&&arg.names$n1%in%names(data)){
        n1 <- data[,arg.names$n1]
      }
      
      if(!is.null(arg.names$x0)&&arg.names$x0%in%names(data)){
        x0 <- data[,arg.names$x0]
      }

      if(!is.null(arg.names$x1)&&arg.names$x1%in%names(data)){
        x1 <- data[,arg.names$x1]
      }

      if(!is.null(arg.names$s20)&&arg.names$s20%in%names(data)){
        s20 <- data[,arg.names$s20]
      }

      if(!is.null(arg.names$s21)&&arg.names$s21%in%names(data)){
        s21 <- data[,arg.names$s21]
      }
      
    }

    if(is.null(s20)){
      s20 <- x0*(1-x0)
      s21 <- x1*(1-x1)
    }
    
    estimates <- lapply(types,function(x){

       ipd.meta.power.base.function(type=x,
	event0=event0,event1=event1,mean0=mean0,mean1=mean1,
        var0=var0,var1=var1,x0=x0,x1=x1,s20=s20,s21=s21,
	n0=n0,n1=n1,effect=effect,alpha=alpha)
    })

    SE <- estimates[[1]]$estimated.se
    SE.l <- estimates[[2]]$estimated.se
    SE.u <- estimates[[3]]$estimated.se

    power <- estimates[[1]]$estimated.power
    power.lower <- estimates[[3]]$estimated.power
    power.upper <- estimates[[2]]$estimated.power

return(
      list(
         estimated.power=power,
         power.lower=power.lower,
         power.upper=power.upper,
         estimated.se=SE,
         se.lower=SE.l,
         se.upper=SE.u,
         sigma=estimates[[1]]$sigma,
         sigma0=estimates[[1]]$sigma0,
         sigma1=estimates[[1]]$sigma1,
         level=alpha)
         )
}


ipd.meta.power.base.function <- function(
                             effect,
                             event0=NULL,
                             event1=NULL,
                             mean0=NULL,
                             mean1=NULL,
                             var0=NULL,
                             var1=NULL,
                             x0,
                             x1,
                             s20=NULL,
                             s21=NULL,
                             n0,
                             n1,
                             alpha=.05,
                             type=c("point","lower","upper"))
{

######DEPENDENT FUNCTIONS

power <- function(beta,se,alpha=.05){
      z <- qnorm(1-alpha/2)
      return(pnorm(-z+beta/se)+pnorm(-z-beta/se))
}

DSL <- function(y,var){

  w <- 1/var
  df <- length(y)-1
  y.bar <- sum(w*y)/sum(w)
  Q <- sum(w*(y-y.bar)^2)

  max(c(0,(Q-df)/(sum(w)-sum(w^2)/sum(w))))
  
}

get.weight <- function(y0,y1,n0,n1,BINARY=TRUE,z){

 if(BINARY){

 theta.bar <- sum(exp(y0)/(1+exp(y0))*n0+exp(y1)/(1+exp(y1))*n1)/sum(n0+n1)
 w.se <- (1-2*theta.bar)*sqrt(theta.bar*(1-theta.bar)/sum(n0+n1))
 weight <- theta.bar*(1-theta.bar)+z*w.se
 return(1/weight)

  }
 else{

 return(1)
 }
}

construct.beta.info.matrix <- function(n0,n1,x0,x1,s20,s21,d,e,f,g){

#RETURNS THE X'SIGMA^(-1)X INFORMATION FOR BETA OF IPD STUDY
#GIVEN COMPONENTS OF WEIGHT MATRIX SIGMA^(-1) COVARIATE SUMMARIES
#d, e, f, g are vectors that are the length of the number of studies K
#n0 control study sample sizes
#n1 treatment study sample sizes
#x0 control study means for covariate
#x1 treatment study means for covariate
#s20 control sample variances for covariate
#s21 treatment sample variances

c.factors <- function(n0,n1,x0,x1,e,f,g){

#control = 0
#treatment = 1
#x are study means

c.ctrl <- e*n0*x0+f*n1*x1
c.trt <- f*n0*x0+g*n1*x1

return(list(c0=c.ctrl,c1=c.trt))
}


######C FACTORS FOR EACH COLUMN VECTOR

c0 <- c.factors(n0,n1,rep(1,length(n0)),rep(1,length(n1)),e,f,g)
c1 <- c.factors(n0,n1,rep(0,length(n0)),rep(1,length(n1)),e,f,g)
c2 <- c.factors(n0,n1,x0,x1,e,f,g)
c3 <- c.factors(n0,n1,rep(0,length(n0)),x1,e,f,g)

#FIRST ROW

X00 <- sum(n0*(d+c0$c0)+n1*(d+c0$c1))
X01 <- sum(n1*(d+c0$c1))
X02 <- sum(n0*x0*(d+c0$c0)+n1*x1*(d+c0$c1))
X03 <- sum(n1*x1*(d+c0$c1))

#SECOND ROW

X11 <- sum(n1*(d+c1$c1))
X12 <- sum(n1*x1*(d+c1$c1)+n0*c1$c0*x0)
X13 <- sum(n1*x1*(d+c1$c1))

#THIRD ROW

X22 <-
sum(d*(n0*(s20+x0^2)+n1*(s21+x1^2)-(s20+s21))+n0*c2$c0*x0+n1*c2$c1*x1)
X23 <- sum(d*(n1*(s21+x1^2)-(s21))+n1*c2$c1*x1)

#LAST COMPONENT
X33 <- sum(d*(n1*(s21+x1^2)-(s21))+n1*c3$c1*x1)

Sigma <- 

matrix(
c(X00,X01,X02,X03,
X01,X11,X12,X13,
X02,X12,X22,X23,
X03,X13,X23,X33),4,4
)
return(Sigma)
}

me.inverse <- function(a,b,c,n,m){

#Determine the study inverse of the marginal variance for a mixed effects model
#mixed effects for intercept and treatment group

#a+b is diagonal for control
#b is control diagonal (covariance)
#c is the off diagonal elements for treatment
#m treated subjects; n control subjects

block.matrix.inverse <- function(n,a,b){

#Returns the 

a.inverse <- 1/a
b.inverse <- -b/(a*(a+n*b))

return(c(a.inverse,off.diag=b.inverse,diag=a.inverse+b.inverse))
}

####

#E matrix (bottom left)

A.inverse.components <- block.matrix.inverse(n,a,b)

a.breve <- A.inverse.components[1]
b.breve <- A.inverse.components[2]

gamma <- b^2*n*(a.breve+n*b.breve)

E.inverse <- block.matrix.inverse(m,a,c-gamma)

#FE.inverse (off diagonals)

a.tilde <- E.inverse[1]
c.tilde <- E.inverse[2]

FE.inverse <- -b*(a.breve+n*b.breve)*(a.tilde+m*c.tilde)

#A.inverse (upper left)

gamma.breve <- b*m*(a.breve+n*b.breve)*(-FE.inverse)
diag=a.breve+b.breve+gamma.breve
off=b.breve+gamma.breve

d=a.breve
e=b.breve+gamma.breve
f=FE.inverse
g=E.inverse[2]

return(
c(d,e,f,g)
)
}

#########PROCESSING

 BINARY <- !is.null(event0)

 if(!BINARY){
  #SAMPLE VARIANCE PER STUDY; RESIDUAL ESTIMATE FOR LINEAR MODEL
 sigma <- (var0*(n0-1)+var1*(n1-1))/(n0+n1-2)
 y0 <- mean0
 y1 <- mean1
  }
else{
  #COMPUTE LOG-ODDS WITH CORRECTION FOR ZERO-EVENTS
  event0 <- ifelse(event0==0,.5,event0)
  event1 <- ifelse(event1==0,.5,event1)
  event0 <- ifelse(event0==n0,event0-.5,event0)
  event1 <- ifelse(event1==n1,event1-.5,event1)
  
  y0 <- log(event0/n0/(1-event0/n0))
  y1 <- log(event1/n1/(1-event1/n1))
  var0 <- 1/(event0/n0*(1-event0/n0))
  var1 <- 1/(event1/n1*(1-event1/n1))
  
  sigma <- 1
 }

tau20 <- DSL(y0,var0/n0)
tau21 <- DSL(y1,var1/n1)

#WALD-TYPE CONFIDENCE INTERVAL
z <- qnorm((1-alpha/2))

lower0 <- max(0,tau20-z*sqrt(2*1/(sum(n0/var0))))
upper0 <- max(0,tau20+z*sqrt(2*1/(sum(n0/var0))))
           
 #OBTAINING BETWEEN-STUDY VARIANCE COMPONENTS

if(type=="point"){
 sigma0 <- tau20
 sigma1 <- tau21-sigma0
 sigma1 <- max(c(0,sigma1))
 weight <- get.weight(y0,y1,n0,n1,BINARY=BINARY,z=0)
 }
else if(type=="lower"){

sigma0 <- lower0
sigma1 <- tau21-sigma0
se.sigma1 <- sqrt(2/(sum(n1/var1))+2/(sum(n1/var1)))
sigma1 <- max(c(0,sigma1-z*se.sigma1))
weight <- get.weight(y0,y1,n0,n1,BINARY=BINARY,z=z)

 }
else{

sigma0 <- upper0
sigma1 <- tau21-upper0
se.sigma1 <- sqrt(2/(sum(n1/var1))+2/(sum(n1/var1)))
sigma1 <- max(c(0,sigma1+z*se.sigma1))
weight <- get.weight(y0,y1,n0,n1,BINARY=BINARY,z=-z)

}

inverse.components <- mapply(me.inverse,n=n0,m=n1,a=sigma,
                             MoreArgs=list(b=sigma0,c=(sigma0+sigma1)))

d <- inverse.components[1,]
e <- inverse.components[2,]
f <- inverse.components[3,]
g <- inverse.components[4,]


Info <- construct.beta.info.matrix(n0,n1,x0,x1,s20,s21,d,e,f,g)

SE <- sqrt(solve(Info)[4,4])*sqrt(weight)
power <- power(beta=effect,se=SE,alpha=alpha)

return(
       list(
       estimated.power=power,
       estimated.se=SE,
       sigma=mean(sigma),
       sigma0=sigma0,
       sigma1=sigma1,
       level=alpha)
       )
}

