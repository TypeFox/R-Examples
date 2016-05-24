##  Odds Ratio and Relative Risks for 2 x 2 Tables Odds Ratio
##  (Case-Control Studies)

## The odds ratio is a useful measure of association for a variety of
## study designs. For a retrospective design called a case- control
## study, the odds ratio can be used to estimate the relative risk
## when the probability of positive response is small (Agresti
## 2002). In a case-control study, two independent samples are
## identified based on a binary (yes-no) response variable, and the
## conditional distribution of a binary explanatory variable is
## examined, within fixed levels of the response variable. See Stokes,
## Davis, and Koch (2000) and Agresti (2007).

oddsratio <- function(x,alpha=0.05,n, ...)
UseMethod("oddsratio")

oddsratio.default <- function(x, alpha=0.05,n,...)
{
  if(length(x) != 2 || length(n) != 2)
    stop("Lengths of 'x' and/or 'n' not equal to two.")
  if(any(x > n) || any(x<0) || any(n<0))
    stop("Invalid data")
  n = as.integer(n); x = as.integer(x)
  y = data.frame(Event = x, Size = n)
  row.names(y) <- c("Sample 1","Sample 2")
  oddsratio(as.data.frame(y),alpha=alpha,...)
}


oddsratio.matrix <- function(x,alpha=0.05,...)
{
  oddsratio(as.data.frame(x),alpha=alpha,...)
}

oddsratio.data.frame <- function(x,alpha=0.05,...)
{
  stopifnot(class(x)=='data.frame')
  stopifnot(alpha>0&&alpha<1)
  n11 = x[1,1]; n1 = x[1,2]; n12 = n1-n11
  n21 = x[2,1]; n2 = x[2,2]; n22 = n2-n21
  if(n12<0||n22<0) stop("Invalid data!")
  OR = n11*n22/(n12*n21)
  ##  var(log(OR))
  CI = NULL
  LL = NULL; UL = NULL;
  if(n11==0 || n12==0 || n21==0||n22 == 0){
    warning("Asymptotic confidence limits for OR are not computed.")
  }else{
    v <- 1/n11+1/n21+1/n12+1/n22
    ## v = sum(1/x); # Thanks mush go to Dr. Mike Schaffer for finding this bug. 
    z = abs(qnorm(alpha/2))
    LL = OR * exp(-z*sqrt(v)); UL = OR * exp(z*sqrt(v));
    Asymptotic = c(LL = LL, UL = UL)
    CI = rbind(CI, Asymptotic)
  }
  ## compute the exact CI
  if(n11==0 || n22 == 0){
    LL = 0
    UL = .Fortran(.F_orexactu,
      as.integer(c(n11,n12,n21,n22)),
      as.double(alpha), ll=as.double(OR))$ll;
  }else if(n12 == 0 || n21 == 0){
    UL = Inf
    LL = .Fortran(.F_orexactl,
      as.integer(c(n11,n12,n21,n22)),
      as.double(alpha), ll=as.double(1000))$ll;
  }else{
    LL = .Fortran(.F_orexactl,
      as.integer(c(n11,n12,n21,n22)),
      as.double(alpha), ll=as.double(OR))$ll;
    UL = .Fortran(.F_orexactu,
      as.integer(c(n11,n12,n21,n22)),
      as.double(alpha), ll=as.double(OR))$ll;
  }
  
  Exact = c(LL = LL, UL = UL)
  CI = rbind(CI, Exact)

  tmp = .orscoreci(n11,n1,n21,n2,1-alpha)
  Score = c(LL = tmp[1], UL = tmp[2])
  CI = rbind(CI, Score)
  
  ## we also compute the relative risks and the corresponding
  ## confience limits.
  RR = n11*n2/(n1*n21)
  ## v = var(log(RR))
  CI2 = NULL
  p1 = n21/n2;
  LL = ifelse(is.finite(LL),LL/((1-p1)+p1*LL),1/p1)
  UL = ifelse(is.finite(UL), UL/((1-p1)+p1*UL),1/p1)
  transform = c(LL = LL, UL = UL)
  CI2 = rbind(CI2, transform)
  if(n11==0 || n21==0){
    warning("Asymptotic confidence limits for RR are not computed.")
  }else{
    v = n21/(n11*n1) + n22/(n21*n2); z = abs(qnorm(alpha/2))
    LL = RR * exp(-z*sqrt(v)); UL = RR * exp(z*sqrt(v));
    Asymptotic = c(LL = LL, UL = UL)
    CI2 = rbind(CI2, Asymptotic)
  }
  tmp = .riskscoreci(n11,n1,n21,n2,1-alpha)
  Score = c(LL = tmp[1], UL = tmp[2])
  CI2 = rbind(CI2, Score)

  ##  adverse <- data.frame(r = c(n11,n21), n = c(n1,n2), treat=c(1,0))
  ##  LL = uniroot(.prof.uniroot, interval=c(-2,10),
  ##    adverse=adverse,treat=treat)$root
  ##  Profile = c(LL=exp(LL), UL=Score[2])
  ##  CI2 = rbind(CI2, Profile)

  structure(list(OR=OR, ORCI=CI, RR=RR,RRCI = CI2,
                 alpha=alpha,data = x),
            class = "odds")
}

print.odds <- function (x, digits = NULL, ...) 
{
  cat("\nData:\n")
  print(x$data)
  cat("\nOdds ratio:\t", format(x$OR, just='centre', digits=digits))
  cat("\n", 100*(1-x$alpha),"% confidence intervals\n")
  print(x$ORCI, digits = digits, ...)
  cat("\nRelative risk:\t", format(x$RR, just='centre', digits=digits))
  cat("\n", 100*(1-x$alpha),"% confidence intervals\n")
  print(x$RRCI, digits = digits, ...)
  cat("\n\n")
  invisible(x)
}

############################################
##  The following codes are selected on Dec 3, 2011. Code has been
## published for generating confidence intervals by inverting a score
## test. It is available from http://web.stat.ufl.edu/~aa/cda/R/two_
## sample/R2/

## Large-sample score confidence interval for a difference of
## proportions in 2x2 table (Mee 1984, Miettinen and Nurminen 1985,
## Nurminen 1986)

.diffscoreci <- function(x1,n1,x2,n2,conflev){ 
   px = x1/n1
   py = x2/n2
   z = qchisq(conflev,1)
   proot = px-py
   dp = 1-proot
   niter = 1
   while(niter <= 50){
     dp = 0.5*dp
     up2 = proot+dp
     score = .z2stat(px,n1,py,n2,up2)
     if(score<z){ proot = up2 }
     niter = niter+1
     if((dp<0.0000001) || (abs(z-score)<.000001)){
       niter = 51
       ul = up2}
    } 
   
   proot = px-py
   dp = 1+proot
   niter = 1
   while(niter <= 50){
     dp = 0.5*dp
     low2 = proot-dp
     score = .z2stat(px,n1,py,n2,low2)
     if(score<z){ proot = low2 }
     niter = niter+1
     if((dp<0.0000001) || (abs(z-score)<.000001)){
     ll = low2
     niter = 51}
     }
 c(ll,ul)
}

.z2stat <- function (p1x,nx,p1y,ny,dif){

      diff = p1x-p1y-dif
      if ( abs(diff) == 0 ) {
        fmdiff = 0}
      else{
        t = ny/nx
        a = 1+t
        b = -(1+ t + p1x + t*p1y + dif*(t+2))
        c = dif*dif + dif*(2*p1x + t +1) + p1x + t*p1y
        d = -p1x*dif*(1+dif)
        v = (b/a/3)^3 - b*c/(6*a*a) + d/a/2
        s = sqrt( (b/a/3)^2 - c/a/3)
        if(v>0){u=s}
         else{u=-s}
        w = (3.141592654+acos(v/u^3))/3
        p1d = 2*u*cos(w) - b/a/3
        p2d = p1d - dif
        var = p1d*(1-p1d)/nx + p2d*(1-p2d)/ny
        fmdiff = diff^2/var
      }
  return(fmdiff)
}


# Computes the score two-sample interval for (p1 - p2)
# with success counts x1, x2 and trials n1, n2
# and with confidence coefficient = conflev 
# written by Yongyi Min
#-----------------------------------------------------------------------------
#Wald confidence interval for difference of proportions in 2x2 table
#-----------------------------------------------------------------------------
.wald2ci <- function(x1, n1, x2, n2, conflev){
   p1hat = x1/n1
   p2hat = x2/n2
   z = abs(qnorm((1-conflev)/2))
   ll = (p1hat - p2hat)  - z*sqrt((p1hat*(1-p1hat))/n1 + (p2hat*(1-p2hat))/n2)
   ul = (p1hat - p2hat)  + z*sqrt((p1hat*(1-p1hat))/n1 + (p2hat*(1-p2hat))/n2)
   c(ll,ul)
}

#Computes the Wald two-sample interval for (p1 - p2)
#with success counts x1, x2 and trials n1, n2
#and with confidence coeff conflev.  Returns a
#list with the lower and upper end points
#in that order.  Add 1 to x1 and x2 and add 2 to n1 and n2
#to get Agresti-Caffo adjusted CI (American Statistician, 2000).

##-----------------------------------------------------------------------------
## R Code for large-sample score confidence interval for a relative
## risk in a 2x2 table (Koopman 1984, Miettinen and Nurminen 1985,
## Nurminen 1986).

## -----------------------------------------------------------------------------
.riskscoreci <- function(x1,n1,x2,n2,conflev)
{
  z =  abs(qnorm((1-conflev)/2))
  if ((x2==0) &&(x1==0)){
    ul = Inf
    ll = 0
    }
  else{  
     a1 =  n2*(n2*(n2+n1)*x1+n1*(n2+x1)*(z^2))
     a2 = -n2*(n2*n1*(x2+x1)+2*(n2+n1)*x2*x1+n1*(n2+x2+2*x1)*(z^2))  
     a3 = 2*n2*n1*x2*(x2+x1)+(n2+n1)*(x2^2)*x1+n2*n1*(x2+x1)*(z^2)
     a4 = -n1*(x2^2)*(x2+x1)
     b1 = a2/a1
     b2 = a3/a1
     b3 = a4/a1
     c1 = b2-(b1^2)/3
     c2 = b3-b1*b2/3+2*(b1^3)/27
     ceta = acos(sqrt(27)*c2/(2*c1*sqrt(-c1)))
     t1 = -2*sqrt(-c1/3)*cos(pi/3-ceta/3)
     t2 = -2*sqrt(-c1/3)*cos(pi/3+ceta/3)
     t3 = 2*sqrt(-c1/3)*cos(ceta/3)
     p01 = t1-b1/3
     p02 = t2-b1/3
     p03 = t3-b1/3
     p0sum = p01+p02+p03
     p0up = min(p01,p02,p03)
     p0low = p0sum-p0up-max(p01,p02,p03)
      
     if( (x2==0) && (x1!=0) ){
        ll = (1-(n1-x1)*(1-p0low)/(x2+n1-(n2+n1)*p0low))/p0low 
        ul = Inf 
       }
     else if( (x2!=n2) && (x1==0)){
        ul = (1-(n1-x1)*(1-p0up)/(x2+n1-(n2+n1)*p0up))/p0up
        ll = 0
        }
     else if( (x2==n2) && (x1==n1)){
         ul = (n2+z^2)/n2
         ll =  n1/(n1+z^2)
        }
     else if( (x1==n1) || (x2==n2) ){
         if((x2==n2) && (x1==0)) { ll = 0 }
         if((x2==n2) && (x1!=0)) {
           phat1  = x2/n2
           phat2  =  x1/n1
           phihat = phat2/phat1
           phil = 0.95*phihat
           chi2 = 0
           while (chi2 <= z){
             a = (n2+n1)*phil
             b = -((x2+n1)*phil+x1+n2)
             c = x2+x1
             p1hat = (-b-sqrt(b^2-4*a*c))/(2*a)
             p2hat = p1hat*phil
             q2hat = 1-p2hat
             var = (n2*n1*p2hat)/(n1*(phil-p2hat)+n2*q2hat)
             chi2 = ((x1-n1*p2hat)/q2hat)/sqrt(var)
             ll = phil
             phil = ll/1.0001}} 
         i = x2
         j = x1
         ni = n2 
         nj = n1 
         if( x1==n1 ){               
            i = x1
            j = x2
            ni = n1 
            nj = n2
         } 
         phat1  = i/ni
         phat2  =  j/nj
         phihat = phat2/phat1
         phiu = 1.1*phihat
         if((x2==n2) && (x1==0)) { 
            if(n2<100) {phiu = .01}
            else {phiu=0.001}
           } 
         chi1 = 0
         while (chi1 >= -z){
         a = (ni+nj)*phiu
         b = -((i+nj)*phiu+j+ni)
         c = i+j
         p1hat = (-b-sqrt(b^2-4*a*c))/(2*a)
         p2hat = p1hat*phiu
         q2hat = 1-p2hat
         var = (ni*nj*p2hat)/(nj*(phiu-p2hat)+ni*q2hat)
         chi1  = ((j-nj*p2hat)/q2hat)/sqrt(var)
         phiu1 = phiu
         phiu = 1.0001*phiu1
         }

         if(x1==n1) {
          ul = (1-(n1-x1)*(1-p0up)/(x2+n1-(n2+n1)*p0up))/p0up  
          ll = 1/phiu1       
         }
         else{ ul = phiu1}                        
       }   

     else{
     ul = (1-(n1-x1)*(1-p0up)/(x2+n1-(n2+n1)*p0up))/p0up
     ll = (1-(n1-x1)*(1-p0low)/(x2+n1-(n2+n1)*p0low))/p0low 
      }
   }  
  c(ll,ul)
}


# computes the score CI for the relative risk p1/p2
# with success counts x1, x2 and trials n1, n2
# and with confidence coefficient = conflev 
# written by Yongyi Min
## -----------------------------------------------------------------------------
##  R Code for large-sample score confidence interval for an odds
##  ratio in a 2x2 table (Cornfield 1956, Miettinen and Nurminen 1985)

## -----------------------------------------------------------------------------
.orscoreci <- function(x1,n1,x2,n2,conflev){
  px = x1/n1
  py = x2/n2
  if(((x1==0) && (x2==0)) || ((x1==n1) && (x2==n2))){
      ul = 1/0
      ll = 0   
      } 
  else if((x1==0) || (x2==n2)){
       ll = 0
       theta = 0.01/n2 
       ul = .limit(x1,n1,x2,n2,conflev,theta,1)      
     }
  else if((x1==n1) || (x2==0)){
       ul = 1/0
       theta = 100*n1
       ll = .limit(x1,n1,x2,n2,conflev,theta,0)       
     }
  else{
      theta = px/(1-px)/(py/(1-py))/1.1
      ll = .limit(x1,n1,x2,n2,conflev,theta,0)       
      theta=px/(1-px)/(py/(1-py))*1.1
      ul = .limit(x1,n1,x2,n2,conflev,theta,1)      
   }
  c(ll,ul)  
}

.limit <- function(x,nx,y,ny,conflev,lim,t){
      
      z = qchisq(conflev,1)
      px = x/nx
      score= 0
      while ( score < z){
        a = ny*(lim-1)
        b = nx*lim+ny-(x+y)*(lim-1)
        c = -(x+y)
        p2d = (-b+sqrt(b^2-4*a*c))/(2*a)
        p1d = p2d*lim/(1+p2d*(lim-1))
        score = ((nx*(px-p1d))^2)*(1/(nx*p1d*(1-p1d))+1/(ny*p2d*(1-p2d)))
        ci = lim
        if(t==0) { lim = ci/1.001 }
        else{ lim = ci*1.001 }
        } 
 return(ci)
}
# computes the score CI for the odds ratio  p1(1-p2)/(p2(1-p1))
# with success counts x1, x2 and trials n1, n2
# and with confidence coefficient = conflev 
# written by Yongyi Min


###########################################
##  The following codes are written by Michael E Dewey in "Collated
##  responses from R-help on confidence intervals for risk ratios"
##  (November, 2006).  Section of Epidemiology, Division of
##  Psychological Medicine, PO 60, Institute of Psy- chiatry, De
##  Crespigny Park, London, SE5 8AF, UK. mailto:m.dewey@iop.kcl.ac.uk
### To call, using
## adverse <- data.frame(r = c(8,0), n = c(500,500), treat = c(1,0))
## uniroot(.prof.uniroot,interval=c(-2,2))


###   need to be double checked.

.prof <- function(delta,adverse,treat) { # vectorised to enable plotting
  n.delta <- length(delta)
  dev <- numeric(n.delta)
  for (i in 1:n.delta) {
    dev[i] <- glm(r ~ offset(treat * delta[i] + log(n)),
                  family = poisson, data = adverse)$deviance
  }
  dev
}

.prof.uniroot <- function(delta,adverse,treat)
  .prof(delta,adverse,treat) - qchisq(0.95, 1)


