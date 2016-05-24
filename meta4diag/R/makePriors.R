makePriors <- function(var.prior = "PC", var2.prior="PC", cor.prior = "PC",
                       var.par = c(3, 0.05), var2.par, cor.par = c(1,-0.1,0.5,-0.95,0.05,0.95,0.05),
                       wishart.par = c(4,1,1,0),
                       init = c(0.01, 0.01, -0.1)){
  
  if(requireNamespace("INLA", quietly = TRUE)){
    options(warn=-1)
    var.prior = tolower(var.prior)
    var2.prior = tolower(var2.prior)
    cor.prior = tolower(cor.prior)
    
    if(any(c(var.prior, var2.prior, cor.prior)=="invwishart")){
      if(!is.numeric(wishart.par)){
        stop("Argument \"wishart.par\" can ONLY be a numeric vector with length 4.!!!")
      }
      if(!length(wishart.par)==4){
        stop("Argument \"wishart.par\" can ONLY be a numeric vector with length 4.!!!")
      }
      wishart.flag = TRUE
      
      prec1 = list(prior = "wishart2d", param = wishart.par)
      
      if(missing(var2.par)){
        var2.par = NULL
      }
      original.setting = list(var1 = list(prior=var.prior, param=var.par, initial=init[1]),
                              var2 = list(prior=var2.prior, param=var2.par, initial=init[2]),
                              cor = list(prior=cor.prior, param=cor.par, initial=init[3]),
                              wishart.par = wishart.par)
      
      priors = list(prec1 = prec1,
                    original.setting = original.setting, wishart.flag=wishart.flag)
    }else{
      wishart.flag = FALSE
      
      halfCauchy = function(gamma){
        sig = c(seq(0.00001,3,by=0.001),seq(3,10,by=0.1),seq(10,10000,by=0.1))
        sig = unique(sig)
        tau = sig^(-2)
        ltau = log(tau)
        pi_sig = 2*gamma/(pi*(sig^2+gamma^2))
        pi_variance = pi_sig*abs(0.5/sig)
        pi_precision = pi_sig*abs(0.5*tau^(-1.5))
        pi_logprecision = pi_sig*abs(0.5*sig)
        df_variance = data.frame(x = sig^2, y = pi_variance)
        df_logprecision = data.frame(x = ltau, y = log(pi_logprecision))
        df_distance = data.frame(x = sig, y = pi_variance*sig)
        return(list(variance=df_variance, logprecision=df_logprecision, distance=df_distance))
      }
      
      uniformSig = function(){
        sig = c(seq(0.00001,3,by=0.001),seq(3,10,by=0.1),seq(10,10000,by=0.1))
        sig = unique(sig)
        tau = sig^(-2)
        ltau = log(tau)
        pi_sig = rep(1/max(sig),length(sig))
        pi_variance = pi_sig*abs(0.5/sig)
        pi_precision = pi_sig*abs(0.5*tau^(-1.5))
        pi_logprecision = pi_sig*abs(0.5*sig)
        df_variance = data.frame(x = sig^2, y = pi_variance)
        df_logprecision = data.frame(x = ltau, y = log(pi_logprecision))
        df_distance = data.frame(x = sig, y = pi_variance*sig)
        return(list(variance=df_variance, logprecision=df_logprecision, distance=df_distance))
      }
      
      tablevar2lprec = function(var.par){
        sig2 = var.par[,1]
        pi_variance = var.par[,2]
        tau = 1/sig2
        ltau = log(tau)
        pi_logprecision = pi_variance*abs(-sig2)
        df = data.frame(x = ltau, y = log(pi_logprecision))
        return(df)
      }
      
      table2distV = function(vp){
        x = sqrt(vp[,1])
        y = vp[,2]*x
        return(data.frame(x = x, y = y))
      }
      
      
      ################ check var.prior and cor.prior
      if(!(var.prior %in% c("pc","invgamma","hcauchy","tnormal","unif","invwishart","table"))){
        stop("Argument \"var.prior\" can ONLY be one string. The options are \"invgamma\", \"PC\", \"hcauchy\", \"tnormal\", \"unif\", \"invwishart\" and \"table\"!!!")
      }
      if(length(var.prior)!=1){
        stop("Argument \"var.prior\" can ONLY be one string. The options are \"invgamma\", \"PC\", \"hcauchy\", \"tnormal\", \"unif\", \"invwishart\" and \"table\"!!!")
      }
      if(!(var2.prior %in% c("pc","invgamma","hcauchy","tnormal","unif","invwishart","table"))){
        stop("Argument \"var2.prior\" can ONLY be one string. The options are \"invgamma\", \"PC\", \"hcauchy\", \"tnormal\", \"unif\", \"invwishart\" and \"table\"!!!")
      }
      if(length(var2.prior)!=1){
        stop("Argument \"var2.prior\" can ONLY be one string. The options are \"invgamma\", \"PC\", \"hcauchy\", \"tnormal\", \"unif\", \"invwishart\" and \"table\"!!!")
      }
      if(!(cor.prior %in% c("pc","normal","beta","invwishart","table"))){
        stop("Argument \"cor.prior\" can ONLY be one string. The options are \"normal\", \"PC\", \"beta\", \"invwishart\" and \"table\"!!!")
      }
      if(length(cor.prior)!=1){
        stop("Argument \"cor.prior\" can ONLY be one string. The options are \"normal\", \"PC\", \"beta\", \"invwishart\" and \"table\"!!!")
      }
      
      ################ check var.par and cor.par
      if(var.prior=="pc"){
        if(!is.numeric(var.par)){
          stop("Argument \"var.par\" can ONLY be numerical vector!!!")
        }
        if(length(var.par)!=2){
          stop("Argument \"var.par\" should be c(u, alpha)!!!")
        }
        if(any(var.par<=0)){
          stop("Each element of argument \"var.par\"  MUST be positive!!!")
        }
      }
      
      if(var.prior=="invgamma"){
        if(!is.numeric(var.par)){
          stop("Argument \"var.par\" can ONLY be numerical vector!!!")
        }
        if(length(var.par)!=2){
          stop("Argument \"var.par\" should be c(a, b)!!!")
        }
        if(any(var.par<=0)){
          stop("Each element of argument \"var.par\"  MUST be positive!!!")
        }
      }
      
      if(var.prior=="tnormal"){
        if(!is.numeric(var.par)){
          stop("Argument \"var.par\" can ONLY be numerical vector!!!")
        }
        if(length(var.par)!=2){
          stop("Argument \"var.par\" should be c(m, v)!!!")
        }
        if(var.par[2]<=0){
          stop("The second element of argument \"var.par\"  MUST be positive!!!")
        }
      }
      
      if(var.prior=="hcauchy"){
        if(!is.numeric(var.par)){
          stop("Argument \"var.par\" can ONLY be numerical vector!!!")
        }
        if(length(var.par)!=1){
          stop("Argument \"var.par\" should be c(gamma)!!!")
        }
        if(var.par<=0){
          stop("The element of argument \"var.par\"  MUST be positive!!!")
        }
      }
      
      if(var.prior=="table"){
        if(!is.data.frame(var.par)){
          stop("Argument \"var.par\" can ONLY be a data frame with 2 columns!!!")
        }
        if(any(is.na(var.par))){
          stop("Argument \"var.par\" can not contain NA values!!!")
        }
        if(any(var.par)<0){
          stop("The element of argument \"var.par\"  MUST be non-negative!!!")
        }
      }
      
      
      if(missing(var2.par)){
        var2.par=var.par
      }else{
        if(var2.prior=="pc"){
          if(!is.numeric(var2.par)){
            stop("Argument \"var2.par\" can ONLY be numerical vector!!!")
          }
          if(length(var2.par)!=2){
            stop("Argument \"var2.par\" should be c(u, alpha)!!!")
          }
          if(any(var2.par<=0)){
            stop("Each element of argument \"var2.par\"  MUST be positive!!!")
          }
        }
        
        if(var2.prior=="invgamma"){
          if(!is.numeric(var2.par)){
            stop("Argument \"var2.par\" can ONLY be numerical vector!!!")
          }
          if(length(var2.par)!=2){
            stop("Argument \"var2.par\" should be c(a, b)!!!")
          }
          if(any(var2.par<=0)){
            stop("Each element of argument \"var2.par\"  MUST be positive!!!")
          }
        }
        
        if(var2.prior=="tnormal"){
          if(!is.numeric(var2.par)){
            stop("Argument \"var2.par\" can ONLY be numerical vector!!!")
          }
          if(length(var2.par)!=2){
            stop("Argument \"var2.par\" should be c(m, v)!!!")
          }
          if(var2.par[2]<=0){
            stop("The second element of argument \"var2.par\"  MUST be positive!!!")
          }
        }
        
        if(var2.prior=="hcauchy"){
          if(!is.numeric(var2.par)){
            stop("Argument \"var2.par\" can ONLY be numerical vector!!!")
          }
          if(length(var2.par)!=1){
            stop("Argument \"var2.par\" should be c(gamma)!!!")
          }
          if(var2.par<=0){
            stop("The element of argument \"var2.par\"  MUST be positive!!!")
          }
        }
        
        if(var2.prior=="table"){
          if(!is.data.frame(var2.par)){
            stop("Argument \"var2.par\" can ONLY be a data frame with 2 columns!!!")
          }
          if(any(is.na(var2.par))){
            stop("Argument \"var2.par\" can not contain NA values!!!")
          }
          if(any(var2.par)<0){
            stop("The element of argument \"var2.par\"  MUST be non-negative!!!")
          }
        }
      }
      
      if(cor.prior=="pc"){
        if(!is.numeric(cor.par)){
          stop("Argument \"cor.par\" can ONLY be numerical vector!!!")
        }
        if(length(cor.par)!=7){
          stop("Argument \"cor.par\" can ONLY be numerical vector with length 7 when cor.prior is \"PC\"!!!")
        }
      }
      
      if(cor.prior=="normal"){
        if(!is.numeric(cor.par)){
          stop("Argument \"cor.par\" can ONLY be numerical vector!!!")
        }
        if(length(cor.par)!=2){
          stop("Argument \"cor.par\" can ONLY be numerical vector with length 2 when cor.prior is \"normal\"!!!")
        }
        if(cor.par[2]<=0){
          stop("The second element of cor.par must be positive!!!")
        }
      }
      if(cor.prior=="table"){
        if(!is.data.frame(cor.par)){
          stop("Argument \"cor.par\" can ONLY be a data frame with 2 columns!!!")
        }
        if(any(is.na(cor.par))){
          stop("Argument \"cor.par\" can not contain NA values!!!")
        }
        if(any(cor.par[,2])<0){
          stop("The second column of argument \"cor.par\"  MUST be non-negative!!!")
        }
      }
      
      ################ check init
      if(!is.numeric(init)){
        stop("Argument \"init\" can ONLY be numerical vector with length 3!!!")
      }
      if(length(init)!=3){
        stop("Argument \"init\" can ONLY be numerical vector with length 3!!!")
      }
      if(init[1]<=0){
        stop("Initial value for the first variance component MUST be positive!!!")
      }
      if(init[2]<=0){
        stop("Initial value for the second variance component MUST be positive!!!")
      }
      if(init[3]<=-1 || init[3]>=1){
        stop("Initial value for the correlation component MUST be between -1 and 1 !!!")
      }
      if(init[3]==0){init[3]=0.001}
      
      ################ transform initial value
      to.theta = function(x) log((1+x)/(1-x))
      from.theta = function(x) 2*exp(x)/(1+exp(x))-1
      prec1.init = log(1/init[1])
      prec2.init = log(1/init[2])
      rho.init = .to(init[3])
      
      density = list()
      density.distance = list()
      ############################## start making
      if(var.prior=="invgamma"){
        # cat("Inverse gamma prior for variance!  \n")
        prec1=list(prior = "loggamma",param = var.par, initial = prec1.init)
        density$variance = .priorGammaV(a = var.par[1], b = var.par[2])
        density.distance$variance = .distGammaV(a = var.par[1], b = var.par[2])
      } else if(var.prior=="pc"){
        # cat("PC prior for variance! \n")
        prec1=list(prior = "pc.prec",param = var.par, initial = prec1.init)
        density$variance = .priorPCV(u = var.par[1], alpha = var.par[2])
        density.distance$variance = .distPCV(u = var.par[1], alpha = var.par[2])
      } else if(var.prior=="hcauchy"){
        priorlist = halfCauchy(var.par)
        #prior.prec = INLA::inla.tmarginal(function(x) log(1/x), priorlist$variance, n=2048)
        prior.table = paste(c("table:",cbind(priorlist$logprecision$x, priorlist$logprecision$y)), sep="", collapse=" ")
        prec1 = list(prior = prior.table, initial = prec1.init)
        density$variance = priorlist$variance
        density.distance$variance = priorlist$distance
      } else if(var.prior=="unif"){
        priorlist = uniformSig()
        #prior.prec = INLA::inla.tmarginal(function(x) log(1/x), priorlist$variance, n=2048)
        prior.table = paste(c("table:",cbind(priorlist$logprecision$x, priorlist$logprecision$y)), sep="", collapse=" ")
        prec1 = list(prior = prior.table, initial = prec1.init)
        density$variance = priorlist$variance
        density.distance$variance = priorlist$distance
      } else if(var.prior=="tnormal"){
        prec1=list(prior = "logtnormal",param = c(var.par[1],1/var.par[2]), initial = prec1.init)
        density$variance = .priorTnormV(m = var.par[1], v = var.par[2])
        density.distance$variance = .distTnormV(m = var.par[1], v = var.par[2])
      } else if(var.prior=="table"){
        tau.prior = tablevar2lprec(var.par)
        #tau.prior = INLA::inla.tmarginal(function(x) log(1/x), var.par, n=1024)
        prior.table = paste(c("table:",cbind(tau.prior$x, tau.prior$y)), sep="", collapse=" ")
        prec1 = list(prior = prior.table, initial = prec1.init)
        density$variance = var.par
        density.distance$variance = table2distV(var.par)
      }
      
      if(var2.prior=="invgamma"){
        prec2=list(prior = "loggamma",param = var2.par, initial = prec2.init)
        density$variance2 = .priorGammaV(a = var2.par[1], b = var2.par[2])
        density.distance$variance2 = .distGammaV(a = var2.par[1], b = var2.par[2])
      } else if(var2.prior=="pc"){
        prec2=list(prior = "pc.prec",param = var2.par, initial = prec2.init)
        density$variance2 = .priorPCV(u = var2.par[1], alpha = var2.par[2])
        density.distance$variance2 = .distPCV(u = var2.par[1], alpha = var2.par[2])
      } else if(var2.prior=="hcauchy"){
        priorlist = halfCauchy(var2.par)
        #prior.prec = INLA::inla.tmarginal(function(x) log(1/x), priorlist$variance, n=2048)
        prior.table = paste(c("table:",cbind(priorlist$logprecision$x, priorlist$logprecision$y)), sep="", collapse=" ")
        prec2 = list(prior = prior.table, initial = prec2.init)
        density$variance2 = priorlist$variance
        density.distance$variance2 = priorlist$distance
      } else if(var2.prior=="unif"){
        priorlist = uniformSig()
        #prior.prec = INLA::inla.tmarginal(function(x) log(1/x), priorlist$variance, n=2048)
        prior.table = paste(c("table:",cbind(priorlist$logprecision$x, priorlist$logprecision$y)), sep="", collapse=" ")
        prec2 = list(prior = prior.table, initial = prec2.init)
        density$variance2 = priorlist$variance
        density.distance$variance2 = priorlist$distance
      } else if(var2.prior=="tnormal"){
        prec2=list(prior = "logtnormal",param = c(var2.par[1],1/var2.par[2]), initial = prec2.init)
        density$variance2 = .priorTnormV(m = var2.par[1], v = var2.par[2])
        density.distance$variance2 = .distTnormV(m = var2.par[1], v = var2.par[2])
      } else if(var2.prior=="table"){
        tau.prior = tablevar2lprec(var2.par)
        #tau.prior = INLA::inla.tmarginal(function(x) log(1/x), var2.par, n=1024)
        prior.table = paste(c("table:",cbind(tau.prior$x, tau.prior$y)), sep="", collapse=" ")
        prec2 = list(prior = prior.table, initial = prec2.init)
        density$variance2 = var2.par
        density.distance$variance2 = table2distV(var2.par)
      }
      
      kld =function(rho, rho.ref) {(1-rho*rho.ref)/(1-rho.ref^2)-0.5*log((1-rho^2)/(1-rho.ref^2))-1}
      d = function(rho, rho.ref) {sqrt(2*kld(rho, rho.ref))}
      partial = function(rho,rho.ref){(rho/(1-rho^2)-rho.ref/(1-rho.ref^2))/d(rho,rho.ref)}
      
      if(cor.prior=="normal"){
        if(length(cor.par)==2){
          cor = list(prior = "normal", param = c(cor.par[1],1/cor.par[2]), initial = rho.init)
          density$correlation = .priorRhoNormal(mean = cor.par[1], variance = cor.par[2])
          rho_vec = density$cor$x
          dist = d(rho_vec, cor.par[1])*(rho_vec >= cor.par[1]) - d(rho_vec, cor.par[1])*(rho_vec < cor.par[1])
          pi.d = density$cor$y/abs(partial(rho_vec, cor.par[1]))
          density.distance$cor = data.frame(x=dist, y=pi.d)
        } else{
          stop("Argument \"cor.par\" should be length 2 when \"cor.prior\" is normal !!!")
        }
      } else if(cor.prior=="pc"){
        if(length(cor.par)==7){
          rho.ref = cor.par[2]
          left.portion = cor.par[3]
          Umin = cor.par[4]
          alpha1 = cor.par[5]
          Umax = cor.par[6]
          alpha2 = cor.par[7]
          strategy = cor.par[1]
          
          if(!(strategy %in% c(1,2,3))){
            stop("The first element of \"cor.par\" is the Strategy we use to construct PC prior for correlation. Strategy must be given, and can only be either 1, 2 or 3!!!")
          }
          if(is.na(rho.ref)){
            stop("The second element of \"cor.par\" is the reference value of correlation, which MUST be given!!!")
          }
          if(strategy==1){
            #cat("cor.par = c(1, rho.ref, left.portion, u1, alpha1, NA, NA). \n")
            Umax = NULL
            alpha2 = NULL
            if(is.na(left.portion) | is.na(Umin) | is.na(alpha1)){
              stop("left.portion, u1 and alpha1 must be given!!!")
            }
          }
          if(strategy==2){
            #cat("cor.par = c(2, rho.ref, left.portion, NA, NA, u2, alpha2). \n")
            Umin = NULL
            alpha1 = NULL
            if(is.na(left.portion) | is.na(Umax) | is.na(alpha2)){
              stop("left.portion, u2 and alpha2 must be given!!!")
            }
          }
          if(strategy==3){
            #cat("cor.par = c(3, rho.ref, NA, u1, alpha1, u2, alpha2). \n")
            left.portion = NULL
            if(is.na(Umin) | is.na(alpha1) | is.na(Umax) | is.na(alpha2)){
              stop("u1, alpha1, u2 and alpha2 must be given!!!")
            }
          }
          
          res = .findLambdas(rho.ref=rho.ref,left.portion=left.portion,Umin=Umin,alpha1=alpha1,Umax=Umax,alpha2=alpha2,density.name="exp",plot.flag=FALSE)
          theta.prior = .thetaLPDF(rho.ref=rho.ref, lambda1=res$lambda1, lambda2=res$lambda2, density.name="exp")
          prior.table = paste(c("table:",cbind(theta.prior$theta, theta.prior$tldv)), sep="", collapse=" ")
          cor=list(prior = prior.table,
                   initial = rho.init,
                   to.theta = function(x) log((1+x)/(1-x)),
                   from.theta = function(x) 2*exp(x)/(1+exp(x))-1)
          
          density$correlation = .rhoPDF(rho.ref=rho.ref, lambda1=res$lambda1, lambda2=res$lambda2, density.name="exp", plot.flag=FALSE)
          rho_vec = density$cor$x
          dist = d(rho_vec, rho.ref)*(rho_vec >= rho.ref) - d(rho_vec, rho.ref)*(rho_vec < rho.ref)
          pi.d = density$cor$y/abs(partial(rho_vec, rho.ref))
          density.distance$cor = data.frame(x=dist, y=pi.d)
        }else{
          stop("Argument \"cor.par\" should be length 7 when \"cor.prior\" is PC !!!")
        }
      } else if(cor.prior=="beta"){
        if(length(cor.par)==2){
          cor = list(prior = "betacorrelation", param = c(cor.par[1], cor.par[2]), initial = rho.init)
          density$correlation = .priorRhoBeta(a = cor.par[1], b = cor.par[2])
        } else{
          stop("Argument \"cor.par\" should be length 2 when \"cor.prior\" is beta !!!")
        }
      } else if(cor.prior=="table"){
        theta.prior = INLA::inla.tmarginal(function(x) log((1+x)/(1-x)), cor.par, n=1024)
        prior.table = paste(c("table:",cbind(theta.prior$x, log(theta.prior$y))), sep="", collapse=" ")
        cor = list(prior = prior.table, initial = rho.init)
        density$correlation = cor.par
      }
      original.setting = list(var1 = list(prior=var.prior, param=var.par,initial=init[1]),
                              var2 = list(prior=var2.prior, param=var2.par,initial=init[2]),
                              cor = list(prior=cor.prior, param=cor.par,initial=init[3]),
                              wishart.par = wishart.par)
      if(cor.prior=="pc"){
        priors = list(prec1 = prec1, prec2 = prec2, cor = cor, 
                      lambdas = c(res$lambda1, res$lambda2), density = density, 
                      original.setting = original.setting, wishart.flag=wishart.flag)
      }else{
        priors = list(prec1 = prec1, prec2 = prec2, cor = cor, 
                      density = density, 
                      original.setting = original.setting, wishart.flag=wishart.flag)
      }
    }
    options(warn=0)
    return(priors)
  }else{
    stop("INLA need to be installed and loaded!\n
         Please use the following commants to install and load INLA,\n
         install.packages(\"INLA\", repos=\"http://www.math.ntnu.no/inla/R/testing\")
         library(INLA) \n")
  }
}

.priorRhoNormal <- function(mean,variance){
  # get the density of old normal prior of correlation parameter "rho"
  rho=c(seq(-1,-0.9,len=500),seq(-0.9001,0.8999,len=500),seq(0.9,1,len=500))
  transf = function(rho){.logit(0.5*rho+0.5)}
  z = transf(rho)
  dens = dnorm(z,mean=mean,sd=sqrt(variance))*abs(2/(1-rho^2))
  return(data.frame(x = rho, y=dens))
}

.priorRhoBeta <- function(a,b){
  # get the density of old normal prior of correlation parameter "rho"
  rho=c(seq(-1,-0.9,len=500),seq(-0.9001,0.8999,len=500),seq(0.9,1,len=500))
  z = 0.5*(rho+1)
  dens = 0.5*dbeta(z,shape1=a,shape2=b)
  return(data.frame(x = rho, y=dens))
}

.priorGammaV <- function(a,b){
  x = seq(0,7,len=100000)
  tau = 1/x
  y = dgamma(tau,shape=a,rate=b)*abs(-x^(-2))
  return(data.frame(x=x,y=y))
}

.priorPCV <- function(u,alpha){
  x = c(seq(0.000001,1,len=1000),seq(1.0001,100,len=100))
  tau = 1/x
  theta = -log(alpha)/u
  y = 0.5*theta*tau^(-1.5)*exp(-theta/sqrt(tau))*abs(-x^(-2))
  return(data.frame(x=x,y=y))
}

.distGammaV <- function(a,b){
  x = seq(0,7,len=100000)
  tau = x^(-2)
  y = dgamma(tau,shape=a,rate=b)*abs(-2/(x^3))
  return(data.frame(x=x,y=y))
}

.distPCV <- function(u,alpha){
  x = seq(0,7,len=100)
  theta = -log(alpha)/u
  y = theta*exp(-theta*x)
  return(data.frame(x=x,y=y))
}

.priorTnormV <- function(m, v){
  par.sig = sqrt(v)
  x = seq(0,7,len=1000)
  var.sig = sqrt(x)
  pi.var.sig = 1/par.sig*dnorm((var.sig-m)/par.sig)/(1-pnorm(-m/par.sig))
  y = 0.5*pi.var.sig/var.sig
  return(data.frame(x = x, y = y))
}

.distTnormV <- function(m, v){
  p = .priorTnormV(m, v)
  x = sqrt(p$x)
  y = p$y*x
  return(data.frame(x = x, y = y))
}
