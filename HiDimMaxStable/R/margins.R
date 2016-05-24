
# Different kinds of margin distributions that are used. The
# distributions must have unit mean

setClass("margin",
         slots=list(
             name="character",type="character",
             d="function",p="function",q="function",
             p.min="numeric",p.max="numeric",p.start="numeric"
             ))

marginUnif<-new("margin",
                name="Uniform",type="unif",
                d=function(x,param) dunif(x,min=1-param,max=1+param),
                p=function(q,param) punif(q,min=1-param,max=1+param),
                q=function(p,param) qunif(p,min=1-param,max=1+param),
                p.min=0,p.max=1,p.start=0.5)

marginExp<-new("margin",
               name="Exponential",type="exp",
               d=function(x,param) dexp(x,1),
               p=function(q,param) pexp(q,1),
               q=function(p,param) qexp(p,1),
               p.min=1,p.max=1,p.start=1)

marginLnorm<-new("margin",
                 name="Log-Normal",type="lnorm",
                 d=function(x,param) dlnorm(x,meanlog=-param^2/2,sdlog=param),
                 p=function(x,param) pnorm((log(x)+param^2/2)/param),
                 q=function(p,param) exp(-param^2/2+param*qnorm(p)),
                 p.min=0,p.max=Inf,p.start=0.2)

marginWeibull<-new("margin",
                   name="Weibull",type="weibull",
                   d=function(x,a) dweibull(x, shape = a, scale = 1/gamma(1+1/a)),
                   p=function(q,a) pweibull(q, shape = a, scale = 1/gamma(1+1/a)),
                   q=function(p,a) qweibull(p, shape = a, scale = 1/gamma(1+1/a)),
                   p.min=0,p.max=Inf,p.start=0.5)

# Needs library(VGAM)
marginFrechet<-new("margin",
                   name="Frechet",type="frechet",
                   d=function(x,a) dfrechet(x, location = 0, scale = 1/gamma(1-1/a),  shape = a),
                   p=function(q,a) pfrechet(q, location = 0, scale = 1/gamma(1-1/a),  shape = a),
                   q=function(p,a) qfrechet(p, location = 0, scale = 1/gamma(1-1/a),  shape = a),
                   p.min=1,p.max=Inf,p.start=1.2)

marginGamma<-new("margin",
                 name="Gamma",type="gamma",
                 d=function(x,a) dgamma(x, shape = a, scale = 1/a),
                 p=function(q,a) pgamma(q, shape = a, scale = 1/a),
                 q=function(p,a) qgamma(p, shape = a, scale = 1/a),
                 p.min=0,p.max=Inf,p.start=1)

# Needs library(VGAM)
marginGPD<-new("margin",
               name="GPD",type="gpd",
               d=function(x,a) dgpd(x, xi = a, mu = 0, beta = 1-a),
               p=function(q,a) pgpd(q, xi = a, mu = 0, beta = 1-a),
               q=function(p,a) qgpd(p, xi = a, mu = 0, beta = 1-a),
               p.min=-Inf,p.max=1,p.start=0)
