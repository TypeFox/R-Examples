nlsfit <-
function(data, model=1, start=c(a=1,b=1,c=1,d=1, e=1)){
    s=start
    n=names(data)
    d1 = as.data.frame(data[, 1])
    d2 = as.data.frame(data[, -1])
    f = function(h) {
        data.frame(d1, d2[h])
    }
    
    h = length(d2)
    h = 1:h
    l = lapply(h, f)
    R2 <- function(m) {
        gl <- length(fitted(m)) - 1
        sqt <- var((fitted(m) + resid(m))) * gl
        r1 <- (sqt - deviance(m))/sqt
        r1=round(r1,4)
        return(r1)
    }
    R3 <- function(m) {
        gl <- length(fitted(m)) - 1
        sqt <- var((fitted(m) + resid(m))) * gl
        r1 <- (sqt - deviance(m))/sqt
        p1 <- (gl/((gl + 1) - (length(coef(m) + 1))) * (1 - r1))
        r2 <- 1 - p1
        r2=round(r2,4)
        return(r2)
    }
    
    # linear
    f1=function(data){names(data) = c("x", "y") 
                      ml = lm(data[, 2] ~ data[, 1])
                      c1 = coef(ml)[[1]]
                      c2 = coef(ml)[[2]]
                      m=nls(y~a+b*x, data=data, start=list(a=c1,b=c2),control = nls.control(maxiter = 6000)) 
                      c=coef(m); s=summary(m); a=c[1];b=c[2];l=c(a,b,summary(m)[11][[1]][7], summary(m)[11][[1]][8], R2(m), R3(m),AIC(m), BIC(m)); l=round(l,4);l=as.data.frame(l)
                      rownames(l)=c("coefficient a", "coefficient b", 
                                    "p-value t.test for a", "p-value t.test for b", 
                                    "r-squared", "adjusted r-squared", 
                                    "AIC", "BIC")
                      return(l)}
    
    
    # quadratic
    f2=function(data){names(data) = c("x", "y") 
                      mq = lm(data[, 2] ~ data[, 1] + I(data[, 1]^2))
                      c3 = coef(mq)[[1]]
                      c4 = coef(mq)[[2]]
                      c5 = coef(mq)[[3]]
                      m=nls(y~a+b*x+c*x^2, data=data, start=list(a=c3,b=c4,c=c5),control = nls.control(maxiter = 6000));c=coef(m); s=summary(m); a=c[1];b=c[2];c=c[3];pm = a - (b^2)/(4 * c)
                      pc = -0.5 * b/c; l=c(a,b,c,summary(m)[11][[1]][10], summary(m)[11][[1]][11], summary(m)[11][[1]][12], R2(m), R3(m),AIC(m), BIC(m), pm,pc); l=round(l,4);l=as.data.frame(l)
                      rownames(l)=c("coefficient a", "coefficient b", 
                                    "coefficient c", "p-value t.test for a", "p-value t.test for b", 
                                    "p-value t.test for c", "r-squared", "adjusted r-squared", 
                                    "AIC", "BIC","maximum or minimum value for y", "critical point in x") 
                      return(l)
    }
    
    # linear plateau
    f3=function(data){names(data) = c("x", "y")
                      ml = lm(data[, 2] ~ data[, 1])
                      mq = lm(data[, 2] ~ data[, 1] + I(data[, 1]^2))
                      c1 = coef(ml)[[1]]
                      c2 = coef(ml)[[2]]
                      c3 = coef(mq)[[1]]
                      c4 = coef(mq)[[2]]
                      c5 = coef(mq)[[3]]
                      pc = -0.5 * c4/c5
                      m <- nls(y ~ a + b * (x - c) * (x <= c), start = list(a = c1, 
                                                                            b = c2, c = pc), data = data, control = nls.control(maxiter = 6000))
                      c=coef(m); a=c[1];b=c[2];c=c[3]
                      l=c(a,b,c,summary(m)[11][[1]][10], summary(m)[11][[1]][11], summary(m)[11][[1]][12], R2(m), R3(m),AIC(m), BIC(m), a,c); l=round(l,4);l=as.data.frame(l)
                      rownames(l)=c("coefficient a", "coefficient b", 
                                    "coefficient c", "p-value t.test for a", "p-value t.test for b", 
                                    "p-value t.test for c", "r-squared", "adjusted r-squared", 
                                    "AIC", "BIC","maximum or minimum value for y", "critical point in x")
                      return(l)
    }
    
    
    # quadratic plateau
    f4=function(data){names(data) = c("x", "y")
                      mq = lm(data[, 2] ~ data[, 1] + I(data[, 1]^2))
                      c3 = coef(mq)[[1]]
                      c4 = coef(mq)[[2]]
                      c5 = coef(mq)[[3]]
                      m <- nls(y ~ (a + b * x + c * I(x^2)) * (x <= -0.5 * 
                                                                   b/c) + (a + I(-b^2/(4 * c))) * (x > -0.5 * b/c), 
                               start = list(a = c3, b = c4, c = c5), data = data, 
                               control = nls.control(maxiter = 6000))
                      c=coef(m); a=c[1];b=c[2];c=c[3]
                      pmm = (a + I(-b^2/(4 * c)))
                      pcc = -0.5 * b/c
                      l=c(a,b,c,summary(m)[11][[1]][10], summary(m)[11][[1]][11], summary(m)[11][[1]][12], R2(m), R3(m),AIC(m), BIC(m), pmm,pcc); l=round(l,4);l=as.data.frame(l)
                      rownames(l)=c("coefficient a", "coefficient b", 
                                    "coefficient c", "p-value t.test for a", "p-value t.test for b", 
                                    "p-value t.test for c", "r-squared", "adjusted r-squared", 
                                    "AIC", "BIC","maximum or minimum value for y", "critical point in x") 
                      return(l)
    }
    
    # bi linear
    f5=function(data){names(data) = c("x", "y")     
                      fp=function(a,b,c,x,d){ifelse(x>=d,(a-c*d)+(b+c)*x, a+b*x)}
                      m=nls(y~fp(a,b,c,x,d), start=list(a=s[1],b=s[2],c=s[3], d=s[4]), data=data,control = nls.control(maxiter = 6000));c=coef(m); a=c[1];b=c[2];c1=c;c=c[3];d=c1[4]
                      l=c(a,b,c,d,summary(m)[11][[1]][13], summary(m)[11][[1]][14], summary(m)[11][[1]][15],summary(m)[11][[1]][16], R2(m), R3(m),AIC(m), BIC(m)); l=round(l,4);l=as.data.frame(l)
                      rownames(l)=c("coefficient a", "coefficient b", 
                                    "coefficient c","coefficient d", "p-value t.test for a", "p-value t.test for b", 
                                    "p-value t.test for c", "p-value t.test for d","r-squared", "adjusted r-squared", 
                                    "AIC", "BIC") 
                      return(l)
    }
    
    # exponential
    f6=function(data){names(data) = c("x", "y")
                      m=nls(y~a*exp(b*x) ,start=list(a=s[1],b=s[2]),data=data,control = nls.control(maxiter = 6000));c=coef(m); a=c[1];b=c[2]
                      l=c(a,b,summary(m)[11][[1]][7], summary(m)[11][[1]][8], R2(m), R3(m),AIC(m), BIC(m)); l=round(l,4);l=as.data.frame(l)
                      rownames(l)=c("coefficient a", "coefficient b", 
                                    "p-value t.test for a", "p-value t.test for b", 
                                    "r-squared", "adjusted r-squared", 
                                    "AIC", "BIC")
                      return(l)
    } 
    
    # logistic model
    f7=function(data){names(data) = c("x", "y") 
                      m=nls(y~a*(1+b*(exp(-c*x)))^-1, data=data, start=list(a=s[1],b=s[2],c=s[3]),control = nls.control(maxiter = 6000));c=coef(m); s=summary(m); a=c[1];b=c[2];c=c[3]; l=c(a,b,c,summary(m)[11][[1]][10], summary(m)[11][[1]][11], summary(m)[11][[1]][12], R2(m), R3(m),AIC(m), BIC(m)); l=round(l,4);l=as.data.frame(l)
                      rownames(l)=c("coefficient a", "coefficient b", 
                                    "coefficient c", "p-value t.test for a", "p-value t.test for b", 
                                    "p-value t.test for c", "r-squared", "adjusted r-squared", 
                                    "AIC", "BIC") 
                      return(l)
    }
    
    # Von Bertalanffy model
    f8=function(data){names(data) = c("x", "y") 
                      m=nls(y~a*(1-b*(exp(-c*x)))^3, data=data, start=list(a=s[1],b=s[2],c=s[3]),control = nls.control(maxiter = 6000));c=coef(m); s=summary(m); a=c[1];b=c[2];c=c[3]; l=c(a,b,c,summary(m)[11][[1]][10], summary(m)[11][[1]][11], summary(m)[11][[1]][12], R2(m), R3(m),AIC(m), BIC(m)); l=round(l,4);l=as.data.frame(l)
                      rownames(l)=c("coefficient a", "coefficient b", 
                                    "coefficient c", "p-value t.test for a", "p-value t.test for b", 
                                    "p-value t.test for c", "r-squared", "adjusted r-squared", 
                                    "AIC", "BIC") 
                      return(l)
    }
    
    # Brody model
    f9=function(data){names(data) = c("x", "y") 
                      m=nls(y~a*(1-b*(exp(-c*x))), data=data, start=list(a=s[1],b=s[2],c=s[3]),control = nls.control(maxiter = 6000));c=coef(m); s=summary(m); a=c[1];b=c[2];c=c[3]; l=c(a,b,c,summary(m)[11][[1]][10], summary(m)[11][[1]][11], summary(m)[11][[1]][12], R2(m), R3(m),AIC(m), BIC(m)); l=round(l,4);l=as.data.frame(l)
                      rownames(l)=c("coefficient a", "coefficient b", 
                                    "coefficient c", "p-value t.test for a", "p-value t.test for b", 
                                    "p-value t.test for c", "r-squared", "adjusted r-squared", 
                                    "AIC", "BIC") 
                      return(l)
    }
    
    # Gompertz model
    f10=function(data){names(data) = c("x", "y") 
                       m=nls(y~a*exp(-b*exp(-c*x)), data=data, start=list(a=s[1],b=s[2],c=s[3]),control = nls.control(maxiter = 6000));c=coef(m); s=summary(m); a=c[1];b=c[2];c=c[3]; l=c(a,b,c,summary(m)[11][[1]][10], summary(m)[11][[1]][11], summary(m)[11][[1]][12], R2(m), R3(m),AIC(m), BIC(m)); l=round(l,4);l=as.data.frame(l)
                       rownames(l)=c("coefficient a", "coefficient b", 
                                     "coefficient c", "p-value t.test for a", "p-value t.test for b", 
                                     "p-value t.test for c", "r-squared", "adjusted r-squared", 
                                     "AIC", "BIC") 
                       return(l)
    }
    
    # lactation curve
    f11=function(data){names(data) = c("x", "y") 
                       m=nls(y~(a*x^b)*exp(-c*x), data=data, start=list(a=s[1],b=s[2],c=s[3]),control = nls.control(maxiter = 6000));c=coef(m); s=summary(m); a=c[1];b=c[2];c=c[3]; pm=a*(b/c)^b*exp(-b);pc=b/c;l=c(a,b,c,summary(m)[11][[1]][10], summary(m)[11][[1]][11], summary(m)[11][[1]][12], R2(m), R3(m),AIC(m), BIC(m), pm,pc); l=round(l,4);l=as.data.frame(l)
                       rownames(l)=c("coefficient a", "coefficient b", 
                                     "coefficient c", "p-value t.test for a", "p-value t.test for b", 
                                     "p-value t.test for c", "r-squared", "adjusted r-squared", 
                                     "AIC", "BIC","maximum or minimum value for y", "critical point in x") 
                       return(l)
                       # a=17;b=0.25; c=0.004
    }
    
    y1=c(25,24,26,28,30,31,27,26,25,24,23,24,22,21,22,20,21,19,18,17,18,18,16,17,15,16,14)
    x1=c(15,15,15,75,75,75,135,135,135,195,195,195,255,255,255,315,315,315,375,375,375,435,435,435,
         495,495,495)
    
    dados=data.frame(x1,y1)
    
    
    # ruminal degradation curve
    f12=function(data){names(data) = c("x", "y") 
                       m=nls(y ~ a + b * (1 - exp(-c * x)), data=data, start=list(a=20,b=60,c=0.05),control = nls.control(maxiter = 6000));c=coef(m); s=summary(m); a=c[1];b=c[2];c=c[3]; l=c(a,b,c,summary(m)[11][[1]][10], summary(m)[11][[1]][11], summary(m)[11][[1]][12], R2(m), R3(m),AIC(m), BIC(m)); l=round(l,4);l=as.data.frame(l)
                       rownames(l)=c("coefficient a", "coefficient b", 
                                     "coefficient c", "p-value t.test for a", "p-value t.test for b", 
                                     "p-value t.test for c", "r-squared", "adjusted r-squared", 
                                     "AIC", "BIC") 
                       return(l)
    }
    
    tempo=c(0,12,24,36,48,60,72,84,96,108,120,144,168,192)
    gas=c(0.002,3.8,8,14.5,16,16.5,17,17.4,17.9,18.1,18.8,19,19.2,19.3)
    
    d=data.frame(tempo,gas)
    
    m=nls(gas~(vf1/(1+exp(2-4*k1*(tempo-l))))+(vf2/(1+exp(2-4*k2*(tempo-l)))), start=list(vf1=19,vf2=4,k1=0.05,k2=0.005,l=5), data=d)
    
    
    # logistico bicompartimental
    
    f13=function(data){names(data) = c("x", "y") 
                       m=nls(y~(a/(1+exp(2-4*c*(x-e))))+(b/(1+exp(2-4*d*(x-e)))), start=list(a=s[1],b=s[2],c=s[3],d=s[4],e=s[5]), data=data,control = nls.control(maxiter = 6000));cc=coef(m); s=summary(m); a=cc[1];b=cc[2];c=cc[3];d=cc[4];e=cc[5]; l=c(a,b,c,d,e,summary(m)[11][[1]][16], summary(m)[11][[1]][17], summary(m)[11][[1]][18],summary(m)[11][[1]][19],summary(m)[11][[1]][20], R2(m), R3(m),AIC(m), BIC(m)); l=round(l,4);l=as.data.frame(l)
                       rownames(l)=c("coefficient a", "coefficient b", "coefficient c","coefficient d","coefficient e",
                                     "p-value t.test for a", "p-value t.test for b", 
                                     "p-value t.test for c","p-value t.test for d","p-value t.test for e" ,"r-squared", "adjusted r-squared", 
                                     "AIC", "BIC") 
                       return(l)
    }
    
    fun=list(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13)
    
    fu=fun[[model]]
    r1=as.data.frame(lapply(l, fu))
    names(r1)=n[-1]
    models=c("y~a+b*x","y~a+b*x+c*x^2","y ~ a + b * (x - c) * (x <= c)","y ~ (a + b * x + c * I(x^2)) * (x <= -0.5 * b/c) + (a + I(-b^2/(4 * c))) * (x > -0.5 * b/c)","ifelse(x>=d,(a-c*d)+(b+c)*x, a+b*x)","y~a*exp(b*x)","y~a*(1+b*(exp(-c*x)))^-1","y~a*(1-b*(exp(-c*x)))^3","y~a*(1-b*(exp(-c*x)))",
             "y~a*exp(-b*exp(-c*x)","y~(a*x^b)*exp(-c*x)","y ~ a + b * (1 - exp(-c * x))","y~(a/(1+exp(2-4*c*(x-e))))+(b/(1+exp(2-4*d*(x-e))))")
    r2=models[model]
    r=list(r2,r1);names(r)=c("Model","Parameters")
    return(r)
    
}
