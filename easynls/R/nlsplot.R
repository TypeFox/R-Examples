nlsplot <-
function(data, model=1, start=c(a=1,b=1,c=1,d=1,e=1), xlab="Response Variable", ylab="Explanatory Variable", position=1)
{
    
    res=nlsfit(data=data, model=model, start=c(a=start[1],b=start[2],c=start[3],d=start[4],e=start[5]))
    
    res=res[[2]]
    means = function(data) {
        t = as.factor(data[, 1])
        d = data.frame(t, data[, -1])
        s = split(data.frame(d[, -1]), d$t)
        r = lapply(s, colMeans, na.rm = TRUE)
        r = lapply(r, round, 2)
        rr = t(data.frame(r))
        rr = data.frame(rr)
        rownames(rr) = NULL
        treat = levels(t)
        rr = data.frame(treat, rr)
        colnames(rr) = colnames(data)
        return(rr)
    }
    
    t = list("top", "bottomright", "bottom", "bottomleft", 
             "left", "topleft", "topright", "right", "center")
    p = t[[position]]
    
    datao=means(data)
    
    plot1=function(data)
    {
        f1=function(i)res[1,i]+res[2,i]*se
        f2=function(i)res[1,i]+res[2,i]*se+res[3,i]*se^2
        f3=function(i)res[1,i] + res[2,i] * (se - res[3,i]) * (se <= res[3,i])
        f4=function(i)(res[1,i] + res[2,i] * se + res[3,i] * I(se^2)) * (se <= -0.5 * res[2,i]/res[3,i]) + 
            (res[1,i] + I(-res[2,i]^2/(4 * res[3,i]))) * (se > -0.5 * res[2,i]/res[3,i])
        f5=function(i)ifelse(se>=res[4,i],(res[1,i]-res[3,i]*res[4,i])+(res[2,i] +res[3,i])*se, res[1,i]+res[2,i] *se)
        f6=function(i)res[1,i]*exp(res[2,i]*se)
        f7=function(i)res[1,i]*(1+res[2,i]*(exp(-res[3,i]*se)))^-1
        f8=function(i)res[1,i]*(1-res[2,i]*(exp(-res[3,i]*se)))^3
        f9=function(i)res[1,i]*(1-res[2,i]*(exp(-res[3,i]*se)))
        f10=function(i)res[1,i]*exp(-res[2,i]*exp(-res[3,i]*se))
        f11=function(i)(res[1,i]*se^res[2,i])*exp(-res[3,i]*se)
        f12=function(i)res[1,i] + res[2,i] * (1 - exp(-res[3,i] * se))
        f13=function(i)(res[1,i]/(1+exp(2-4*res[3,i]*(se-res[5,i]))))+(res[2,i]/(1+exp(2-4*res[4,i]*(se-res[5,i]))))
        
        mod=list(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13)
        se=seq(min(data[,1]),max(data[,1]),by=0.01)
        i=1:(ncol(data)-1)
        pred=lapply(i,mod[[model]])
        u=unlist(pred)
        mat=matrix(u, ncol=length(pred))
        pred=as.data.frame(mat)
        dd=data.frame(se,pred)
        n=names(data)
        names=n[-1]
        matplot(dd[,1], dd[,-1], type = "l", bty = "n", xlab = xlab, 
                ylab = ylab, col = c(1:(ncol(data)-1)), lty = c(1:(ncol(data)-1)))
        legend(p, names, bty = "n", col = c(1:(ncol(data)-1)), 
               lty = c(1:(ncol(data)-1)))
    }
    
    
    
    plot2=function(data)
    {
        ff1=function(se)res[1,i]+res[2,i]*se
        ff2=function(se)res[1,i]+res[2,i]*se+res[3,i]*se^2
        ff3=function(se)res[1,i] + res[2,i] * (se - res[3,i]) * (se <= res[3,i])
        ff4=function(se)(res[1,i] + res[2,i] * se + res[3,i] * I(se^2)) * (se <= -0.5 * res[2,i]/res[3,i]) + 
            (res[1,i] + I(-res[2,i]^2/(4 * res[3,i]))) * (se > -0.5 * res[2,i]/res[3,i])
        ff5=function(se)ifelse(se>=res[4,i],(res[1,i]-res[3,i]*res[4,i])+(res[2,i] +res[3,i])*se, res[1,i]+res[2,i] *se)
        ff6=function(se)res[1,i]*exp(res[2,i]*se)
        ff7=function(se)res[1,i]*(1+res[2,i]*(exp(-res[3,i]*se)))^-1
        ff8=function(se)res[1,i]*(1-res[2,i]*(exp(-res[3,i]*se)))^3
        ff9=function(se)res[1,i]*(1-res[2,i]*(exp(-res[3,i]*se)))
        ff10=function(se)res[1,i]*exp(-res[2,i]*exp(-res[3,i]*se))
        ff11=function(se)(res[1,i]*se^res[2,i])*exp(-res[3,i]*se)
        ff12=function(se)res[1,i] + res[2,i] * (1 - exp(-res[3,i] * se))
        ff13=function(se)  (res[1,i]/(1+exp(2-4*res[3,i]*(se-res[5,i]))))+(res[2,i]/(1+exp(2-4*res[4,i]*(se-res[5,i]))))
        mod2=list(ff1,ff2,ff3,ff4,ff5,ff6,ff7,ff8,ff9,ff10,ff11,ff12,ff13)
        i=1:(ncol(data)-1)
        minx = min(data[, 1]) - sd(data[, 1])/2
        maxx = max(data[, 1]) + sd(data[, 1])/2
        miny = min(data[, 2]) - sd(data[, 2])/2
        maxy = max(data[, 2]) + sd(data[, 2])/2
        
        c1=res[1,];c2=res[2,];c3=res[3,];c4=res[4,];c5=res[5,]
        
        sin1 = ifelse(c1 > 0, "+", "")
        sin2 = ifelse(c2 > 0, "+", "")
        sin3 = ifelse(c3 > 0, "+", "")
        sin4 = ifelse(c4 > 0, "+", "")
        sin5 = ifelse(c5 > 0, "+", "")
        
        r2=res["r-squared",]
        
        e1 = substitute(y == c1 * sin1 * c2 * x * "  " * R^2 * " = " * r2, list(c1 = c1, c2 = c2, 
                                                                                r2 = r2, sin1 = sin1))
        
        e2 = substitute(y == c1 * sin2 * c2 * x * sin3 * c3 * 
                            x^2 * "  " * R^2 * " = " * r2, list(c1 = c1, c2 = c2, 
                                                                c3 = c3, r2 = r2, sin2 = sin2, sin3 = sin3))
        e3 = substitute(y == c1 * sin2 * c2 * (x -c3) * 
                            "  " * R^2 *" = " * r2, list(c1 = c1, c2 = c2, r2 = r2, c3 = c3, 
                                                         sin2 = sin2))
        e4 = substitute(y == c1 * sin2 * c2 * x * sin3 * c3 * 
                            x^2 *  "  " * R^2 * " = " * r2, list(c1 = c1, c2 = c2, c3 = c3,sin2 =sin2, sin3=sin3, r2=r2)) 

cn1=c1-(c3*c4)
cn2=c2+c3

        sin6 = ifelse(cn2 > 0, "+", "")                                             
        e5 = substitute(atop(y1==c1*sin2*c2*x, 
     y2==cn1*sin6*cn2*x )*"    "* R^2 * " = " * r2, list(c1 = c1, c2 = c2, c3 = c3,c4 = c4,cn1=cn1,cn2=cn2, r2 = r2, sin2 = sin2, sin3 = sin3,sin4 = sin4, sin6 = sin6))
        
        e6 = substitute(y==c1*e^{c2*x}*"  " * R^2 * " = " * r2, list(c1 = c1, c2 = c2, 
                                                                     r2 = r2)) 
        e7 = substitute(y==c1*(1*sin2*c2*e^{-c3*x})^-1*"  " * R^2 * " = " * r2, list(c1 = c1, c2 = c2, 
                                                                                     sin2 = sin2, c3 = c3, r2 = r2)) 
        e8 = substitute(y==c1*(1-c2*e^{-c3*x})^3*"  " * R^2 * " = " * r2, list(c1 = c1, c2 = c2, 
                                                                               c3 = c3, r2 = r2)) 
        e9 = substitute(y==c1*(1-c2*e^{-c3*x})*"  " * R^2 * " = " * r2, list(c1 = c1, c2 = c2, 
                                                                             c3 = c3, r2 = r2)) 
        e10 = substitute(y==c1*e^(-c2*e^{-c3*x})*"  " * R^2 * " = " * r2, list(c1 = c1, c2 = c2, 
                                                                               c3 = c3, r2 = r2)) 
        e11 = substitute(y==c1*x^c2*e^{-c3*x}*"  " * R^2 * " = " * r2, list(c1 = c1, c2 = c2, 
                                                                            c3 = c3, r2 = r2))
        e12 = substitute(y==c1 + c2 * (1 - e^{-c3 * x})*"  " * R^2 * " = " * r2, list(c1 = c1, c2 = c2, 
                                                                                      c3 = c3, r2 = r2))
        cc3=4*c3
        cc4=4*c4
        e13 = substitute(y==frac(c1,1+e^(2-cc3*(x-c5)))+frac(c2,1+e^(2-cc4*(x-c5)))*"  " * R^2 * " = " * r2, list(c1 = c1, c2 = c2, 
                                                                                                                  cc3 = cc3, r2 = r2, cc4=cc4, c5=c5))
        
        # demo(plotmath)
        
        ee=list(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13)
        
        eee=ee[[model]]
        
        plot(datao[,2]~as.numeric(as.character(datao[,1])), data=datao, xlab = xlab, 
             ylab = ylab, col="dark red", bty="n",xlim = c(minx, 
                                                           maxx), ylim = c(miny, maxy))
        plot(mod2[[model]], min(data[,1]), max(data[,1]), add = TRUE, 
             col = "dark blue", lty = 2)
        legend(p,legend=eee,bty = "n", cex=0.9)
    }
    ppp=ifelse(ncol(data)==2,2,1)
    lll=list(plot1,plot2)
    lll[[ppp]](data)
    
}
