mywhiskers<-function(x, y, nclass=10, limits=NA,
                     add=FALSE, se=TRUE, main="", xlab="x", ylab="y", 
                     ylim=NA, lwd=1, highlight="red") {
            away<-is.na(x+y)
            x<-x[!away]
            y<-y[!away]
			means<-ses<-NA
			nclass<-nclass+1
			if (is.na(limits[1])) { 
			while (sum(is.na(means+ses))&nclass>1) { # decrease the number of classes if there are too few observations
			   nclass<-nclass-1
               limits<-seq(min(x),max(x)+1e-10,length=nclass+1) 
 
               means<-sapply(1:nclass,function(i) mean(y[x>=limits[i]&x<limits[i+1]]))
               if (se) { 
                  ses<-sapply(1:nclass,function(i) sd(y[x>=limits[i]&x<limits[i+1]])/
                                                      sqrt(sum(x>=limits[i]&x<limits[i+1])))
                  } else {
                  ses<-sapply(1:nclass,function(i) sd(y[x>=limits[i]&x<limits[i+1]]))
                  }
            }
		    } else { 
				nclass=length(limits)-1
				means<-sapply(1:nclass,function(i) mean(y[x>=limits[i]&x<limits[i+1]]))
				if (se) { 
					ses<-sapply(1:nclass,function(i) sd(y[x>=limits[i]&x<limits[i+1]])/
										sqrt(sum(x>=limits[i]&x<limits[i+1])))
				} else {
					ses<-sapply(1:nclass,function(i) sd(y[x>=limits[i]&x<limits[i+1]]))
				}
			}
            lb<-means-1.96*ses
            ub<-means+1.96*ses
            xclass<-1/2*(limits[-1]+limits[-nclass-1])
            if (add) {
               points(xclass,means) 
               } else {
               if(is.na(ylim[1])) ylim<-c(min(lb),max(ub))
               plot(xclass,means,ylim=ylim,main=main,xlab=xlab,ylab=ylab,xlim=range(x))
               }
            color<-rep("black",nclass)
            if (sum(ub*lb>0)) color[ub*lb>0]<-highlight
            sapply(1:nclass,function(i) lines(xclass[c(i,i)],c(lb[i],ub[i]),lwd=lwd,col=color[i],lend="butt"))
            list(x=xclass,m=means,s=ses,lb,ub)
            }
            

# This is for plotting the response against age and connecting the observations of the same subject
# in a longitudinal study.
# in a longitudinal study.
linesplot<-function(x, y, group,
		xlab="x", ylab="y", main="", cex=0.5, pch=19, 
		col=1, col.lin=1, lw=FALSE, ylim=NULL, xlim=NULL,
		add=FALSE, lty="solid",lwd=1) {
	group<-group[order(x)]
	y<-y[order(x)]
	
	col<-rep(col,length(x))
	if (length(col.lin)==1) col.lin=rep(col.lin,length(x))
	col.lin<-col.lin[order(x)]
	
	if (length(lty)==1) lty=rep(lty,length(x))
	lty<-lty[order(x)]
	
	if (length(lwd)==1) lwd=rep(lwd,length(x))
	lwd<-lwd[order(x)]
	
	x<-x[order(x)]
	if (!add) {
		plot(x,y,type="n",xlab=xlab,ylab=ylab,main=main,cex=cex,pch=pch,col=col,ylim=ylim,xlim=xlim)
	}
	apu<-unique(group)
	for (i in 1:length(apu)) {
		print(col.lin[group==apu[i]])
		lines(x[group==apu[i]],y[group==apu[i]],col=col.lin[group==apu[i]],lty=lty[group==apu[i]],
				lwd=lwd[group==apu[i]])
	}
	points(x,y,col=col,cex=cex)
	if (lw) {
		apu<-unique(col)
		for (i in 1:length(apu))
			lines(lowess(x[col==apu[i]],y[col==apu[i]],f=0.5,iter=5),col=apu[i],lwd=3)
	}
}
# this function plots circles of radius r at locations specified by x and y
circle<-function(x, y, r, col="black", lty="solid", lwd=1, grayfill=FALSE) {
        xapu<-sin(seq(0,pi,length=50)-pi/2)
        for (i in 1:length(x)) {
            if (grayfill) {
               xv1<-x[i]+xapu*r[i]/2
               xv2<-x[i]-xapu*r[i]/2
               yv1<-sqrt(pmax(0,(r[i]/2)^2-(xv1-x[i])^2))+y[i]
               yv2<--sqrt(pmax(0,(r[i]/2)^2-(xv2-x[i])^2))+y[i]
               yv<-c(yv1,yv2)
               xv<-c(xv1,xv2)
               lines(xv,yv,col="gray",lwd=5*r[i])               
               }
            xv1<-x[i]+xapu*r[i]
            xv2<-x[i]-xapu*r[i]
            yv1<-sqrt(pmax(0,r[i]^2-(xv1-x[i])^2))+y[i]
            yv2<--sqrt(pmax(0,r[i]^2-(xv2-x[i])^2))+y[i]
            yv<-c(yv1,yv2)
            xv<-c(xv1,xv2)
            lines(xv,yv,col=col,lty=lty,lwd=lwd)
            }
        }
        
# Newton-Raphson algorithm for n multidimensional functions using numeric differentiation
NRnum<-function(init, fnlist, crit=6,...) {
        par<-init
#        sapply(par,function(x) cat(x," "))
#        cat("\n")
        value<-sapply(fnlist,function(f) f(par,...))
        j<-1
        while (round(sum(abs(value)),crit)!=0&j<100) {
#              cat(j,"\n")
              grad<-t(sapply(fnlist,function(f) attributes(numericDeriv(quote(f(par)),c("par")))$gradient))
              deltax<-solve(grad,-value)
              par<-par+deltax
#              sapply(par,function(x) cat(x," "))
#              cat("\n")
              value<-sapply(fnlist,function(f) f(par,...))
              j<-j+1
              }
        if (j<100) {
           list(par=par,value=value)
           } else {
           list(par=NA,value=NA)
           }
        }
        
# An up-and-down algorithm
# Solves fn(x)=0 for x between l (lower bound) and u (upper bound).
updown<-function(l, u, fn, crit=6) {
        fnu<-fn(u)
        fnl<-fn(l)
        if (fnl*fnu>0) return(NA)
        value<-fn((u+l)/2)
        while (round(value,crit)!=0&(u-l)>10^(-crit)) {
              if (fnu*value>0) {
                 u<-(u+l)/2 
                 fnu<-value 
                 } else {
                 l<-(u+l)/2
                 fnl<-value
                 }
              value<-fn((u+l)/2) 
#              cat(l,fnl,u,fnu,"\n")
              }
        nollakohta<-(l+u)/2
        nollakohta
        }


# Solve unidimensional equation using Newton-Rapson algorithm
# init=Initial guess, 
# fn=function, 
# gr=gradient  
# crit= convergence criteria 
# range= range of the possible values  
NR<-function(init, fn, gr, crit=6, range=c(-Inf,Inf)) {
    par<-init
    value<-fn(par)
    j<-1
    while (round(value,crit)!=0&j<100000) {
          grad<-gr(par)
          a<-value-grad*par
          par<--a/grad
          par<-min(max(par,range[1]),range[2])
          value<-fn(par)
          cat(par,value,"\n")
          j<-j+1
          }
    if (j<100000) {
       list(par=par,value=value)
       } else {
       list(par=NA,value=NA)
       }
    }

	HDnaslund<-function(d,a,b,bh=1.3) {
		d^2/(a+b*d)^2+bh
	}
	
	startHDnaslund<-function(d,h,bh=1.3) {
		start<-coef(lm(I(d/sqrt(h-bh))~d))
		names(start)<-c("a","b")
		start
	}	

	
	
# Michailoff 1943. H=a exp{-bD^{-1}}
# Lundqvist-Korf with c=1
	HDmichailoff<-function(d,a,b,bh=1.3) {
		a*exp(-b*d^(-1))+bh
	}
	
	
	startHDmichailoff<-function(d,h,bh=1.3) {
		a<-1.3*max(h-bh)
		mod<-lm(log((h-bh))~I(1/d))
		start<-c(exp(coef(mod)[1]),-coef(mod)[2])
		names(start)<-c("a","b")
		start
	}

	
# Curtis says it can be found in Prodan(book, 1965) or Strand (article, 1951).
	
	HDcurtis<-function(d,a,b,bh=1.3) {
		a*(d/(1+d))^b+bh
	}
	
	startHDcurtis<-function(d,h,bh=1.3) {
		start<-coef(lm(log(h)~log(d/(1+d))))
		start[1]<-exp(start[1])
		names(start)<-c("a","b")
		start
	}	
	

         
         
         
# Meyer: H=a(1-exp{-bD}) Meyer 1940
HDmeyer<-function(d,a,b,bh=1.3) {
         a*(1-exp(-b*d))+bh
         }
         
startHDmeyer<-function(d,h,bh=1.3) {
         a<-1.3*max(h-bh)
	 b<--coef(lm(log(1-(h-bh)/a)~d-1))
         start<-c(a,b) 
         names(start)<-c("a","b")
         start
         }  
         
         
# Power function
# Stoffels 1950
		 HDpower<-function(d,a,b,bh=1.3) {
			 a*d^b+bh
		 }
		 
		 startHDpower<-function(d,h,bh=1.3) {
			 start<-coef(lm(log(h-bh)~log(d)))  # includes starting values estimates of a and b
			 start[1]<-exp(start[1])
			 names(start)<-c("a","b")
			 start
		 }
         

HDnaslund2<-function(d,a,b,bh=1.3) {
           d^2/(a+exp(b)*d)^2+bh
           }
           
startHDnaslund2<-function(d,h,bh=1.3) {
         start<-coef(lm(I(d/sqrt(h-bh))~d))
         start[2]<-log(start[2])
         names(start)<-c("a","b")
         start
         }         

HDnaslund3<-function (d, a, b, bh = 1.3) {
			 d^2/(exp(a) + b * d)^2 + bh
		 }
		 
startHDnaslund3<-function (d, h, bh = 1.3) {
         start <- coef(lm(I(d/sqrt(h - bh)) ~ d))
         start[1]<-log(max(start[1],0.1))
         names(start) <- c("a", "b")
         start
		 }		 
		  
		 
HDnaslund4<-function (d, a, b, bh = 1.3) {
			 d^2/(exp(a) +exp( b) * d)^2 + bh
		 }
		 
startHDnaslund4<-function (d, h, bh = 1.3) {
			 start <- coef(lm(I(d/sqrt(h - bh)) ~ d))
			 start[1]<-log(max(start[1],0.1))
			 start[2]<-log(max(start[2],0.1))
			 names(start) <- c("a", "b")
			 start
		 }
		          
# Bates and Watts aD/(b+d) (Calama and Montero 2004)
HDmicment<-function(d,a,b,bh=1.3) {
           a*d/(b+d)+bh
           }
           
startHDmicment<-function(d,h,bh=1.3) {
           tmp<-coef(lm(I(d/(h-bh))~d))
           start<-c(1/tmp[2],tmp[1]/tmp[2])
		   names(start)<-c("a","b")
           start
           }		   
		   
           

# Bates and Watts aD/(b+d) (Calama and Montero 2004)
HDmicment2<-function(d,a,b,bh=1.3) {
           d/(a+b*d)+bh
           }
           
startHDmicment2<-function(d,h,bh=1.3) {
         start<-coef(lm(I(d/(h-bh))~d))
         names(start)<-c("a","b")
         start
         }
         
         
           
HDwykoff<-function(d,a,b,bh=1.3) {
          exp(a+b/(d+1))
          }
          
startHDwykoff<-function(d,h,bh=1.3) {
         mod<-lm(log((h-bh))~I(1/d))
         start<-c(coef(mod)[1],coef(mod)[2])
         names(start)<-c("a","b")
         start
         }


         
# 3- parameter models
HDprodan<-function(d,a,b,c,bh=1.3) {
           d^2/(a+b*d+c*d^2)+bh
           }
           
startHDprodan<-function(d,h,bh=1.3) {
         start<-coef(lm(I(d^2/(h-bh))~d+I(d^2)))
         names(start)<-c("a","b","c")
         start
         }
         
HDlogistic<-function(d,a,b,c,bh=1.3) {
           a/(1+b*exp(-c*d))+bh
           }
           
startHDlogistic<-function(d,h,bh=1.3) {
         a<-1.3*max(h-bh)
         start<-c(a,coef(lm(log((a-h+bh)/(h-bh))~d)))
         start[2]<-exp(start[2])
         start[3]<--start[3]
         names(start)<-c("a","b","c")
         start
         }
         
HDrichards<-function(d,a,b,c,bh=1.3) {
           a*(1-exp(-b*d))^c+bh
           }
           
startHDrichards<-function(d,h,bh=1.3,b=0.04) {
         start<-coef(lm(log(h-bh)~log(1-exp(-b*d))))
         start<-c(exp(start[1]),b,start[2])
         names(start)<-c("a","b","c")
         start
         }

HDweibull<-function(d,a,b,c,bh=1.3) {
           a*(1-exp(-b*d^c))+bh
           }
           
startHDweibull<-function(d,h,bh=1.3) {
         a<-1.3*max(h-bh)
         start<-c(a,coef(lm(log(-log(1-(h-bh)/a))~log(d))))
         start[2]<-exp(start[2])
         names(start)<-c("a","b","c")
         start
         }   


HDgomperz<-function(d,a,b,c,bh=1.3) {
           a*exp(-b*exp(-c*d))+bh
           }
           
startHDgomperz<-function(d,h,bh=1.3) {
         a<-1.3*max(h-bh)
         start<-c(a,coef(lm(log(-log((h-bh)/a))~d)))
         start[2]<-exp(start[2])
         start[3]<--start[3]
         names(start)<-c("a","b","c")
         start
         }  
         
HDsibbesen<-function(d,a,b,c,bh=1.3) {
           a*d^(b*d^(-c))+bh
           }
           
startHDsibbesen<-function(d,h,bh=1.3,a=0.5) {
         start<-c(a,coef(lm(log(log((h-bh)/a)/log(d))~log(d))))
         start[2]<-exp(start[2])
         start[3]<--start[3]
         names(start)<-c("a","b","c")
         start
         }  
         

# Korf: H=a exp{-bD^{-c}}
HDkorf<-function(d,a,b,c,bh=1.3) {
         val<-a*exp(-b*d^(-c))+bh
    #     print(a,b,c)
#         points(d,val,col="green")
         val
         }
         
startHDkorf<-function(d,h,bh=1.3) {
         a<-1.3*max(h-bh)
         y<-log(-log((h-bh)/a))
         startmod<-lm(y[is.finite(y)]~log(d)[is.finite(y)],na.action=na.omit)
         start<-c(a,exp(coef(startmod)[1]),-coef(startmod)[2])
         names(start)<-c("a","b","c")
         start
         }


# Korf: H=a exp{-bD^{-c}}
HDratkowsky<-function(d,a,b,c,bh=1.3) {
         val<-a*exp(-b/(d+c))+bh
    #     print(a,b,c)
#         points(d,val,col="green")
         val
         }
         
startHDratkowsky<-function(d,h,bh=1.3,c=5) {
         start<-c(coef(lm(log(h-bh)~I(1/(d+c)))),c)
         start[1]<-exp(start[1])
         start[2]<--start[2]
         names(start)<-c("a","b","c")
         start
         }
         
## Korf: H=a exp{-bD^{-c}}
#HDhossfeldIV<-function(d,a,b,c,bh=1.3) {
#         val<-a/(1+1/b*d^(-c))+bh
#         val
#         }
#         
#startHDhossfeldIV<-function(d,h,bh=1.3,c=5) {
#         a<-1.3*max(h-bh)
#         start<-c(a,coef(lm(log(a/(h-bh)-1)~log(d))))
#         start[2]<-1/exp(start[2])
#         start[3]<--start[3]
#         names(start)<-c("a","b","c")
#         start
#         }         

HDhossfeldIV<-function(d,a,b,c,bh=1.3) {
         val<-a/(1+1/(b*d^c))+bh
         val
         }
         
startHDhossfeldIV<-function(d,h,bh=1.3,c=5) {
         a<-1.3*max(h-bh)
         start<-c(a,coef(lm(log(1/(a/h-1))~log(d))))
         start[2]<-exp(start[2])
         names(start)<-c("a","b","c")
         start
         }
         
         
HDmehtatalo<-function(d,a,b,c,bh=1.3) {
           val<-d^c/(a+exp(b)*d)^c+bh
           points(d,val,col="red")
           val
           }

startHDmehtatalo<-function(d,h,bh=1.3) {
         start<-c(coef(lm(I(d/sqrt(h-bh))~d)),2)
         start[2]<-log(start[2])
         names(start)<-c("a","b","c")
         start
         }



         
#startHDkorf<-function(d,h,bh=1.3) {
#         a<-1.3*max(h-bh)
#         mod<-lm(log((h-bh)/a)~I(1/d)-1)
#         start<-c(a,-coef(mod))
#         names(start)<-c("a","b")
#         start
#         }


         
         

         

           

         

                 
#startHDkorf3<-function(d,h,bh=1.3) {
#         a<-max(h-bh)
#         b<--coef(lm(log((h-bh)/a)~I(1/d)-1))
#         start<-c(a,b,1)
#         names(start)<-c("a","b","c")
#         start
#         }
         
# Lineaariset
# h=a+b*log(D)
# a*D/(D+1)+b*D
# 2. asteen polynomi
# spliniregressio

# Fits one of the commonly used models to H-D data and returns an
# object of class hdmod. 
fithd<-function(d, h, plot=c(), modelName="naslund",
		        nranp=2, random=NA, varf=0, na.omit=TRUE, 
				start=NA, bh=1.3, control = list(), SubModels=NA, vfstart=0) {

            if (!is.na(list(random))) nranp<-1

			varf<-as.numeric(varf)
			
            if (na.omit) {
               sel<-is.na(h*d)
               d<-d[!sel]
               plot<-plot[!sel]
               h<-h[!sel]
               }
			   
#			if (!is.na(SubModels)[1]) {
               dmean<-tapply(d,plot,mean)
			   dmean<-dmean[match(plot,names(dmean))]
			   
			   dsd<-tapply(d,plot,sd)
			   dsd<-dsd[match(plot,names(dsd))]
			   
			   dstd<-(d-dmean)/dsd
#		       }
	
	        if (varf==1)  w<-d  else if (varf==2) w<-pmax(1,dstd+3)

# Linear model fitting              
            if (class(modelName)=="formula") { 
			   if (!is.na(SubModels)[1]) warning("Argument 'SubModels' not implemented for linear model fit")
               modEq<-modelName
               if (!is.na(random)) {
                  # If a formula is provided in ranp, use it.
                  ranEq=as.formula(paste(deparse(random),"|plot"))
                  } else if (nranp>0) {
                  # Otherwise take the given number of trems of model from the beginning. 
                  # Constant is always assumed to be the first if model just has it.
                  ranTerms<-paste("+",attributes(terms(modelName))$term.labels,sep="")
                  if (attributes(terms(modelName))$intercept) {
                     ranTerms<-c("1",ranTerms)
                     ranp<-1:nranp
                     } else {
                     ranTerms<-c("-1",ranTerms)
                     ranp<-1:(nranp+1)
                     }
                  ranEq<-"~"
                  for (j in ranp) {
                      ranEq<-paste(ranEq,ranTerms[j],sep="")
                      } 
                  ranEq<-as.formula(paste(ranEq,"|plot",sep=""))
                  } 
                 
               
               if (!nranp) {
                  if (varf) {
                     mod<-gls(modEq, 
                         weights=varPower(vfstart,~w),
						 control=control)         
                     } else { 
                     mod<-gls(modEq,
							 control=control)
                     }
                  } else {
                  if (varf) {
                     mod<-lme(fixed=modEq, 
                         random=ranEq,
                         weights=varPower(vfstart,~w),
						 control=control)         
                     } else { 
                     mod<-lme(fixed=modEq, 
                         random=ranEq,
						 control=control)
                     }
                  }
               class(mod)<-c("hdmod",class(mod))
               attributes(mod)$model<-modelName
               attributes(mod)$d<-d
               attributes(mod)$h<-h
               attributes(mod)$plot<-plot
               mod$call[[2]]<-modEq
               if (nranp) mod$call[[3]]<-ranEq
               
               } else { 
               
# Nonlinear model fitting
               if (is.na(start[1])) {
				  start<-eval(call(paste("startHD",modelName,sep=""),d=d,h=h,bh=bh))
				  if (!is.na(SubModels)[1]) {
					 start2<-c()
					 for (i in 1:length(SubModels)) {
					     start2<-c(start2,start[i])
					     start2<-c(start2,0[as.character(SubModels[i])!="1"])
					 }
					 start<-start2
			      }
			  }
             
       #        # The # of random parameters should not exceed the # of parameters
           #    ranp<-min(ranp,length(formals(paste("HD",model,sep="")))-2)
               
               if (!is.na(list(random))) {
                  ranEq<-random
                  } else if (nranp) {
                  ranEq<-lapply(paste(names(formals(paste("HD",modelName,sep="")))[1+1:nranp],"~1",sep=""),as.formula)
                  }

   #            if (nranp==1) {
   #               ranEq<-list(a~1)
   #               } else if (nranp==2) {
   #               ranEq<-list(a~1,b~1)
   #               } else if (nranp==3) {
   #               ranEq<-list(a~1,b~1,c~1)
   #               }
               
               if (length(formals(paste("HD",modelName,sep="")))==5) {
            #     if (modelName=="korf3") {
                    modEq<-as.formula(paste("h~HD",modelName,"(d,a,b,c,bh=",bh,")",sep=""))
					if (is.na(SubModels)[1]) {
                       fixEq<-list(as.formula("a~1"),as.formula("b~1"),as.formula("c~1"))
				       } else if (length(SubModels==3)) {
					   fixEq<-lapply(paste(c("a~","b~","c~"),SubModels),as.formula)
                     #  fixEq<-list(as.formula("a~1"),as.formula("b~1"),as.formula("c~1"))   
					   } else {
						stop("length of SubModels is not 3")
					   }
                    } else {
                    modEq<-as.formula(paste("h~HD",modelName,"(d,a,b,bh=",bh,")",sep=""))
					if (is.na(SubModels[1])) {
						fixEq<-list(as.formula("a~1"),as.formula("b~1"))
					} else if (length(SubModels==2)) {
						fixEq<-lapply(paste(c("a~","b~"),SubModels),as.formula)
						#  fixEq<-list(as.formula("a~1"),as.formula("b~1"),as.formula("c~1"))   
					} else {
						stop("length of SubModels is not 2")
					}
                    
                    }
                 

                
                 if (!nranp) {
                    if (!varf) {
                       mod<-gnls(modEq,
							  params=fixEq,
                              start=start,
							  control=control
                              )
                    } else {
                       mod<-gnls(modEq,
							  params=fixEq,
                              start=start,
                              weights=varPower(vfstart,~w),
							  control=control
                              )                 
                    }
                 } else {              
                    if (!varf) {
                       mod<-nlme(modEq,
                              fixed=fixEq,
                              random=ranEq,
                              groups=~plot,
                              start=start,
							  control=control)
                       } else {
                       mod<-nlme(modEq,
                              fixed=fixEq,
                              random=ranEq,
                              groups=~plot,
                              start=start,
                              weights=varPower(vfstart,~w),
							  control=control)
                       }
                    }
                 class(mod)<-c("hdmod",class(mod))
                 attributes(mod)$model<-modelName
                 attributes(mod)$d<-d
                 attributes(mod)$h<-h
                 attributes(mod)$plot<-plot
                 mod$call[[2]]<-modEq
                 mod$call[[3]]<-fixEq
                 if (nranp) mod$call[[4]]<-ranEq
                 }
                 mod
                 }



# Plot a height-diameter model
plot.hdmod<-function(x,col.point="blue",highlight="red",standd=TRUE,cex=1,corD=FALSE,ask=TRUE,...) {
            model<-x
            res<-resid(model,type="p")
            d<-attributes(model)$d
            plot<-attributes(model)$plot
            plots<-unique(plot)
            dstd<-d

            # compute diameter that is standardized for each plot
            if (standd) {
               for (i in 1:length(plots)) {
                   # thisplot will include observations only from i:th plot
                   dthis<-d[plot==plots[i]]
                   dstd[plot==plots[i]]<-(dthis-mean(dthis))/sd(dthis)
                   }
               plot(dstd,
                    res,
                    col=col.point,cex=cex,
                    main=paste(attributes(model)$model,"s.e.=",round(sd(resid(model)),3)),
                    xlab="Relative diameter",
                    ylab="Standardized residual")
            } else {
               plot(dstd,
                    res,
                    col=col.point,cex=cex,
                    main=paste(attributes(model)$model,"s.e.=",round(sd(resid(model)),3)),
                    xlab="Diameter, cm",
                    ylab="Standardized residual")
            }
            mywhiskers(dstd,res,add=TRUE,highlight="black",se=FALSE)
            mywhiskers(dstd,res,add=TRUE,highlight=highlight,lwd=3)
            abline(0,0)
			if (inherits(model,"nlme")&corD) {
				re<-ranef(model)
				re$plot<-as.numeric(rownames(ranef(model)))
				dmean<-tapply(d,plot,mean)
				re$dmean<-dmean[match(re$plot,names(dmean))]
				nranp<-dim(ranef(model))[2]
			    for (i in 1:nranp) {
			        devAskNewPage(ask)
			        plot(log(re$d),re[,i],xlab="log(dmean)",ylab=paste(i,"th random effect",sep=""))
					abline(lm(re[,i]~log(re$d)))
				    }
				}
            }
            
# Plot a height-diameter model
qqplotHD<-function(model,startnew=TRUE) {
            res<-resid(model,type="p")
 
            if (sum(class(model)=="lme")) {
               rem<-ranef(model)
               nplots<-dim(rem)[2]+1
               if (startnew) { 
                  par(mfcol=c(2,ceiling(nplots/2)))

                  }                  
               qqnorm(res,main="N-QQ-plot of standardized residuals")
               qqline(res)
               sapply(1:dim(rem)[2],function(x) {qqnorm(rem[,x],main=paste("QQ-plot of",colnames(rem)[x])); qqline(rem[,x])})
               } else {
               qqnorm(res,main="N-QQ-plot of standardized residuals")
               qqline(res)
               }
            }

		
			
ImputeHeights<-function(d, h, plot=c(),modelName="naslund",nranp=2, varf=TRUE,  
                        addResidual=FALSE, makeplot=TRUE, level=1,
						start=NA,bh=1.3,control=list(),random=NA) {

               data<-data.frame(d=d,h=h,plot=plot)
               data1<-data[!is.na(data$h*data$d),]

               mod<-fithd(data1$d,data1$h,data1$plot,
                          modelName=modelName,nranp=nranp,random=random,
						  varf=varf,start=start,bh=bh,control=control)                                       

               hplots<-unique(data$plot[!is.na(data$h*data$d)])
               type<-rep(2,dim(data)[1])
               for (i in 1:length(hplots)) {
                   type[data$plot==hplots[i]]<-1
                   } 
               nohplots<-unique(data$plot[type==2])
               
               if (makeplot) plot(mod)  

               # To initialize, compute the predictions using the fixed part only. 
               # This will remain as the final predioction 
               # a) if nranp=0 (no randiom effects) and addResidual==0 
               # b) if level=0 (fixed part prediction is requested) and addResidual==0. 
               
               predall<-predict(mod,newdata=data, level=0)
	                                    
               # If a higher level prediction is requested, then the realized random effects are  used.
               # If Some plots did not have measured heights for the preiction of random effects,
               # Then random effects are sampled from among the predictions. (Not using MVN(0,D)!)
               if (level>0 && nranp>0) {            
                  predall[type==1]<-predict(mod,newdata=data[type==1,])
                 
                  D<-as.matrix(mod$modelStruct$reStruct[[1]])*mod$sigma^2
                  if (length(nohplots)) {
                     for (i in 1:length(nohplots)) {
                  #       thiscoef<-fixef(mod)+mvrnorm(1,mu=rep(0,dim(D)[1]),Sigma=D)
                         thiscoef<-coef(mod)[sample(1:mod$dims$ngrps[1],1),]
                         if (length(formals(paste("HD",modelName,sep="")))==5) {
                            predall[data$plot==nohplots[i]]<-
                                eval(call(paste("HD",modelName,sep=""),
                                     d=data$d[data$plot==nohplots[i]],a=thiscoef[[1]],b=thiscoef[[2]],c=thiscoef[[3]],bh=bh))
                            } else {
                            predall[data$plot==nohplots[i]]<-
                                eval(call(paste("HD",modelName,sep=""),
                                     d=data$d[data$plot==nohplots[i]],a=thiscoef[[1]],b=thiscoef[[2]],bh=bh))
                            }         
         #                predall[data$plot==nohplots[i]] <- HDnaslund(d=data$d[data$plot==nohplots[i]],a=thiscoef[[1]],b=thiscoef[[2]])
                         }
                     }
                  }
          
               # If level = 0 and AddResidual=True, then random effects are randomly assigned for each plot
               # before adding the tree level residual
               if (level==0 & nranp>0 & addResidual) {
                  plots<-unique(data$plot)
                  for (i in 1:length(plots)) {
               #       thiscoef<-fixef(mod)+mvrnorm(1,mu=rep(0,dim(D)[1]),Sigma=D)
                      thiscoef<-coef(mod)[sample(1:mod$dims$ngrps[1],1),]
                      if (length(formals(paste("HD",modelName,sep="")))==5) {
                         predall[data$plot==plots[i]]<-
                             eval(call(paste("HD",modelName,sep=""),
                                  d=data$d[data$plot==plots[i]],a=thiscoef[[1]],b=thiscoef[[2]],c=thiscoef[[3]],bh=bh))
                         } else {
                         predall[data$plot==plots[i]]<-
                             eval(call(paste("HD",modelName,sep=""),
                                  d=data$d[data$plot==plots[i]],a=thiscoef[[1]],b=thiscoef[[2]],bh=bh))
                         }
                                  
            #          predall[data$plot==plots[i]] <- HDnaslund(d=data$d[data$plot==plots[i]],a=thiscoef[[1]],b=thiscoef[[2]])
                      }
                  }
                      
               # If addResidual==True, a random residual is added from a normal distribution 
               # using the estimated variance function.
               if (addResidual) {
                  if (varf) {
                     ressd<-mod$sigma*data$d^mod$modelStruct$varStruct
                     } else if (!varf) {
                     ressd<-mod$sigma
                     }
                  predall<-predall+rnorm(length(data$d),0,ressd)
                  }
                  
               type[!is.na(h)]<-0
               h[is.na(h)]<-predall[is.na(h)]
              
               list(h=h,imputed=is.na(data$h),model=mod,predType=type,hpred=predall)
               
               }
			   
			   
# Predicts tree volumes (liters) using the functions of Laasasenaho.
# pl=species (1= pine)
# malli mean model. Value of 1 uses only d; balue of 2 uses d and h.
predvol<-function(species,d,h=0,model=1) {
	               pl<-species
				   malli<-model
				   pl2<-pl
				   pl2[pl2>4|pl2<1]<-5
				   if (malli==1) {
					   pars<-matrix(c(-5.39417,3.48060,0.039884,-5.39934,3.46468,0.0273199,
									   -5.41948,3.57630,0.0395855,-5.41948,3.57630,0.0395855,
									   rep(NA,3)),ncol=3,byrow=TRUE)  
					   exp(pars[pl2,1]+pars[pl2,2]*log(2+1.25*d)-pars[pl2,3]*d)
					   
				   } else if (malli==2) {
					   pars<-matrix(c(0.036089,2.01395,0.99676,2.07025,-1.07209,
									   0.022927,1.91505,0.99146,2.82541,-1.53547,
									   0.011197,2.10253,0.98600,3.98519,-2.65900,
									   0.011197,2.10253,0.98600,3.98519,-2.65900,
									   rep(NA,5)),ncol=5,byrow=TRUE)
					   pars[pl2,1]*d^pars[pl2,2]*pars[pl2,3]^d*h^pars[pl2,4]*(h-1.3)^pars[pl2,5]
				   }
			   }
			   
# TODO: Add comment
# 
# Author: Lauri Mehtatalo and Jouni Siipilehto 27.1.2014
			   ###############################################################################
			   
# Siipilehto & Meht2talo 2013. "Parameter recovery vs. parameter prediction for the Weibull distribution
# validated for Scots pine stands in Finland". Silva Fennica 47(4), article id 1057; 
			   
# R-script for parameter recovery for the unweighted 2-parameter Weibull function assumed for diameter distribution
# Assume we know one of the following mean characteristics:
# A: D = arithmetic mean diameter i.e. the first moment
# B: DG = basal area weighted mean diameter i.e. mean of the weighted second order Weibull distribution 
# C: DM = arithmetic median diameter i.e. 0.5 quantile
# D: DGM = basal area median diameter i.e. 0.5 quantile of the weighted second order Weibull distribution
# The second moment (Dq^2) is calculated from the stem number (N_ha) and basal area per hectare (G_ha)
			   
# Basal area (G) is given as m^2/ha
# The number of stems (N) is given as 1/ha
# The mean diameter is given in cm.
			   
			   
# Functions returning the Weibull scale parameter for given shape and mean diameter, 
# for the four altermnative mean dianmeters.
# Functions DXXX return the first derivative of these functions, 
# but they are not currently in use and may have mistakes. 
			   
# A: Unweighted mean
			   scaleDMean<-function(D,shape) {
				   D/gamma(1/shape+1) 
			   }
			   
# B Basal-area weighted mean    
			   scaleDGMean<-function(D,shape) {
				   D*gamma(2/shape+1)/gamma(3/shape+1) 
			   }
			   
			   
# C: Unweighted median
			   scaleDMed<-function(D,shape) {
				   D/(log(2)^(1/shape))   
			   }  
			   
			   
# D: Basal-area weighted median
			   scaleDGMed<-function(D,shape) {
				   D/(qgamma(0.5,2/shape+1))^{1/shape}      
			   }
			   
# The function to be set to zero: 
# Returns the difference of the given shape parameter and 
# the value corresponding to given stand chaacteristics.
# Dtype should be any of letters "A", "B", "C" and "D",
# And specifies which type of mean or median diameter was given 
# in argument "D".
			   
			   func.recweib<-function(lshape,G,N,D,Dtype, trace=FALSE) {  
				   shape<-exp(lshape)+0.01
				   if (Dtype=="A") {
					   scale<-scaleDMean(D,shape)
				   } else if (Dtype=="B") {
					   scale<-scaleDGMean(D,shape)
				   } else if (Dtype=="C") {
					   scale<-scaleDMed(D,shape)
				   } else if (Dtype=="D") {
					   scale<-scaleDGMed(D,shape)
				   } else stop("Dtype should be any of 'A', 'B', 'C', 'D'")
				   
				   sol<-40000*G/(pi*N)-scale^2*gamma(2/shape+1) 
				   attributes(sol)$scale<-scale
				   if (trace) cat(shape,attributes(sol)$scale,sol,"\n")
				   sol
			   }
			   
			   
			   
# The function for recovery of shape and scale.
# G is Basal area, m^2/ha
# N Number of stems, 1/ha
# Weighted or unweighted D is mean or median diameter, cm
# Dtype the type of the given mean diameter, given as "A", "B", "C" and "D".
# A: The arithmetic mean 
# B: The basal-area weighted mean
# C: The median 
# D: The median weighted by basal area
# trace: if TRUE, then intermediate outbut during estimation is printed on the screen.
# init: the initial guess for the shape parameter. 
# If not given, the initial guesses are based on the PPM equations given 
# in the Appendix of Siipilehto and Mehtatalo (2013).
# 
# The recovery is based on the solution of the equation 
# shape(G,N,D)-shape = 0
# Where shape(G,N,D) expresses the shape parameter of The weibull distribution
# for the given values of the stand characteristics and shape is the estimated
# value of the shape parameter.
# 
# The function returns a list with the following elements
# shape: the shape parameter (c) in the solution
# scale: the scale parameter (b) in the solution 
# G, N, D, Dtype: the values of the input arguments
# val: the value of the equation shape(G,N,D)-shape at the solution
# 
			   recweib<-function(G,N,D,Dtype,init=NA,trace=FALSE) {
				   if (!sum(Dtype!=c("A","B","C","D"))) stop("Dtype should be any of 'A', 'B', 'C', 'D'")
				   
				   DQM<-sqrt(40000*G/(pi*N)) 
				   
				   if (is.na(init)) {
					   if (Dtype=="A") {
						   init<-1.0/log(DQM^4/D^4 + 0.1)
					   } else if (Dtype=="B") {
						   init<-1/log(D^2/DQM^2 + 0.05)
					   } else if (Dtype=="C") {
						   init<-1/log(DQM^4/D^4 + 0.1)
					   } else if (Dtype=="D") {
						   init<-1/log(D^2/DQM^2 + 0.05)
					   } 
				   }
				   
				   lshape<-NRnum(init=log(init-0.01), 
						   fnlist=list(function(lshape) func.recweib(lshape,G=G,N=N,D=D,Dtype=Dtype,trace=trace)),
						   crit=10)
				   shape=exp(lshape$par)+0.01
				   if (Dtype=="A") {
					   scale<-scaleDMean(D,shape)
				   } else if (Dtype=="B") {
					   scale<-scaleDGMean(D,shape)
				   } else if (Dtype=="C") {
					   scale<-scaleDMed(D,shape)
				   } else if (Dtype=="D") {
					   scale<-scaleDGMed(D,shape)
				   } 
				   
				   list(shape=exp(lshape$par)+0.01, scale=scale, G=G, N=N, D=D, Dtype=Dtype, val=lshape$value)
			   }
			   