
pass<-function(data, base="LASSO", lambda.grid=NULL, num.grid=20, num.split=20, alpha=0.1){
           
           n<-nrow(data); p<-ncol(data)-1
           centerized=function(x){x-mean(x)}
           data=apply(data, 2, centerized) ### Dataset is centerized

           lambda.max<-2*max(t(data[, 1:p])%*%data[,(p+1)])
           n1<-as.integer(n/2); n2<-n-n1

	     if (length(lambda.grid)==0) lambda.grid<-10^(seq(log10(0.001), log10(lambda.max), length.out=num.grid))  
           num.grid<-length(lambda.grid)
           
           pass.vec<-rep(0, num.grid); kappa.vec<-rep(0, num.grid); cv.vec<-rep(0, num.grid)
           kappa.mat<-matrix(0, num.split, num.grid); cv.mat<-matrix(0,num.split, num.grid)

           for(i.split in 1:num.split){
				set.seed(i.split)  
              for(i.grid in 1:num.grid){ 
                data.perm=data[sample(1:n, n, replace=FALSE), ] 	
                data1=data.perm[1:n1, ]
                data2=data.perm[-(1:n1), ]
                
                beta1.hat<-rep(0, p); beta2.hat<-rep(0, p)

                if (base=="LASSO"){ ### Get penalized coeffecies using LASSO based on sub-datasets
                   beta1.hat<-predict(lars(data1[, 1:p], data1[, (p+1)], use.Gram=FALSE, intercept=FALSE),
                                      s=lambda.grid[i.grid], type="coef", mode="lambda")$coefficients
                   beta2.hat<-predict(lars(data2[, 1:p], data2[,(p+1)], use.Gram=FALSE, intercept=FALSE),
                                      s=lambda.grid[i.grid], type="coef", mode="lambda")$coefficients
	          }

                if (base=="aLASSO"){ ### Get penalized coeffecies using adaptive LASSO based on sub-datasets
                   out.lm1=lm.ridge(data1[,(p+1)]~data1[, 1:p])
                   out.lm2=lm.ridge(data2[,(p+1)]~data2[, 1:p])
                   coef.lm1=out.lm1$coef; weight1=1/abs(coef.lm1[abs(coef.lm1)>0])
                   coef.lm2=out.lm2$coef; weight2=1/abs(coef.lm2[abs(coef.lm2)>0])
                   data1.xx=scale(data1[, 1:p], center=FALSE, scale=weight1)
                   data2.xx=scale(data2[, 1:p], center=FALSE, scale=weight2)   
                   beta1.hat=predict(lars(data1.xx, data1[, (p+1)], normalize=FALSE, intercept=FALSE, use.Gram=FALSE),
                                     s=lambda.grid[i.grid], type="coef", mode="lambda")$coefficients
                   beta2.hat=predict(lars(data2.xx, data2[, (p+1)], normalize=FALSE, intercept=FALSE, use.Gram=FALSE),
                                     s=lambda.grid[i.grid], type="coef", mode="lambda")$coefficients
                   beta1.hat<-beta1.hat*coef.lm1; beta2.hat<-beta2.hat*coef.lm2
                 }
                
                if (base=="SCAD"){ ### Get penalized coeffecies using SCAD based on sub-datasets  
                   fit1=ncvreg(data1[, 1:p], data1[,(p+1)], family="gaussian", penalty="SCAD", lambda=lambda.grid)
                   fit2=ncvreg(data2[, 1:p], data2[,(p+1)], family="gaussian", penalty="SCAD", lambda=lambda.grid)
                   beta1.hat=predict(fit1,data1[, 1:p],which=i.grid,type="coefficients")[-1]
                   beta2.hat=predict(fit2,data2[, 1:p],which=i.grid,type="coefficients")[-1]
                 }
                 
                 kappa.mat[i.split, i.grid]<-agree.twosets(which(as.vector(beta1.hat)!=0), which(as.vector(beta2.hat)!=0), p)
                 cv.mat[i.split, i.grid]<-cv.twosets(data1, beta1.hat, data2, beta2.hat)     
             }   
           }

           kappa.vec<-apply(kappa.mat, 2, mean)
           cv.vec<-apply(cv.mat, 2, mean)
           pass.vec<-kappa.vec/cv.vec     

           i.lambda.opt.kappa<-which(kappa.vec>=max(kappa.vec)*(1-alpha))[1]
           i.lambda.opt.pass<-which(pass.vec==max(pass.vec))[1]

           beta.hat.kappa<-rep(0, p)
           beta.hat.pass<-rep(0, p)


                  if (base=="LASSO"){ ### Get penalized coeffecies using LASSO based on sub-datasets
                             beta.hat.kappa<-predict(lars(data[, 1:p], data[, (p+1)], use.Gram=FALSE, intercept=FALSE),
                                      s=lambda.grid[i.lambda.opt.kappa], type="coef", mode="lambda")$coefficients
                             beta.hat.pass<-predict(lars(data[, 1:p], data[,(p+1)], use.Gram=FALSE, intercept=FALSE),
                                      s=lambda.grid[i.lambda.opt.pass], type="coef", mode="lambda")$coefficients
	                }

                  if (base=="aLASSO"){ ### Get penalized coeffecies using adaptive LASSO based on sub-datasets
                      out.lm=lm.ridge(data[,(p+1)]~data[, 1:p])
                      coef.lm=out.lm$coef; weight=1/abs(coef.lm[abs(coef.lm)>0])
                      data.xx=scale(data[, 1:p], center=FALSE, scale=weight) 
                      beta.hat.kappa=predict(lars(data.xx, data[, (p+1)], normalize=FALSE, intercept=FALSE, use.Gram=FALSE),
                                        s=lambda.grid[i.lambda.opt.kappa], type="coef", mode="lambda")$coefficients
                      beta.hat.pass=predict(lars(data.xx, data[, (p+1)], normalize=FALSE, intercept=FALSE, use.Gram=FALSE),
                                        s=lambda.grid[i.lambda.opt.pass], type="coef", mode="lambda")$coefficients
                      beta.hat.kappa<-beta.hat.kappa*coef.lm; beta.hat.pass<-beta.hat.pass*coef.lm
                     }
                
                  if (base=="SCAD"){ ### Get penalized coeffecies using SCAD based on sub-datasets  
                     fit=ncvreg(data[, 1:p], data[,(p+1)], family="gaussian", penalty="SCAD", lambda=lambda.grid)
                     beta.hat.kappa=predict(fit,data[, 1:p],which=i.lambda.opt.kappa,type="coefficients")[-1]
                     beta.hat.pass=predict(fit,data[, 1:p],which=i.lambda.opt.pass,type="coefficients")[-1]  
                     }



           results<-list(lambda.grid=lambda.grid, alpha=alpha, base=base, 
                          pass.values=pass.vec, kappa.values=kappa.vec,
                          lambda.pass=lambda.grid[i.lambda.opt.pass], 
                          lambda.kappa=lambda.grid[i.lambda.opt.kappa], 
                          beta.kappa=as.numeric(beta.hat.kappa), beta.pass=as.numeric(beta.hat.pass),   
                          subset.kappa=which(as.vector(beta.hat.kappa)!=0), subset.pass=which(as.vector(beta.hat.pass)!=0) 
                          )

           class(results)<-"pass"
           return(results)
          }

print.pass<-function(x, ...){
      cat("Lambda is evaluated over a grid of ", length(x$lambda.grid), 
           " values from ", min(x$lambda.grid), " to ", max(x$lambda.grid), ".",  fill=TRUE)
      cat("\n")
	cat("Optimal lambda selected by Kappa criterion with alpha=", x$alpha, "is: \n") 
      print(x$lambda.kappa) 
      cat("\n")
      cat("Optimal lambda selected by PASS criterion is: \n") 
      print(x$lambda.pass)  
      cat("\n")
      cat("The estimated coefficients using selected lambda by Kappa criterion are: \n")
      print(x$beta.kappa) 
      cat("\n")
      cat("The estimated coefficients using selected lambda by PASS criterion are: \n")
      print(x$beta.pass) 
      cat("\n")
      cat("The selected submodel by Kappa criterion is: \n")
      print(x$subset.kappa)  
      cat("\n")
      cat("The selected submodel by PASS criterion is: \n")
      print(x$subset.pass)
}

plot.pass<-function(x, ...){
     par(mfrow=c(2,1)) 
     plot(log(x$lambda.grid), x$kappa.values, type="l", ylab="Kappa scores", xlab=expression(paste("log(", lambda, ")")),  
           main=paste("Kappa Criterion for", x$base, "Procedure")) 
     plot(log(x$lambda.grid), x$pass.values, type="l", ylab="PASS scores", xlab=expression(paste("log(", lambda, ")")), 
           main=paste("PASS Criterion for", x$base, "Procedure"))
}



###########################
###### Sub-functions ######
###########################


### Some libraries needed

# library(MASS)
# library(lars)
# library(ncvreg)


agree.twosets=function(aset1,aset2,p.tot){ ### Caluculate Kappa coefficient of two sets
    if(length(aset1)==0||length(aset2)==0||length(aset1)+length(aset2)==2*p.tot) ratio=-1 else 
    {n11=length(intersect(aset1,aset2))
     n22=length(intersect(setdiff(c(1:p.tot),aset1),setdiff(c(1:p.tot),aset2)))
     n12=length(setdiff(aset1,intersect(aset1,aset2)))
     n21=length(setdiff(aset2,intersect(aset1,aset2)))
     ratio=((n11+n22)/p.tot-((n11+n12)*(n11+n21)+(n12+n22)*(n21+n22))/(p.tot*p.tot))/(1-((n11+n12)*(n11+n21)+(n12+n22)*(n21+n22))/(p.tot*p.tot))
    }
    return(ratio)
}

cv.twosets=function(data1, beta1.hat, data2, beta2.hat){ ### Caluculate two-fold cross-validation 
        n1<-nrow(data1); n2<-nrow(data2); ncol<-ncol(data1)
        data1.x<-data1[,-ncol]; data1.y<-data1[,ncol]
        data2.x<-data2[,-ncol]; data2.y<-data2[,ncol]

        RSS=function(y,yhat){sum((y-yhat)^2)} ### Residual sum of squares 
 
        cv.value=(RSS(data1.y, data1.x%*%beta2.hat)/n1 + RSS(data2.y, data2.x%*%beta1.hat)/n2)/2
        return(cv.value)
}









