
wald.test<-function(out,statistic="diff",quiet=FALSE,alpha=0.01)
{	
  y<-out$y
  gname=out$gname
  f.z.y<-out$f.z.y
  a<-out$a
  lambda<-out$lambda
  w<-out$w
  cr=out$cr
  nj<-table(cr)
  p<-dim(y)[1]
  d<-max(out$cr)
  ycond<-matrix(0,p,d)
  K=out$K
  
  
  for(j in 1:d) ycond[,j]<-rowSums(y[,cr==j,drop=FALSE])
  
  ycond<-ifelse(ycond==0,1e-10,ycond)
  lambda<-ifelse(lambda==0,1e-10,lambda)
  
  var.y<-matrix(0,p,d) ##  var(Y_{ij+})
	for(j in 1:d) var.y[,j]<-nj[j]*lambda[,j]*(1+rowSums(f.z.y/matrix(a,p,K,byrow=TRUE))*lambda[,j])*nj[j]/(nj[j]-1)
  var.lambda<-var.y*(1/nj[j]^2)
  
  var.loglambda<-NULL

##############
# Difference #
##############
  if(statistic=="diff"){
    var.somma<-rowSums(var.lambda)
    varianza<-var.somma
    stat<-(lambda[,1]-lambda[,2])/sqrt(var.somma)
    pvalue<-2*(pnorm(abs(stat),lower.tail=F))
                    
    }

#########
# Ratio #
#########
  if(statistic=="ratio"){
    var.ratio<-rep(0,p)
    Elambda<-lambda
    var.ratio<-(var.lambda[,1]/Elambda[,2]^2)+(Elambda[,1]^2*var.lambda[,2]/Elambda[,2]^4)
    varianza<-var.ratio
    stat<-(lambda[,1]/lambda[,2]-1)/sqrt(var.ratio)
    pvalue<-2*(pnorm(abs(stat),lower.tail=F))
      }
  
############
# LogRatio #
############
  if(statistic=="logratio"){
    var.loglambda<-var.y*(1/ycond^2)
    var.sommalog<-rowSums(var.loglambda)
    varianza<-var.sommalog
    stat<-(log(lambda[,1]/lambda[,2]))/sqrt(var.sommalog)
    pvalue<-2*(pnorm(abs(stat),lower.tail=F))
    }
  
  pvalueadj=p.adjust(pvalue, method = "fdr")
    
if(! quiet){
  if (is.null(rownames(y))) rownames(y)=paste("g",gname,sep="")
	nomi.colonna= c("Statistic","p-value","adj. p-value","Variance")
        df=data.frame(cbind(round(stat,2),round(pvalue,4),round(pvalueadj,4),round(varianza,2)),row.names=rownames(y))
        index=pvalueadj<alpha
        colnames(df)=nomi.colonna
        message("")
        message("DE genes:")
        print(df[index,])
   }
out=(list(stat=stat,pvalue=pvalue,pvalueadj=pvalueadj,var=varianza,gname=gname))
invisible(out)
}


