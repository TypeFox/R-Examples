`haplo.interaction` <-
function(formula, data, SNPs.sel, quantitative = is.quantitative(formula, data), 
             haplo.freq.min=0.05,...)

{
   if (!inherits(data, "setupSNP")) 
        stop("data must be an object of class 'setupSNP'")
 
   control.SNPs<-sum(!is.na(match(names(data),SNPs.sel)))
   if (control.SNPs!=length(SNPs.sel))
     stop("Some of the SNPs selected are not in the data set")
    
   cl <- match.call()
   mf <- match.call(expand.dots = FALSE)
   m0 <- match(c("formula", "data", "subset"), names(mf), 0)
   mf <- mf[c(1, m0)]
   mf[[1]] <- as.name("model.frame")
   mf <- eval(mf, parent.frame())
   mt <- attr(mf, "terms")
   
   special <- c("int")
   Terms <- if (missing(data)) 
        terms(formula, special)
    else terms(formula, special, data = data)
   posInt <- attr(Terms, "specials")$int
   
   if(length(posInt))
    {
      var2<-mf[,posInt]
      if(!length(levels(var2))) 
       {
        stop("interaction variable must be a factor")
       }
    }  
   else 
      stop("formula needs an 'interaction' term")   
   
   control.missing<-dimnames(mf)[[1]]
   geno <- make.geno(data[dimnames(data)[[1]]%in%control.missing,], 
                   SNPs.sel)

   dep <- mf[, 1]
   if (ncol(mf) > 2)
    adj <- data.frame(mf[, -c(1,posInt)])
   else
    adj <- NULL 
   varAdj <- attr(mt, "term.labels")[-(posInt-1)]
   varInt <- attr(mt, "term.labels")[posInt-1]
   varInt <-gsub("int\\(","",varInt)
   varInt <-gsub("\\)","",varInt)

   out<-haplo.inter.fit(geno, var2, dep, adj , ifelse(quantitative,"gaussian","binomial"), haplo.freq.min, ...)
   
   res.corner<-out[[1]]
   xx<-dimnames(res.corner)[[2]] 
   xx[xx=="li"]<-"lower"
   xx[xx=="ls"]<-"upper"
   dimnames(res.corner)[[2]]<-xx

   temp<-out[[2]]   
   etiq1<-dimnames(temp)[[1]]
   aux0<-dimnames(temp)[[2]]
   etiq2<-aux0[seq(2,length(aux0),3)]
      
   ans<-list(NA)
   for (i in 1:length(etiq2))
       {
         ans[[i]]<-temp[,c(1,(2+3*(i-1)):(4+3*(i-1)))]
         ans[[i]][1,2]<-ifelse(quantitative,0,1)  
         ans[[i]]<-ans[[i]][,-1]  
         if (!quantitative)
           dimnames(ans[[i]])[[2]]<-c("OR","lower","upper")               
         else
            dimnames(ans[[i]])[[2]]<-c("diff","lower","upper")               
       }
   names(ans)<-etiq2      
   res.int1<-ans      

   temp<-out[[3]]   
   etiq1<-dimnames(temp)[[1]]
   aux0<-dimnames(temp)[[2]]
   etiq2<-aux0[seq(2,length(aux0),3)]
      
   ans2<-list(NA)
   for (i in 1:length(etiq1))
       {
         ans.i<- matrix(temp[i,][-1],nrow=length(etiq2) ,ncol=3,byrow=TRUE)
         ans2[[i]]<-data.frame(ans.i)
         dimnames(ans2[[i]])[[1]]<-etiq2 
         ans2[[i]][1,1]<-ifelse(quantitative,0,1)  
         if (!quantitative)
           dimnames(ans2[[i]])[[2]]<-c("OR","lower","upper")               
         else
           dimnames(ans2[[i]])[[2]]<-c("diff","lower","upper")               
       }
   names(ans2)<-etiq1      
   res.int2<-ans2      
  
   res<-list(res.corner,res.int1,res.int2,out$pval)
   attr(res,"label.snp")<-SNPs.sel
   attr(res,"varAdj")<-varAdj 
   attr(res,"varInt")<-varInt
   attr(res,"quantitative")<-quantitative
   class(res)<-"haploOut"
   res

}

