curvetest <-
function (formula,   data1=NULL, data2=NULL,  
   equal.var = TRUE, alpha =  0.5, bw=NULL, plotit = TRUE, conf.level = 0.05,  
    kernel = c("Triangle",  "Gaussian", "Trio","Uniform", "Triweight",  "Epanechnikov", 
    "Quartic"), nn=100, myx = NULL, bcorrect="simple",...) {
   if (missing(kernel)) kernel = "Quartic"
   if(length(kernel)>1) kernel = "Quartic"
   types<-c("Gaussian", "Trio",  "Uniform", "Triweight", "Triangle", "Epanechnikov", "Quartic")
   pm = pmatch(kernel, types)
   ww<-get.weight.function(types[pm]) 
   data1 <-model.frame(formula=formula, data=data1); 
   if(!is.null(data2)) data2 <-model.frame(formula=formula, data=data2); 
   if(missing(alpha)) alpha=0.5 
   if(is.null(data2)) myx<-seq(min(data1[,2]), max(data1[,2]), length=nn) else 
   if(is.null(myx)) myx<-seq(min(data1[,2], data2[,2]), max(data1[,2], data2[,2]), length=nn)    
   if(!is.null(alpha)){
     if(is.numeric(alpha)) alpha1= alpha else  
     if(alpha=="optimal")  alpha1= as.numeric(getoptimalalpha(formula=formula, data1)) 
   }else  alpha1=NULL
   f1<-curvefit(formula=formula, data=data1, kernel=kernel, alpha=alpha1,bw=bw,bcorrect=bcorrect, myx=myx)   
   f2<-NULL;
   if(!is.null(data2)) {
     if(!is.null(alpha)) {
     if(is.numeric(alpha)) alpha2= alpha else  
     if(alpha=="optimal")  alpha2= as.numeric(getoptimalalpha(formula=formula, data2)) 
     }else  alpha2=NULL
   f2<-curvefit(formula=formula, data=data2, kernel=kernel,  alpha=alpha2, bw=bw,myx=myx, bcorrect=bcorrect)
     out<-curvetest.raw(fits1=f1, fits2=f2,equal.var =equal.var,   
       conf.level = conf.level, plotit = plotit)
   }
   if(plotit) plot(x=f1, col=2,...)
   if(plotit&!is.null(data2)) plot(x=f2,add=T, col=3, ...)
   return(out)
}
