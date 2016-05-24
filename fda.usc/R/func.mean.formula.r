func.mean.formula<-function(formula, data = NULL,...,drop=FALSE)
{
 tf <- terms.formula(formula)
 fac <- attr(tf, "term.labels")
 if (attr(tf, "response") > 0)  response <- as.character(attr(tf, "variables")[2])
 if (missing(formula) || (length(formula) != 3L))
        stop("'formula' missing or incorrect")
 ldata<-data
 data<-ldata$df
 if (is.null(ldata$df))   stop("'df' element is missing or incorrect")
# drop=TRUE
 if (is.vector(data)) data<-as.data.frame(data)
 if (is.matrix(data)) data<-as.data.frame(data)
   f<-ldata$df[[fac]]
   dat<-ldata[[response]]
   if (!is.factor(f)) f<-factor(f)
   nlev<-nlevels(f)
   lev<-levels(f)
   if (is.matrix(dat$data)) dat$data<-data.frame(dat$data)
#   out<-do.call("split.fdata",list("x"=dat,"f"=f,"drop"=drop))
#   for (i in 1:nlev) out[[lev[i]]]<-func.mean(out[[lev[i]]])
      out<-split(dat$data,f,drop=drop,...)    
 
     out2<-func.mean(fdata(out[[lev[1]]],dat$argvals,dat$rangeval,dat$names))     

     for (i in 2:nlev) out2<-c(out2,func.mean(fdata(out[[lev[i]]],dat$argvals,dat$rangeval,dat$names)))
#    for (i in 1:nlev) out[[lev[i]]]<-func.mean(fdata(out[[lev[i]]],dat$argvals,dat$rangeval,dat$names))
rownames(out2$data)<-lev
out2
}
     
