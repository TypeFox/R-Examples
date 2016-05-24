hca.test <-
function(g,o,dendcut=2,method="correlation",link="ward",test="chisq",workspace=2e7){
   require(amap)
   if(is.data.frame(o)){
      if(!identical(rownames(o),colnames(g))){
        stop("Colnames of g are not the same as rownames of o")}
          }
    if(!is.data.frame(o)){
      o<-data.frame(o,row.names=colnames(g))}

    classes<-unlist(lapply(unclass(o),class))
   if(all(classes%in%c("factor","numeric","integer"))==F){stop("o can only contain factors and numeric")}
   if(test=="fisher"&ncol(g)>50){warning("more than 50 samples for Fisher's exact may get time-consuming, method chisq can be used instead")}
      

fit<-hcluster(t(g),method=method,link=link)
dendc<-as.factor(cutree(fit,dendcut))

pvals<-rep(NA,ncol(o))
names(pvals)<-colnames(o)
for (i in 1:ncol(o)){
  if(classes[i]=="factor"){
    if (test=="fisher"){pvals[i]<-fisher.test(table(dendc,o[,i]),workspace=workspace)$p.value}
    if (test=="chisq"){pvals[i]<-chisq.test(table(dendc,o[,i]))$p.value}
    }
  if(classes[i]%in%c("numeric","integer")){
          s<-summary(lm(o[,i]~dendc))
          pvals[i]<-pf(s$fstatistic[1],s$fstatistic[2],s$fstatistic[3],lower.tail=FALSE)
    }
  }
  return(list(p.values=pvals,classes=classes))
  }

