feature.assoc <-
function(g,y,method="correlation",g1=NULL,exact=1){
        if(is.matrix(g)!=T){stop("g is not a matrix")}
        if(length(y)!=ncol(g)){stop("y does not have as many elements as g has columns")}
        if(class(y)=="character") {stop("y cannot be a character")}
        if(is.factor(y)&length(levels(y))<1.5){stop("y has to have more than one level")}
        if(!is.null(g1)){if(!is.matrix(g1)){stop("g1 is not a matrix")} } 
        if(!is.null(g1)){if(!identical(dim(g),dim(g1))){stop("g1 does not have the same dimensions as g")} }   
    
    if(is.null(g1)){
    g1<-g
    for (i in 1:nrow(g)){
    g1[i,]<-sample(g[i,],ncol(g),replace=F)}
    }
    
    
    if(is.factor(y)){
          if(length(levels(y))==2){
              if(method=="correlation"){
              ycor<-as.numeric(y)
              observed<-apply(g,1,cor,y=ycor,use="complete.obs")
              observed.p<-apply(g,1,function(x) cor.test(x,y=ycor,use="complete.obs")$p.value)
              permuted<-apply(g1,1,cor,y=ycor,use="complete.obs")
              permuted.p<-apply(g1,1,function(x) cor.test(x,y=ycor,use="complete.obs")$p.value)
              }
              
              if(method=="t.test"){
              lev1<-levels(y)[2]
              lev2<-levels(y)[1]
              observed<-vector(mode = "numeric", length = length(rownames(g)))
              observed.p<-vector(mode = "numeric", length = length(rownames(g)))
              names(observed)<-rownames(g)
              names(observed.p)<-rownames(g)
              for (vad in 1:length(rownames(g))){
              x1 = g[vad,y==lev1]
              x2 = g[vad,y==lev2]
              tested<-t.test(x1,x2)
              observed[vad] = tested$statistic
              observed.p[vad] = tested$p.value
              }                           
              permuted<-vector(mode = "numeric", length = length(rownames(g)))
              permuted.p<-vector(mode = "numeric", length = length(rownames(g)))
              names(permuted)<-rownames(g)
              names(permuted.p)<-rownames(g)
              for (vad in 1:length(rownames(g))){
              x1 = g1[vad,y==lev1]
              x2 = g1[vad,y==lev2]
              tested<-t.test(x1,x2)
              permuted[vad] = tested$statistic
              permuted.p[vad] = tested$p.value
              } 
              }
              
              if(method=="AUC"){
              lev1<-levels(y)[2]
              lev2<-levels(y)[1]
              observed<-vector(mode = "numeric", length = length(rownames(g)))
              observed.p<-vector(mode = "numeric", length = length(rownames(g)))
              names(observed)<-rownames(g)
              names(observed.p)<-rownames(g)
              for (v1 in 1:length(rownames(g))){
              x1 = g[v1,y==lev1]
              n1 = sum(!is.na(x1))
              x2 = g[v1,y==lev2]
              n2 = sum(!is.na(x2))
              r = rank(c(x1,x2),na.last=NA)
              observed[v1] = (sum(r[1:n1]) - n1*(n1+1)/2) / (n1*n2)
              observed.p[v1] = wilcox.test(na.omit(x1),na.omit(x2),exact=exact)$p.value
              }
              
              permuted<-vector(mode = "numeric", length = length(rownames(g)))
              permuted.p<-vector(mode = "numeric", length = length(rownames(g)))
              names(permuted)<-rownames(g)
              names(permuted.p)<-rownames(g)
              for (v1 in 1:length(rownames(g))){
              x1 = g1[v1,y==lev1]
              n1 = sum(!is.na(x1))
              x2 = g1[v1,y==lev2]
              n2 = sum(!is.na(x2))
              r = rank(c(x1,x2),na.last=NA)
              permuted[v1] = (sum(r[1:n1]) - n1*(n1+1)/2) / (n1*n2)
              permuted.p[v1] = wilcox.test(na.omit(x1),na.omit(x2),exact=exact)$p.value
              }  
              }
              }
                    
          if(length(levels(y))>2){
              method<-"linear.model"
              observed<-vector(mode = "numeric", length = length(rownames(g)))
              observed.p<-vector(mode = "numeric", length = length(rownames(g)))
              names(observed)<-rownames(g)
              names(observed.p)<-rownames(g)
              for (v1 in 1:length(rownames(g))){
              s<-summary(lm(g[v1,]~y))
              observed[v1]<-s$r.squared
              observed.p[v1]<-pf(s$fstatistic[1],s$fstatistic[2],s$fstatistic[3],lower.tail=FALSE)
              }
              
              permuted<-vector(mode = "numeric", length = length(rownames(g)))
              permuted.p<-vector(mode = "numeric", length = length(rownames(g)))
              names(permuted)<-rownames(g)
              names(permuted.p)<-rownames(g)
              for (v1 in 1:length(rownames(g))){
              s<-summary(lm(g1[v1,]~y))
              permuted[v1]<-s$r.squared
              permuted.p[v1]<-pf(s$fstatistic[1],s$fstatistic[2],s$fstatistic[3],lower.tail=FALSE)
              }             
              }
          
          }
              
    if(is.numeric(y)){    
          method<-"correlation"
              observed<-apply(g,1,cor,y=y,use="complete.obs")
              observed.p<-apply(g,1,function(x) cor.test(x,y=y,use="complete.obs")$p.value)
              permuted<-apply(g1,1,cor,y=y,use="complete.obs")
              permuted.p<-apply(g1,1,function(x) cor.test(x,y=y,use="complete.obs")$p.value)
          }
    return(list(observed=observed,permuted=permuted,observed.p=observed.p, permuted.p=permuted.p,method=method,class.of.y=class(y),permuted.data=g1))      
    }

