SCD<-function(fac)   {
 cnames<-paste("x",1:fac, sep="")
 b<-c(rep(0,fac))
  for (i in 1:fac) {
    cb<-combn(fac,i)
    t<-dim(cb)
    e<-t[1]
    r<-t[2]
        if (i==1) {
      SC<-rbind(b,diag(r))
                  } else {
       for (j in 1:r) {
         v<-b
         for (k in 1:e) {
           v[cb[k,j]]<-1/e
                        }
           SC<-rbind(SC,v)
                      }
                  }
                   }
dimS<-dim(SC)
rows<-dimS[1]
SC<-SC[2:rows, ]
colnames(SC)<-cnames
rnames<-paste(1:(rows-1))
rownames(SC)<-rnames
SC<-data.frame(SC)
return(SC)
                     }

