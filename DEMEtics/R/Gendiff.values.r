# Function used within calc.r to calculate final measures of genetic differentiation/fixation
Gendiff.values <- function(Hs.values2.1tab,Ht.values2.1tab,sample.sizes2.1tab,x){
Hs.actual <- as.numeric(as.vector(Hs.values2.1tab$Hs.value))
Ht.actual <- as.numeric(as.vector(Ht.values2.1tab$Ht.value))
sample.size<-as.numeric(as.vector(sample.sizes2.1tab$sample.size))

if (x=="D"){one.locus<-D(Hs.actual,Ht.actual,sample.size)}
if (x=="Dest"){one.locus<-Dest(Hs.actual,Ht.actual,sample.size)}
if (x=="Gst"){one.locus<-Gst(Hs.actual,Ht.actual,sample.size)}
if (x=="Gst.est"){one.locus<-Gst.est(Hs.actual,Ht.actual,sample.size)} 
one.locus
}
