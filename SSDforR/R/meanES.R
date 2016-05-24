meanES <-
  function(es,lab,esmain){
    
    
esmean<-mean(es,na.rm=T)
essd<-sd(es,na.rm=T)
writeLines("-----------mean-------------") 
print(esmean)
writeLines("------------SD--------------")
print(essd)
writeLines("------------% change--------------")
dchange=pnorm(esmean)-.5
l6<-c("% change=",round(dchange,4)*100)
print(c(l6))
writeLines("********************************************************")
l1<-c("small effect size: <.87")
l2<-c("medium effect size: .87 to 2.67 ")
l3<-c("large effect size: >2.67")
writeLines(l1)
writeLines(l2)
writeLines(l3)
e<-data.frame(es,lab)
e<-e[ order(e$es), ]

#dotchart(es,groups=lab,color="red",cex=.8, xlab="Cohen's D",main=esmain)
dotchart(e$es,labels=e$lab,color="red",cex=.8, xlab="Cohen's D",main=esmain)
    
  }   
