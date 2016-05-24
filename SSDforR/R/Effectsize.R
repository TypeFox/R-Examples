Effectsize <-
function(behavior,phaseX,v1,v2) {
  mean1 <-tapply(behavior,phaseX,mean,na.rm=T)
  s1 <-tapply(behavior,phaseX, sd, na.rm=T)
 
  tphase<-table(phaseX)
  n1<- tphase[names(tphase)==v1]
  n2<- tphase[names(tphase)==v2]
  
  n1<-n1-1 
  n2<-n2-1
  totaln<-n1+n2
  DIFF<-mean1[names(mean1)==v1]-mean1[names(mean1)==v2]
 
  S<-sqrt(((n1*(s1[names(s1)==v1]^2)+(n2*s1[names(s1)==v2]^2))/totaln))
  CD<-(DIFF/S)
 
  cm=gamma((totaln/2))/(sqrt(totaln/2)*gamma((totaln-1)/2))
  G<-CD*cm
   t1<-table(phaseX)
  tmaxA<-t1[names(t1)==v1]
  startA<-match(v1,phaseX)
  endA<-tmaxA+startA-1
  A<-behavior[startA:endA]
  PA<-phaseX[startA:endA]
         
  tmaxA<-t1[names(t1)==v2]
  startA<-match(v2,phaseX)
  endA<-tmaxA+startA-1
  B<-behavior[startA:endA]
   PB<-phaseX[startA:endA]
  
  
IV<-c(PA,PB)
  DV<-c(A,B)
  
  
 reg1<- summary(lm(DV~IV))
  rvalue2<-reg1$r.squared
  rvalue<-sqrt(reg1$r.squared)
  
  
  DIFF<-mean1[names(mean1)==v2]-mean1[names(mean1)==v1]
  S<-s1[names(s1)==v1]
  es<-(DIFF/S)
  
  l1<-c("small effect size: <.87")
  l2<-c("medium effect size: .87 to 2.67 ")
  l3<-c("large effect size: >2.67")
  writeLines(l1)
  writeLines(l2)
  writeLines(l3)
  writeLines("********************************************************")
  writeLines("********************ES**********************************")
  
  es1<-round(es,5)
  pes1<-c("ES=        ",es1)
  eschange=pnorm(es)-.5
  l5<-c("% change=",round(eschange,4)*100)
  print(c(pes1,l5))
  
  
  
  writeLines("*****************d-index*******************************")
  cd1<-(round(abs(CD),5))
  pcd1<-c("d-index=   ",cd1)
  dchange=pnorm(CD)-.5
  l6<-c("% change=",round(dchange,4)*100)
  
  print(c(pcd1,l6))
  
 writeLines("*****************Hedges's g****************************")
  hchange=pnorm(G)-.5
G1<-(round(abs(G),5))
PG1<-c("Hedges's g=",G1)
l7<-c("% change=",round(hchange,4)*100)

print(c(PG1,l7))

  writeLines("*****************Pearson's r***************************")
  print(round(rvalue,3))
  writeLines("*****************R-squared*****************************")
  print(round(rvalue2,3))
  a<-readline("(s)ave, (a)ppend, or (n)either results? (s/a or n) ")
  
  ES<-cd1
 
  ES=data.frame(ES)
  
  if (a=="s")
    {Label<-readline("Enter a behavior variable label ")
    ES<-data.frame(ES,Label)
  write.csv(ES,file = tclvalue(tcl("tk_getSaveFile")),row.names=FALSE)
  
   } 

   if (a=="a")
   { Label<-readline("Enter a behavior variable label ")
   ES<-cd1
     ES<-data.frame(ES,Label)
     writeLines("*****************open file to append to***************************")
     effA<-read.table(file.choose(),header=TRUE,sep=',') 
            out=rbind(effA,ES)
            writeLines(" ")
            writeLines(" ")
            writeLines(" ")
            writeLines("*****************save appended file***************************")            
      write.csv(out,file = tclvalue(tcl("tk_getSaveFile")),row.names=FALSE) }
   
  
}
