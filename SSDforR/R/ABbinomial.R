ABbinomial <-
function(phaseX,v1,v2,successA,successB){
  t1<-table(phaseX)
  tmaxA<-t1[names(t1)==v1]
  tmaxB<-t1[names(t1)==v2]
  psucessA=successA/tmaxA
  b=binom.test(x=successB, n=tmaxB, p=psucessA)
  print(b)
 
}
