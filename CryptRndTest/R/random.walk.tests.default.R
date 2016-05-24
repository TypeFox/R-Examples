random.walk.tests.default=function(x,B=64,Excursion=TRUE,Expansion=TRUE,Height=TRUE,alpha=0.05){
  options(warn=-1)
  check(test=5,x=x,B=B,alpha=alpha,Excursion=Excursion,Expansion=Expansion,Height=Height)
  
  res.tbl=random.walk.tests.main(x,B=64,Excursion=TRUE,Expansion=TRUE,Height=TRUE,alpha=0.05)
  res.tbl$call = match.call()
  class(res.tbl) = c("random.walk.tests","CryptRndTest")
  res.tbl
}