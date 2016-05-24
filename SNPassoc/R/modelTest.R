`modelTest` <-
function(X,Y,quantitative,type,genotypingRate)
{
  control<-ifelse(6%in%type,5,length(type)) 
  controlGeno <- GenotypeRate(X)
  if (genotypingRate > controlGeno)
   {
    ans<-c("Genot error",rep(NA,control))
   }
  else
   {
    if (is.Monomorphic(X))
      ans<-c("Monomorphic",rep(NA,control))
    else { 
     ans<-NA  
     if (1%in%type | 6%in%type) { 
       mco<-assoc(Y,codominant(X),quantitative=quantitative) 
       ans<-c(ans,mco)
     }
     if (2%in%type | 6%in%type) {
       mdo<-assoc(Y,dominant(X),quantitative=quantitative) 
       ans<-c(ans,mdo)
     }
     if (3%in%type | 6%in%type) {
       mre<-assoc(Y,recessive(X),quantitative=quantitative) 
       ans<-c(ans,mre)
     }
     if (4%in%type | 6%in%type) {
       mov<-assoc(Y,overdominant(X),quantitative=quantitative) 
       ans<-c(ans,mov)
     }
     if (5%in%type | 6%in%type) {
      mad<-assoc(Y,additive(X),quantitative=quantitative) 
      ans<-c(ans,mad)
     }
    }
   }
 ans  
}

