effect.size <-
function(lm.out){
   n<-length(lm.out$resid)
   p<-length(lm.out$coef)
   rsquared<-summary(lm.out)$r.sq 
   es<-matrix(ncol=2,nrow=6)
   rowNames=format(c("Wherry1","Claudy3","Smith","Wherry2","Olkin & Pratt","Pratt"),justify="left")
   colNames=format(c("Effect Size","Recommended"),justify="left")
   dimnames(es)<-list(rowNames,colNames)

####Compute Effect Sizes using Correction Formulas
####Wherry:
  es[1,1]= (round(1 - (((n-1)/(n-p-1)) * (1 - rsquared)),digits=4))
####Claudy3:
  es[2,1]=(round(1 - (((n-4)*(1-rsquared))/(n-p-1)) * (1 + (2*(1-rsquared))/(n-p+1)),digits=4))
####Smith:
  es[3,1]=(round(1-((n/(n-p)) * (1-rsquared)),digits=4))
####Wherry2:
  es[4,1]=(round(1- (((n-1)/(n-p)) * (1-rsquared)),digits=4))
####Olin & Pratt:
  es[5,1]=(round(1-(((n-3)*(1-rsquared))/(n-p-1))*(1 + (2*(1-rsquared))/(n-p+1)),digits=4))
####Pratt:
  es[6,1]=(round(1-(((n-3)*(1-rsquared))/(n-p-1))*(1 + (2*(1-rsquared))/(n-p-2.3)),digits=4))

####Identify recommended correction formula:
  es[,2]="No" 
  if (n<=30)es[6,2]="Yes"
  else 
    {
    if (n <=50) 
      {
      if (p<=2) es[6,2]="Yes"
      else 
        {
        if (p<=6) es[5,2]="Yes"
        else es[5,2]=es[6,2]="Yes"
        }
      }
    else
      {
      if (n<=80)
        {
        if (p<=2)es[6,2]=es[2,2]="Yes"
          else  ifelse((p<=6),es[1,2]<-"Yes", es[6,2]<-"Yes")
        }
      else
        {
        if (n <=150)
          {
          if (p<=2) es[4,2]<-"Yes"
          else es[6,2]<-"Yes"
          }
        else
          {
          if (p<=2)es[1,2]=es[3,2]="Yes"
          else ifelse((p<=6),es[2,2]<-"Yes", es[4,2]<-"Yes")
          }
        }
      }  
    }
  es<-data.frame(es)
  return(es)
}

