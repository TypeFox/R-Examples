print.SemiParSampleSel <- function(x,...){

  if(x$margins[2]=="N")      {nn <- "Gaussian"                 ; m2l <- "identity"}
  if(x$margins[2]=="G")      {nn <- "gamma"                    ; m2l <- "log"}
  if(x$margins[2]=="P")      {nn <- "Poisson"                  ; m2l <- "log"}
  if(x$margins[2]=="NB")     {nn <- "negative binomial type I" ; m2l <- "log"}
  if(x$margins[2]=="D")      {nn <- "Delaporte"                ; m2l <- "log"}
  if(x$margins[2]=="PIG")    {nn <- "Poisson inverse Gaussian" ; m2l <- "log"}
  if(x$margins[2]=="S")      {nn <- "Sichel"                   ; m2l <- "log"} 
  if(x$margins[2]=="BB")      {nn <- "beta binomial"              ; m2l <- "logit"}
  if(x$margins[2]=="BI")      {nn <- "binomial"                   ; m2l <- "logit"}
  if(x$margins[2]=="GEOM")    {nn <- "geometric"                  ; m2l <- "log"}
  if(x$margins[2]=="LG")      {nn <- "logarithmic"                ; m2l <- "logit"}
  if(x$margins[2]=="NBII")    {nn <- "negative binomial type II"  ; m2l <- "log"}
  if(x$margins[2]=="WARING")  {nn <- "Waring"                     ; m2l <- "log"}
  if(x$margins[2]=="YULE")    {nn <- "Yule"                       ; m2l <- "log"}
  if(x$margins[2]=="ZIBB")    {nn <- "zero inflated beta binomial"  ; m2l <- "logit"}
  if(x$margins[2]=="ZABB")    {nn <- "zero altered beta binomial"   ; m2l <- "logit"}
  if(x$margins[2]=="ZABI")    {nn <- "zero inflated binomial"       ; m2l <- "logit"}
  if(x$margins[2]=="ZIBI")    {nn <- "zero altered binomial"        ; m2l <- "logit"}
  if(x$margins[2]=="ZALG")    {nn <- "zero altered logarithmic"        ; m2l <- "logit"}
  if(x$margins[2]=="ZINBI")   {nn <- "zero inflated negative binomial"   ; m2l <- "log"}
  if(x$margins[2]=="ZANBI")   {nn <- "zero altered negative binomial"    ; m2l <- "log"}
  if(x$margins[2]=="ZIP")     {nn <- "zero inflated Poisson type I"      ; m2l <- "log"}
  if(x$margins[2]=="ZAP")     {nn <- "zero altered Poisson"              ; m2l <- "log"}
  if(x$margins[2]=="ZIP2")    {nn <- "zero inflated Poisson type II"     ; m2l <- "log"}
  if(x$margins[2]=="ZIPIG")   {nn <- "zero inflated Poisson inverse Gaussian"      ; m2l <- "log"}
  
  
  if(x$BivD=="N")   {cop <- "Bivariate Normal"                      ;lind <- "atanh(theta)"}
  if(x$BivD=="C0")  {cop <- "Clayton Copula"                        ;lind <- "log(theta)"}
  if(x$BivD=="C90") {cop <- "Rotated Clayton Copula (90 degrees)"   ;lind <- "log(-theta)"}
  if(x$BivD=="C180"){cop <- "Survival Clayton Copula"		    ;lind <- "log(theta)"}
  if(x$BivD=="C270"){cop <- "Rotated Clayton Copula (270 degrees)"  ;lind <- "log(-theta)"}
  if(x$BivD=="J0")  {cop <- "Joe Copula"                            ;lind <- "log(theta-1)"}
  if(x$BivD=="J90") {cop <- "Rotated Joe Copula (90 degrees)"	    ;lind <- "log(-theta-1)"}
  if(x$BivD=="J180"){cop <- "Survival Joe Copula"	            ;lind <- "log(theta-1)"}
  if(x$BivD=="J270"){cop <- "Rotated Joe Copula (270 degrees)"      ;lind <- "log(-theta-1)"}
  if(x$BivD=="FGM") {cop <- "FGM Copula"                            ;lind <- "atanh(theta)" }
  if(x$BivD=="F")   {cop <- "Frank Copula"                          ;lind <- "theta" }
  if(x$BivD=="AMH") {cop <- "AMH Copula"                            ;lind <- "atanh(theta)"}
  if(x$BivD=="G0")  {cop <- "Gumbel Copula"                         ;lind <- "log(theta-1)"}
  if(x$BivD=="G90") {cop <- "Rotated Gumbel Copula (90 degrees)"    ;lind <- "log(-theta-1)"}
  if(x$BivD=="G180"){cop <- "Survival Gumbel Copula"	            ;lind <- "log(theta-1)"}
  if(x$BivD=="G270"){cop <- "Rotated Gumbel Copula (270 degrees)"   ;lind <- "log(-theta-1)"}


  cat("\nERRORS' DISTRIBUTION:",cop)
  
  cat("\n\nSELECTION EQ.")
  cat("\nFamily: Bernoulli") 
  cat("\nLink function: probit")
  cat("\nFormula: "); print(x$gam1$formula)

  cat("\n") 
  cat("OUTCOME EQ.")
  cat("\nFamily:",nn) 
  cat("\nLink function:",m2l)
  cat("\nFormula: "); print(x$gam2$formula)
  
  
  
  
  if(length(x$formula) == 3){
  cat("\n") 
  cat("EQUATION 3")
  cat("\nLink function:",lind)
  cat("\nFormula: "); print(x$gam3$formula)
  
  } 
  
  if(length(x$formula) == 4){
  cat("\n") 
  cat("EQUATION 3")
  if(x$margins[2] %in% c("N", "G", "NB", "PIG", "WARING", "BB", "NBII")) cat("\nLink function:","log(sigma)")
  if(x$margins[2] %in% c("ZIP", "ZAP", "ZIP2", "ZALG", "ZIBI", "ZABI")) cat("\nLink function:","plogis(sigma)")
  cat("\nFormula: "); print(x$gam3$formula)
  
  cat("\n") 
  cat("EQUATION 4")
  cat("\nLink function:",lind)
  cat("\nFormula: "); print(x$gam4$formula)
  
  }  
  
    if(length(x$formula) == 5){
    cat("\n") 
    cat("EQUATION 3")
    cat("\nLink function:","log(sigma)")
    cat("\nFormula: "); print(x$gam3$formula)
    
    cat("\n") 
    cat("EQUATION 4")
    
    if(x$margins[2]=="S") cat("\nLink function:","identity(nu)")
    if(x$margins[2] %in% c("D", "ZINBI", "ZANBI", "ZIPIG", "ZIBB", "ZABB")) cat("\nLink function:","qlogis(nu)")
    cat("\nFormula: "); print(x$gam4$formula)
    
    cat("\n") 
    cat("EQUATION 5")
    cat("\nLink function:",lind)
    cat("\nFormula: "); print(x$gam4$formula)    
    
  } 
  
      
            
    
  cat("\n") 


  cp <- ("theta = ")

if(x$margins[2] %in% c("N", "NB", "PIG", "BB", "NBII", "WARING", "ZABI", "ZIBI", "ZALG", "ZAP", "ZIP", "ZIP2", "D", "S", "ZIBB", "ZABB", "ZANBI", "ZINBI", "ZIPIG"))
    {dis <-  "  sigma = "; diss <- x$sigma.a}
if(x$margins[2]=="G"){dis <-  "  sigma = "; diss <- x$sigma.a}
if(x$margins[2] %in% c("D", "S", "ZIBB", "ZABB", "ZANBI", "ZINBI", "ZIPIG")){dis1 <-  "  nu = "; diss1 <- x$nu.a}


if(x$margins[2] %in% c("G", "N", "NB", "PIG", "BB", "NBII", "WARING", "ZABI", "ZIBI", "ZALG", "ZAP", "ZIP", "ZIP2")){
cat("n = ",x$n,"  n.sel = ",x$n.sel,dis,round(diss,3),"\n",cp,round(x$theta.a,3),"  total edf = ",round(x$t.edf,3),"\n\n",sep="")
} 

if(x$margins[2] %in% c("P", "BI", "GEOM", "LG", "YULE")){
cat("n = ",x$n,"  n.sel = ",x$n.sel,"\n",cp,round(x$theta.a,3),"  total edf = ",round(x$t.edf,3),"\n\n",sep="")
}


if(x$margins[2] %in% c("D", "S", "ZIBB", "ZABB", "ZANBI", "ZINBI", "ZIPIG")){
cat("n = ",x$n,"  n.sel = ",x$n.sel,dis,round(diss,3),dis1,round(diss1,3),"\n",cp,round(x$theta.a,3),"  total edf = ",round(x$t.edf,3),"\n\n",sep="")
}  

  invisible(x)

}
















