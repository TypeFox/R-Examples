#'@title Checks Assumptions for Constructing Internal Instruments 
#'@description The internal instruments are constructed as in Lewbel(1997). See \code{\link{hmlewbel}}.
#'@keywords internal
#'@keywords lewbel 
checkAssumptions <- function(y,X,P,IIV,EIV=NULL, data=NULL)
{

 # check if model error is symetrically distributed
  
  e1 <- stats::lm(y~., data= data.frame(cbind(X,P)))$residuals

  
  if (abs(e1071::skewness(e1)) > 0.10 & IIV == "y2")  
    warning(gettextf("Model error not symetrically distributed"), domain = NA)
 

# check assumption E(qp) !=0
# save residuals of P regressed on 1 and X

  p1 <- stats::lm(P~., data=data.frame(X))$residuals

 IV <- internalIV(y,X,P,IIV)

 IVS <- cbind(X,IV)
# compute the mean of the product btw residuals and Instruments
 A12 <- round(apply(IVS*p1,2,mean),3)

for (i in length(A12)){
      if (abs(A12[i])<=0.1) 
       warning(gettextf("Assumptions E(qp)!=0 is not met"), domain = NA)
   }
}
