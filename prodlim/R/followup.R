### followup.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Sep 22 2015 (10:29) 
## Version: 
## last-updated: Sep 25 2015 (06:19) 
##           By: Thomas Alexander Gerds
##     Update #: 2
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
followup <- function(formula,data,...){
    G <- prodlim(formula,data,reverse=TRUE)
    quantile(G,...)
}


#----------------------------------------------------------------------
### followup.R ends here
