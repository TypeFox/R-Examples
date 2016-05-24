summary.ClaimsRule <-
function(object, ...)
{
   x<-object
   cat("\n")
   cat("Assignment according to the", x$Method, "rule for an Endowment of",x$E,  "\n")
   M<-cbind(x$Claims,x$R)
   colnames(M)<-c("Claims",x$Short)
   rownames(M)<-x$Names
   cat("\n")
   print(M)
   }
