.getSB <- function(IC,neighbor)
   list(s = getRiskIC(IC,risk=trAsCov())$trAsCov$value^.5,
        b = getRiskIC(IC,risk=asBias(),neighbor=neighbor)$asBias$value)

getReq <- function(Risk,neighbor,IC1,IC2,n=1,upper=15){
            if(!is(IC1,"IC")||!is(IC2,"IC")) 
               stop("Arguments IC1, IC2 must be of class 'IC'.")
            if(!identical(IC1@CallL2Fam,IC2@CallL2Fam))
               stop("Arguments IC1, IC2 must be of defined for the same model.")
            sb1 <- .getSB(IC1,neighbor)
            sb2 <- .getSB(IC2,neighbor)
            if(abs(sb1$s-sb2$s)+ abs(sb1$b-sb2$b)<1e-6){
               cat(gettext("IC1 is just as good as IC2.\n"))
               return(c(0,Inf))
            }
            if((sb1$s<=sb2$s && sb1$b<sb2$b)||(sb1$s<=sb2$s && sb1$b<=sb2$b)){
               cat(gettext("IC1 is strictly better than IC2.\n"))
               return(c(0,Inf))}
            if((sb2$s<=sb1$s && sb2$b<sb1$b)||(sb2$s<sb1$s && sb2$b<=sb1$b)){
               cat(gettext("IC2 is strictly better than IC1.\n"))
               return(NA)}
            dRisk <- function(r){
                 get.asGRisk.fct(Risk)(r,s=sb1$s,b=sb1$b)-
                 get.asGRisk.fct(Risk)(r,s=sb2$s,b=sb2$b)
            }
            r0 <- uniroot(dRisk,lower=0, upper=upper)$root/n^.5
            if(sb1$s<=sb2$s)
               return(c(0,r0))
            else   
               return(c(r0,Inf))
            }

             

  