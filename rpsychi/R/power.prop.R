power.prop <- function(n, h=NULL, sig.level=0.05, one.sided=FALSE, prop=NULL){

#-----------------------------------------#
    prop.es <- function(p1, p2,one.sided=FALSE){
      phi1 <- 2 * asin(sqrt(p1)) #equation 6.2.1
      phi2 <- 2 * asin(sqrt(p2)) #equation 6.2.2
      h <- phi1 - phi2
      if(one.sided==FALSE){
      h <- abs(h)
      }
      return(h)
    }
#-----------------------------------------#

  if(one.sided==TRUE){
               both <- 1
                }else{
               both <- 2
  }

#-------one sample case----------#
if(length(n)==1){
  if(!is.null(h)){
           z1a <- qnorm(p=sig.level/both,lower.tail=FALSE)
           h <- h * sqrt(2)
           z1b <- h * sqrt(n/2) - z1a
           return(pnorm(z1b))
           }
   if(!is.null(prop)){
           p1 <- prop[1]; p2 <- prop[2]
           h <- prop.es(p1, p2, one.sided=one.sided)
           h <- h * sqrt(2)
           z1a <- qnorm(p=sig.level/both,lower.tail=FALSE)
           z1b <- h * sqrt(n/2) - z1a
           return(pnorm(z1b))
           }
}#end one sample

#-------two sample case----------#
if(length(n)==2){
  n1 <- n[1]; n2 <- n[2]
  n <- (2 * n1 * n2) / (n1 + n2)

  if(!is.null(h)){
           z1a <- qnorm(p=sig.level/both,lower.tail=FALSE)
           z1b <- h * sqrt(n/2) - z1a
           return(pnorm(z1b))
           }
   if(!is.null(prop)){
           p1 <- prop[1]; p2 <- prop[2]
           h <- prop.es(p1, p2, one.sided=one.sided)
           z1a <- qnorm(p=sig.level/both,lower.tail=FALSE)
           z1b <- h * sqrt(n/2) - z1a
           return(pnorm(z1b))
           }
}#end two sample

}
