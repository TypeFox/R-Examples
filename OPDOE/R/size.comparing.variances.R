size.comparing.variances <- function(ratio,alpha,power)
 {
   g <- function(n){
     qf(1-alpha/2,n-1,n-1)*qf(power,n-1,n-1)-ratio
   }
   u <- uniroot(g,c(ratio,1000*ratio))$root
   ceiling(u)
 }

# Examples
#size.comparing.variances(5.5,0.05,0.8)
#[1] 13
#size.comparing.variances(2,0.05,0.9)
#[1] 90
#size.comparing.variances(2,0.05,0.8)
#[1] 68
