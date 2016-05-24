"bsearch" <-
function(last,nints,i,pd,stdv,ya,yb){
   tol <- 10^(-7)
   del <- 10
   uppr <- yb[i-1]
   q <- qp(uppr,last,nints[i-1],ya[i-1],yb[i-1],stdv)
   while (abs(q-pd) > tol){
      del <- del/10
      incr <- 2*as.integer(q > pd+tol)-1
      j <- 1
      while (j <= 50){
         uppr <- uppr+incr*del
         q <- qp(uppr,last,nints[i-1],ya[i-1],yb[i-1],stdv)
         if ({abs(q-pd) > tol}&{j==50}){
            stop("Error in search: not converging")
          }
         else if ({{incr==1}&{q <= pd+tol}}|{{incr==-1}&{q >= pd-tol}}){
            j <- 50
          }
         j <- j+1
       }
    }
   ybval <- uppr
   return(ybval)
 }

