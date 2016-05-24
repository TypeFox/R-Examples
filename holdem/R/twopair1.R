twopair1 <-
function(x){
## 225*highpair + 15*lowpair + next
a1 = mycount1(x)
a2 = a1$v
a3 = a1$ct
a4 = sum(a3>1.5)
if(a4<1.5) return(0)
a5 = sort(a2[a3>1.5],decreasing=TRUE)
a6 = max(c(0,a2[(a2 != a5[1]) & (a2 != a5[2])]))
225*a5[1] + 15*a5[2] + a6
}    ## end of twopair1

