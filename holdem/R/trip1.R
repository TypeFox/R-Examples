trip1 <-
function(x){
## 225*triple + 15*next + next
a1 = mycount1(x)
a2 = a1$v
a3 = a1$ct
a4 = sum(a3 > 2.5)
if(a4 < 0.5) return(0)
a5 = sort(a2[a3>2.5],decreasing=TRUE)
a6 = sort(c(0,0,x[(x != a5[1])]),decreasing=TRUE)
225*a5[1] + 15*a6[1] + a6[2]
}    ## end of trip1

