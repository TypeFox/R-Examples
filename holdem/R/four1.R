four1 <-
function(x){
## 15*number of the foursome + next
a1 = mycount1(x)
a2 = a1$v
a3 = a1$ct
a4 = sum(a3 > 3.5)
if(a4 < 0.5) return(0)
a5 = sort(a2[a3>3.5],decreasing=TRUE)
a6 = sort(c(0,x[(x != a5[1])]),decreasing=TRUE)
15*a5[1] + a6[1]
}    ## end of four1

