onepair1 <-
function(x){
## 15*15*15*pair + 15*15*next + 15*next + next
a1 = unique(x)
a2 = length(a1)
if(a2 == length(x)) return(0)
a3 = rep(0,a2)
for(i in 1:a2) a3[i] = sum(x == a1[i])
a4 = max(a1[a3>1.5])
a5 = sort(c(x[x != a4],rep(0,3)),decreasing=TRUE)
15*15*15*a4 + 225*a5[1] + 15*a5[2] + a5[3]
}

