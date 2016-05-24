nothing1 <-
function(x){
## 15^4*highest + 15^3*next + 225*next + 15*next + next
y = c(x,rep(0,5))  ## this is in case x has length < 5
a1 = sort(y,decreasing=TRUE)
15*15*15*15*a1[1] + 15*15*15*a1[2] + 225*a1[3] + 15*a1[4] + a1[5]
}	## end of nothing1

