flush1 <-
function(x,y){
    ## Can break down if dealing with 10 cards or more, i.e. if 2 flushes are possible
    ## returns 15^4 times top number + 15^3 times 2nd-highest + ...
a1 = mycount1(y)
if(max(a1$ct)<4.5) return(0)
a2 = c(1:length(a1$ct))[a1$ct > 4.5]
a3 = a1$v[a2] ## this is the suit
a4 = sort(x[y == a3],decreasing=TRUE)
sum(15^c(4:0) * a4[1:5])
} ## end of flush1

