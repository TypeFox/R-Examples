strflsh1 <-
function(x,y){
a1 = mycount1(y)
if(max(a1$ct)<4.5) return(0)
a2 = c(1:length(a1$ct))[a1$ct > 4.5]
a3 = a1$v[a2] ## this is the suit
a4 = sort(x[y == a3],decreasing=TRUE)
straight1(a4)
} ## end of strflush1

