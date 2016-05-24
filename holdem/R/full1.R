full1 <-
function(x){
a1 = mycount1(x)
a2 = sort(a1$ct,decreasing=TRUE)
if(length(a2) < 1.5) return(0)
if(a2[1] < 2.5) return(0)
if(a2[2] < 1.5) return(0)
a3 = min(c(1:length(a1$ct))[a1$ct > 2.5])
a4 = a1$v[a3] ## the number of the trip
a5 = a1$ct[-a3]
a6 = a1$v[-a3]
a7 = min(c(1:length(a5))[a5 > 1.5]) ## the number of the pair
a8 = a6[a7]
15*a4 + a8
}  ## end of full1

