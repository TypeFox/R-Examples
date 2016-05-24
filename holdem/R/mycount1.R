mycount1 <-
function(x){
    ## returns sorted unique values of x and how many times each appears
b1 = sort(unique(x),decreasing=TRUE)
b2 = length(b1)
b3 = rep(0,b2)
for(i in 1:b2) b3[i] = sum(x == b1[i])
list(v=b1, ct = b3)
}

