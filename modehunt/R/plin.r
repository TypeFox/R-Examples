plin <- function(q, a, b, s){

## Values of a density that is uniform on [0,1], except on
## intervals given by vectors a (= start points) and b (= endpoints) 
## where it is equal to a linear function with slopes s.

if (length(a) == 1){
    if (b < a){cat('Interval is not proper!')}} else {
        for (i in 2:length(a)){if (b[i-1] > a[i]){cat('Vectors a and b are not proper!')}}}
    
for (i in 1:length(a)){if (s[i] > 2/(b[i]-a[i])){cat('Not a density for this choice of a, b, s!')}}

    y <- q
    
    for (i in 1:length(q)){
        for (j in 1:length(a)){
            if (q[i] >= a[j] & q[i] < b[j]){y[i] <- q[i] + s[j] / 2 * (q[i]^2 - a[j]^2+(a[j]-q[i])*(a[j]+b[j]))}
            }
        }
return(y)
}
