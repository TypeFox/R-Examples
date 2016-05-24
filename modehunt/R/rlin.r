rlin <- function(n, a, b, s){

## Values of a density that is uniform on [0,1], except on
## intervals given by vectors a (= start points) and b (= endpoints) 
## where it is equal to a linear function with slopes s.

if (length(a) == 1){if (b < a){cat('Interval is not proper!')}} else {
    for (i in 2:length(a)){if (b[i-1] > a[i]){cat('Vectors a and b are not proper!')}}}
for (i in 1:length(a)){if (s[i] > 2/(b[i]-a[i])){cat('Not a density for this choice of a, b, s!')}}

u <- sort(runif(n))
r <- qlin(u, a, b, s)
return(r)
}
