dlin <- function(x, a, b, s){

## Values of a density that is uniform on [0,1], except on
## intervals given by vectors a (= start points) and b (= endpoints)
## where it is equal to a linear function with slopes s.
## The resulting function is right-continuous.

if (length(a) == 1){if (b < a){cat('Interval is not proper!')}} else {
    for (i in 2:length(a)){if (b[i-1] > a[i]){cat('Vectors a and b are not proper!')}}}

for (i in 1:length(a)){if (s[i] > 2/(b[i]-a[i])){cat('Not a density for this choice of a, b, s!')}}

    y <- rep(1, length(x))

    for (i in 1:length(x)){
        for (j in 1:length(a)){
            if (x[i] >= a[j] & x[i] < b[j]){y[i] <- s[j] * x[i] - s[j] * (a[j] + b[j]) / 2 + 1}
            }
        }
return(y)
}
