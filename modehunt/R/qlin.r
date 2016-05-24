qlin <- function(p, a, b, s){

## Values of a density that is uniform on [0,1], except on
## intervals given by vectors a (= start points) and b (= endpoints) 
## where it is equal to a linear function with slopes s.

if (length(a) == 1){if (b < a){cat('Interval is not proper!')}} else {
    for (i in 2:length(a)){if (b[i-1] > a[i]){cat('Vectors a and b are not proper!')}}}
for (i in 1:length(a)){if (s[i] > 2/(b[i]-a[i])){cat('Not a density for this choice of a, b, s!')}}

    y <- p
    
    for (i in 1:length(p)){
        for (j in 1:length(a)){
            if (p[i] >= a[j] & p[i] < b[j]){y[i] <- - 1 / s[j] + (a[j] + b[j]) / 2 + sign(s[j]) * 1/2*sqrt((a[j]-b[j])^2+4/s[j]*(1/s[j]-(a[j]+b[j])+2*p[i]))}
            }
        }
return(y)
}
