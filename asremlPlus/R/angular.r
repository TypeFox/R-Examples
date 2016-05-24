"angular" <- function(proportions, n)
{180/pi*asin(sqrt(proportions + (proportions==0)/(4*n) + (proportions==1)*((n-0.25)/n-1)))}
"angular.mod" <- function(count, n)
{180/pi*asin(sqrt((count+0.375)/(n+0.75)))}

