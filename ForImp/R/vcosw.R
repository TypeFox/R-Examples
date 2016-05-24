vcosw <-
function(v,w)
{
# inner product divided by the product of modules
sum(v*w)/(sqrt(sum(v^2)*sum(w^2)))
}

