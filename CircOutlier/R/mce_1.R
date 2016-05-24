MCe <-
function (u) 
{
    v = matrix(matrix(u, length(u), length(u))[-seq(1, length(u)^2, 
        length(u) + 1)], length(u) - 1)
    apply(v, 2, function(b) 1 - mean(b))
}
