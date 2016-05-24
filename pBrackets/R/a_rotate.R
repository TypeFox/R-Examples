a_rotate <-
function(M, p)
{
    M %*% t(matrix(c(cos(p), sin(p), -sin(p), cos(p)), 2, 2))
}
