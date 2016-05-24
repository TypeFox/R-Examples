getKcutoff <-
function (m, n, s = 0.45014) 
{
    mu = (sqrt(m - 0.5) + sqrt(n - 0.5))^2
    sigma = (sqrt(m - 0.5) + sqrt(n - 0.5)) * (1/sqrt(m - 0.5) + 
        1/sqrt(n - 0.5))^(1/3)
    return(sqrt(mu/n + s * sigma/n))
}
