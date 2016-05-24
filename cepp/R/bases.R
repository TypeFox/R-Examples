basis_random <- function (n, d = 2) 
{
    mvn <- matrix(rnorm(n * d), ncol = d)
    orthonormalise(mvn)
}

basis_nearby <- function(alpha = 0.75, method = 'geodesic', d = 2)
{
    function (current)
        {
            current <- matrix(current, ncol = d)
            new <- basis_random(nrow(current), d)
            switch(method,
                   linear = orthonormalise((1 - alpha) * current + alpha * new),
                   geodesic = step_fraction(geodesic_info(current, new), alpha))
        }
}
