## 'glm.control1' is a minor modification of 'glm.conrol'
## The only different is the addition of a ... argument
## Ioannis Kosmidis <I.Kosmidis@warwick.ac.uk> [15/02/2008]
`glm.control1` <-
function (epsilon = 1e-08, maxit = 25, trace = FALSE, ...) 
{
    glm.control(epsilon, maxit, trace)
}
