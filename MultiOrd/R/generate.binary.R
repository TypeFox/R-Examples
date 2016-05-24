generate.binary <-
function(no.rows, prop.vec.bin, corr.mat){
validation.CorrMat(prop.vec.bin, corr.mat) # conformity, symmetry, Pd and range is checked here
sigma_star=compute.sigma.star(prop.vec.bin, corr.mat) 

d = ncol(sigma_star)
data = rmvnorm(no.rows, mean = rep(0, d), sigma = sigma_star)
    p = prop.vec.bin
    q = 1 - p

for (i in 1:no.rows) {
            for (j in 1:d) {
                if (data[i, j] <= qnorm(1 - p[j])) 
                  data[i, j] = 0
                else data[i, j] = 1
        }
    }
    return(data)
}
