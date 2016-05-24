##' Least-Committed  Principle for creating bbas
##'
##' @export
##' @param Mat matrix, \eqn{m \times k}, \eqn{m} is the number of sources, \eqn{k} is the length of probability vectors. If the number of sources is 1, the input probability could be a vector. 
##' @return mass_bba matrix, \eqn{m \times 2^k}, each column is a bba. If there is only one source, the output is a bba vector.
##' @examples
##' pro1 = c(0.25, 0.25, 0.25, 0.25);
##' pro2 = c(0.3, 0.2, 0.2, 0.1);
##' pro3 = rbind(pro1, pro2);
##' 
##' LCPrincple(pro1)
##' LCPrincple(pro2)
##' LCPrincple(pro3)
##'

LCPrincple <- function(Mat) {
    
    # From probability to bbas using Least-Committed principle
    
    # Input: Mat, matrix, m*k, m is the number of sources, k is the length of pro.vectors.
    
    # Output: mass_bba, matrix, 2^k * m, each column is a bba.
    
    if (is.null(nrow(Mat))) {
        Mat = matrix(Mat, 1)
    }
    
    m = nrow(Mat)
    k = ncol(Mat)
    
    masses = c()
    ords = c()
    Label_vec = rep(0, k)
    mass_bba = matrix(0, m, 2^k)
    
    for (j in 1:m) {
        
        pro = Mat[j, ]
        mass = c()
        ord = c()
        step = 0
        # browser()
        while (step < k) {
            min_val = min(pro)
            min_order = which(pro == min_val)
            ord = c(ord, min_order)
            Label_vec[step + (1:length(min_order))] = step + 1
            value = min_val * (k - step)
            mass = c(mass, value)
            pro = pro - min_val
            pro[min_order] = Inf
            step = step + length(min_order)
        }
        
        delete_times = sort(unique(Label_vec))
        credal_order = 2^k
        
        if (length(delete_times) > 1) {
            temp1 = 1:k
            for (jj in 1:(length(delete_times) - 1)) {
                temp1 = setdiff(temp1, ord[Label_vec == delete_times[jj]])
                temp = sum(2^(temp1 - 1)) + 1
                credal_order = c(credal_order, temp)
            }
        }
        mass_bba[j, credal_order] = mass
    }
	mass_bba = t(mass_bba)
    if(ncol(mass_bba) == 1) mass_bba = as.vector(mass_bba) 
    return(mass_bba)
} 
