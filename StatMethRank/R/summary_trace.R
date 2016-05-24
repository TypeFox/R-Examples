# Create a output list by beta_trace and Sigma_trace
summary_trace <- function(beta_trace, Sigma_trace, item_name = NULL)
{
    # Input  beta_trace(nTrace * p), Sigma_trace(nTrace * (k-1)^2), X((nTrace*(k-1)) * p)
    # Return beta_trace, Sigma_trace, rs_trace, V_trace, rv_trace
    nTrace = dim(beta_trace)[1]
    p = dim(beta_trace)[2]
    k = sqrt(dim(Sigma_trace)[2]) + 1
    if (is.null(item_name))
    {
        item_name = 1:k
    }
    # Generate rs_trace
    rs_trace = t(apply(Sigma_trace, 1, sigma2corr))
    # Generate V_trace from beta_trace, Sigma_trace    
    V_trace = t(apply(Sigma_trace, 1, sigma2v))
    rV_trace = t(apply(V_trace, 1, sigma2corr))


    beta_trace = data.frame(beta_trace)
    for (i in 1:p)
    {
        names(beta_trace)[i] = paste("beta", item_name[i], sep = "")
    }

    temp_ind = c(t(upper.tri(matrix(1:(k-1)^2, ncol = k - 1), diag = TRUE)))
    Sigma_trace = data.frame(Sigma_trace[, temp_ind])
    for (i in 1:(k-1))
    {
        for (j in i:(k-1))
        {
            names(Sigma_trace)[(i - 1) * (k - 1) + j - (i - 1) * i / 2] = paste("s", i, j, sep = "")
        }
    }
    
    temp_ind = c(t(upper.tri(matrix(1:(k-1)^2, ncol = k - 1), diag = FALSE)))
    rs_trace = data.frame(rs_trace[, temp_ind])
    for (i in 1:(k-2))
    {
        for (j in (i+1):(k-1))
        {
            names(rs_trace)[(i - 1) * (k - 2) + (j - 1) - (i - 1) * i / 2] = paste("r", i, j, sep = "")
        }
    }    

    temp_ind = c(t(upper.tri(matrix(1:k^2, ncol = k), diag = TRUE)))
    V_trace = data.frame(V_trace[, temp_ind])
    for (i in 1:k)
    {
        for (j in i:k)
        {
            names(V_trace)[(i - 1) * k + j - (i - 1) * i / 2] = paste("v", item_name[i], item_name[j], sep = "")
        }
    }

    temp_ind = c(t(upper.tri(matrix(1:k^2, ncol = k), diag = FALSE)))
    rV_trace = data.frame(rV_trace[, temp_ind])
    for (i in 1:(k-1))
    {
        for (j in (i+1):k)
        {
            names(rV_trace)[(i - 1) * (k - 1) + (j - 1) - (i - 1) * i / 2] = paste("r", item_name[i], item_name[j], sep = "")
        }
    }

    return(
            list(
                    beta_trace = beta_trace,
                    beta_msci = msci(beta_trace),
                    Sigma_trace = Sigma_trace,
                    Sigma_msci = msci(Sigma_trace),
                    rs_trace = rs_trace,
                    rs_msci = msci(rs_trace),
                    V_trace = V_trace,
                    V_msci = msci(V_trace),
                    rV_trace = rV_trace,  
                    rV_msci = msci(rV_trace)                  
                )   
        )
} # end of summary.trace