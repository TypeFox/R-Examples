nrr <- function (data_x, logband = TRUE)
{
    dm = ncol(data_x)
    data_num = nrow(data_x)
    h = vector(, dm)
    for (j in 1:dm)
    {
        temp = sum(data_x[, j])
        temp2 = sum(data_x[, j] * data_x[, j])
        sigma = sqrt(temp2/data_num - (temp/data_num) * (temp/data_num))
        temp = exp(1/(dm + 4) * log(4/(dm + 2)))
        h[j] = temp * sigma * exp(-1/(dm + 4) * log(data_num))
    }
    ifelse(logband == TRUE, return(log(h)), return(h))
}
