chiMerge <-
function (data, varcon, alpha = 0.1,out=c("symb","num")) 
{
    n = dim(data)[1]
    p = dim(data)[2]
    class = dim(table(data[, p]))
    discredata = data
    threshold = qchisq(1 - alpha, class - 1)
    for (j in varcon) {
        z = sort(unique(data[, j]))
        midpoint = rep(0, (length(z) + 1))
        midpoint[1] = min(z)
        midpoint[length(z) + 1] = max(z)
        for (k in 2:length(z)) {
            midpoint[k] = (z[k - 1] + z[k])/2
        }
        cutpoint = midpoint[-c(1, length(z) + 1)]
        for (l in 1:(length(z) + 5)) {
            a = cut(data[, j], breaks = midpoint, include.lowest = TRUE)
            b = table(a, data[, p])
            m = dim(b)[1]
            if (m == 1) 
                break
            test = rep(0, (m - 1))
            d = matrix(0, 2, class)
            for (k in 1:(m - 1)) {
                d[1, ] = as.vector(b[k, ])
                d[2, ] = as.vector(b[k + 1, ])
                test[k] = tchisq(d)
            }
            if (min(test) < threshold) {
                index = rep(0, m)
                serie = seq(1, (m - 1), by = 1)
                for (k in 1:(m - 1)) {
                  if (test[k] == min(test)) {
                    index[k] = serie[k]
                  }
                }
                index = unique(sort(index))[-1]
                outpoint = cutpoint[min(index)]
                for (k in 1:length(midpoint)) {
                  if (midpoint[k] == outpoint) {
                    midpoint[k] = 0
                  }
                }
                midpoint = sort(midpoint)[-1]
                cutpoint = cutpoint[-min(index)]
            }
            else {
                break
            }
            if(out=="num")
            discredata[, j] = cut(data[, j], breaks = midpoint, 
                include.lowest = TRUE, label = FALSE)
            else {
                max=length(midpoint)
                discredata[,j]=cut(data[,j],breaks=c(-Inf,midpoint[-c(1,max)],Inf))
            }
        }
    }
    discredata
}
