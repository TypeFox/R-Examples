Huberized <-
function (t) 
{
    par(mfrow = c(1, 2))
    n = length(t)
    m = mean.circular(circular(t))
    s = sd.circular(circular(t))
    k = t > m + (1.5 * s)
    k1 = t < m - (1.5 * s)
    plot(t, type = "l", sub = "plot1", main = paste("Outliers:", 
        paste((1:n)[k], (1:n)[k1], collapse = ",")))
    abline(h = m - (1.5 * s), col = 6)
    abline(h = m + (1.5 * s), col = 5)
    t[k] = m + (1.5 * s)
    t[k1] = m - (1.5 * s)
    plot(t, type = "l", sub = "plot2", main = paste("Outliers:", 
        paste((1:n)[k], (1:n)[k1], collapse = ",")))
    abline(h = m - (1.5 * s), col = 6)
    abline(h = m + (1.5 * s), col = 5)
    m1 = mean.circular(circular(t))
    s1 = sd.circular(circular(t)) * 1.134
    output = cbind(m = m, s = s, m1 = m1, s1 = s1)
    list(output = output)
}
