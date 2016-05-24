Predict <-
function (x, y) 
{
    alpha = beta = u = c()
    alpha1 = 0.15
    beta1 = 1
    repeat {
        k1 = sum(x * sin(y - alpha1 - beta1 * x))
        k2 = sum((x^2) * cos(y - alpha1 - beta1 * x))
        beta2 = beta1 + (k1/k2)
        s = sum(sin(y - beta2 * x))
        c = sum(cos(y - beta2 * x))
        if (s > 0 & c > 0) {
            alpha2 = atan(s/c)
        }
        else if (c < 0) {
            alpha2 = atan(s/c) + pi
        }
        else if (s < 0 && c > 0) {
            alpha2 = atan(s/c) + 2 * pi
        }
        if (isTRUE(all.equal(c(alpha2, beta2), c(alpha1, beta1)))) 
            break
        else alpha1 = alpha2
        beta1 = beta2
    }
    output = cbind(alpha1, beta1)
    list(output = output)
}
