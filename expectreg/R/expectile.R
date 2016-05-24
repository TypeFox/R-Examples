expectile <-
function (x, probs = seq(0, 1, 0.25), dec = 4) 
{
    if (!is.vector(x)) 
        stop("observations are needed in vector form.")
    if (min(probs) < 0 || max(probs) > 1) 
        stop("only asymmetries between 0 and 1 allowed.")
    e = mean(x)
    ee = 0 * probs
    g = max(abs(x)) * 1e-06
    for (k in 1:length(probs)) {
        p = probs[k]
        if (p == 0) 
            ee[k] = min(x, na.rm = TRUE)
        else if (p == 1) 
            ee[k] = max(x, na.rm = TRUE)
        else {
            for (it in 1:20) {
                w = ifelse(x < e, 1 - p, p)
                enew = sum(w * x)/sum(w)
                de = max(abs(enew - e))
                e = enew
                if (de < g) 
                  break
            }
            ee[k] = e
        }
    }
    names(ee) = probs
    ee = round(ee, dec)
    return(ee)
}
