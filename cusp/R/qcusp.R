`qcusp` <- 
function (p, alpha, beta) 
{
    P = p
    p = p[(not1 <- p < 1) & (not0 <- p > 0)]
    f = function(t) pcusp(t, alpha, beta) - p
    b = a = p
    a[] = -10
    b = 10 
    fa = f(a)
    fb = f(b)
    c = co = b
    co[] = 0
    while (max(abs(c - co) > 1e-08)) {
        co = c
        c = (a + b)/2
        fc = f(c)
        iac = fa * fc < 0
        a = ifelse(iac, a, c)
        fa = ifelse(iac, fa, fc)
        b = ifelse(iac, c, b)
        fb = ifelse(iac, fc, fb)
    }
    Q = P
    Q[not1 & not0] = (a + b)/2
    Q[!not1] = Inf
    Q[!not0] = -Inf
    Q
}