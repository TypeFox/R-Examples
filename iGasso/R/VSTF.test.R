VSTF.test <-
function(G, method="Anscombe")
{
    n1 = sum(G[1,])
    n2 = sum(G[2,])

    F1 = function(X, n) acos(1-2*X/n)/2    # = asin(2*X/n-1)/2 + pi/2
    F2 = function(X, n)  F1(X+3/8, n+3/4)
    F3 = function(X, n)  F1(X, n+1) + F1(X+1, n+1)
    F4 = function(X, n)  { c = 1/(4*n); asin(sqrt(c+(1-2*c)*X/n)) }

    if (method == "arcsine")  S = 4*(F1(G[1,2], n1) - F1(G[2,2], n2))^2/(1/n1 + 1/n2)
    else if (method == "Anscombe")  S = 4*(F2(G[1,2], n1)-F2(G[2,2], n2))^2/(1/(n1+0.5)+1/(n2+0.5))
    else if (method == "Freeman-Tukey")  S = (F3(G[1,2], n1)-F3(G[2,2], n2))^2/(1/(n1+0.5)+1/(n2+0.5))
    else if (method == "Chanter")  S = 4*(F4(G[1,2], n1) - F4(G[2,2], n2))^2/(1/n1 + 1/n2)
    else stop("Valid transformation is arcsine, Anscombe, Freeman-Tukey, or Chanter")

	names(S) = "statistic"
    structure(list(statistic = S, p.value = 1-pchisq(S, 1),
                   method = paste("Comparing Two Proportions Using the", method, "Transformation"), 
                   data.name = deparse(substitute(G))), class = "htest") 
}

