"dgt" <-
function(x, mu = 0.0, sigma = 1.0, df = stop("no df arg"))
    dt((x - mu) / sigma, df = df) / sigma

