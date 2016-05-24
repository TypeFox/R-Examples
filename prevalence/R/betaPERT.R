betaPERT <-
function(a, m, b, k = 4, method = c("classic", "vose")) {
  ## check input
  if (!exists("a")) stop("'a' is missing")
  if (!exists("m")) stop("'m' is missing")
  if (!exists("b")) stop("'b' is missing")
  if (!is.numeric(a)) stop("'a' must be a numeric value")
  if (!is.numeric(m)) stop("'m' must be a numeric value")
  if (!is.numeric(b)) stop("'b' must be a numeric value")

  if (!exists("method")) stop("'method' is missing")
  method <- match.arg(method)

  if (method == "classic") {
    if (!exists("k"))
      stop("'k' is missing")
    if (!is.numeric(k))
      stop("'k' must be a numeric value")
    mu <- (a + k * m + b) / (k + 2)
    sdev <- (b - a) / (k + 2)
    alpha <- ((mu - a) / (b - a)) * ( ((mu - a) * (b - mu) / (sdev^ 2 )) - 1 )
    beta <- alpha * (b - mu) / (mu - a)
  }

  if (method == "vose") {
    if (!exists("k"))
      stop("'k' is missing")
    if (!is.numeric(k))
      stop("'k' must be a numeric value")
    mu <- (a + k * m + b) / (k + 2)
    alpha <- ifelse(mu == m,
                    1 + k / 2,
                    ((mu - a) * (2 * m - a - b)) / ((m - mu) * (b - a)))
    beta <- alpha * (b - mu) / (mu - a)
  }

  out <- list(alpha = alpha, beta = beta,
              a = a, m = m, b = b,
              method = method)
  class(out) <- "betaPERT"

  return(out)
}