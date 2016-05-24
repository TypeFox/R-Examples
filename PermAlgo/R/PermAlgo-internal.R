`.form.data` <-
function (ordered.info, m, n, Xmat, XmatNames = NULL) 
{
    data <- data.frame(cbind(sort(rep(seq(1, n), m)), 
  rep(0, m * n), rep(ordered.info$obs.t, each = m), rep(seq(0, (m-1)), n), 
  rep(seq(0, (m-1)), n) + 1, Xmat))
    if (is.null(XmatNames) == TRUE) {
        XmatNames <- paste("X", seq(1, dim(Xmat)[2]), sep = "")
    }
    colnames(data) <- c("Id", "Event", "Fup", "Start", "Stop", 
        XmatNames)
    data$Event[data$Stop == data$Fup & rep(ordered.info$d, each = m) == 
        1] <- 1
    data <- subset(data, data$Stop <= data$Fup)
}


`.partialHazards` <-
function (t, v, covArray, betas) 
{
if (length(betas) > 1) {return(exp(covArray[v, t, ] %*% betas))} else {return(exp(covArray[v, t, ] * betas))}
}

`.PermuteCovariateVectors` <-
function (t, d, count, I, covArray, betas) 
{
    n = sum(count[I])
    p <- integer(n)
    v <- seq(n)
    ip = 1
    for (k in I) {
        if (d[k]) {
            J = sample(length(v), size = count[k], replace = FALSE, prob = .partialHazards(t[k], v, covArray, betas))
        }
        else {
            J = sample(length(v), size = count[k], replace = FALSE, 
                prob = NULL)
        }
        p[seq(from = ip, along.with = J)] <- v[J]
        ip <- ip + count[k]
        v = v[-J]
    }
    return(p)
}

