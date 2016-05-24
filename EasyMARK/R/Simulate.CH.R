Simulate.CH <-
function (surv.form, p.form, p.constant = NULL, surv.constant = NULL, 
    N = 100, max.occ = 100, noise = 0.2) 
{
    if (is.null(surv.constant)) {
        phi = 1/(1 + exp(-(surv.form + rnorm(N, mean = 0, sd = noise))))
    }
    else {
        cat("Ignoring surv.form, since a constant is supplied\n")
        phi = rep(surv.constant, N)
    }
    if (is.null(p.constant)) {
        p = 1/(1 + exp(-(p.form + rnorm(N, mean = 0, sd = noise))))
    }
    else {
        cat("Ignoring p.form, since a constant is supplied\n")
        p = rep(p.constant, N)
    }
    h.list = list()
    for (i in 1:N) {
        j = 1
        h = c(1)
        alive = h
        while (alive[j] == 1) {
            alive[j + 1] <- rbinom(1, 1, phi[i] * alive[j])
            h[j + 1] <- rbinom(1, 1, p[i] * alive[j + 1])
            if (length(h) >= max.occ) 
                break
            j = j + 1
        }
        h.list[[i]] = h
    }
    n.obs = sapply(h.list, length)
    seq.max = seq_len(max(n.obs))
    ch_split = t(sapply(h.list, "[", i = seq.max))
    ch_split[is.na(ch_split)] = 0
    ch = apply(ch_split, 1, paste, collapse = "")
    ch = data.frame(ch = ch)
    ch["ch"] = as.character(ch$ch)
    return(list(ch = ch, ch_split = ch_split, p = p, phi = phi))
}
