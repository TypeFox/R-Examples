`subsets` <-
function (n, r, v = 1:n)
if (r <= 0) vector(mode(v), 0) else if (r >= n) v[1:n] else {
    rbind(cbind(v[1], Recall(n - 1, r - 1, v[-1])), Recall(n -
        1, r, v[-1]))
}

