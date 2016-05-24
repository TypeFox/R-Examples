seats<-function (N, M, r0 = 2.5) {
    radii <- seq(r0, 1, len = M)
    counts <- numeric(M)
    pts = do.call(rbind, lapply(1:M, function(i) {
        counts[i] <<- round(N * radii[i]/sum(radii[i:M]))
        theta <- seq(0, pi, len = counts[i])
        N <<- N - counts[i]
        data.frame(x = radii[i] * cos(theta), y = radii[i] * 
            sin(theta), r = i, theta = theta)
    }))
    pts = pts[order(-pts$theta, -pts$r), ]
    pts
}

election<-function (seats, result, formula,
 colours = sample(rainbow(length(counts)))) {
    result = model.frame(formula, result)
    counts = result[, 2]
    stopifnot(sum(counts) == nrow(seats))
    seats$party = factor(rep(result[, 1], counts))
    seats$colour = colours[as.numeric(seats$party)]
    seats
}
