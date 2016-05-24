"samplesize.orig" <-
function(n, p., MC, theta, sigma, m, a, dist.t, dist.x, grouping.mech){
    if(missing(dist.x)) {p. - samplesize.or(round(n), MC, theta, sigma, m, a, dist.t, grouping.mech)} else {
    p. - samplesize.or(round(n), MC, theta, sigma, m, a, dist.t, dist.x, grouping.mech)}
}
