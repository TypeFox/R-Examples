"samplesize.approx" <-
function(n, p., MC, theta, sigma, m, a, dist.t, dist.x, grouping.mech){
    if(missing(dist.x)) {p. - samplesize.ap(round(n), MC, theta, sigma, m, a, dist.t, grouping.mech)} else{
    p. - samplesize.ap(round(n), MC, theta, sigma, m, a, dist.t, dist.x, grouping.mech)}
}
