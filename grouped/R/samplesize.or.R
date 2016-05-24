"samplesize.or" <-
function(n, MC, theta, sigma, m, a, dist.t, dist.x, grouping.mech){
    n <- if(n %% 2) n + 1 else n
    if(missing(dist.x)) marg.power.orig(MC, n, m, theta, sigma, a, dist.t, grouping.mech)[[1]] else
    marg.power.orig(MC, n, m, theta, sigma, a, dist.t, dist.x, grouping.mech)[[1]]
}
