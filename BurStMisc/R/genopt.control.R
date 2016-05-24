"genopt.control" <-
function (births = 100, random.n = 0, jitters.n = 3, trace = TRUE, 
    eps = 0.1, prob = 0.4, scale.min = 1e-12, maxeval = Inf) 
{
    dcon <- c(eps = eps, prob = prob, scale.min = scale.min)
    icon <- c(births = births, random.n = random.n, jitters.n = jitters.n, 
        trace = trace, maxeval = maxeval)
    list(icontrol = icon, dcontrol = dcon)
}
