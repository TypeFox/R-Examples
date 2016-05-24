"f.B" <-
function(t, mu = 0.341, tau = 0.247, c1 = 2.58)
{
# COMPUTES THE SURVIVAL FUNCTION B, FORMULA (6)
#
### 	pnorm((c1 - mu * t)/sqrt(t^2 * tau^2 + t)) - exp(2 * c1 * mu + 2 * c1^2 *
###		tau^2) * pnorm(( - c1 - 2 * c1 * t * tau^2 - mu * t)/sqrt(t^2 * 
###		tau^2 + t))
 	pnorm((c1 - mu * t)/sqrt(t^2 * tau^2 + t)) - exp(2 * c1 * mu + 2 * c1^2 *
		tau^2 + pnorm(( - c1 - 2 * c1 * t * tau^2 - mu * t)/sqrt(t^2 * 
		tau^2 + t), log.p = T) )
}

