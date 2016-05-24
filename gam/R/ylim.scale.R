"ylim.scale" <-
function(ylim, scale = 0.)
{
	scale2 <- diff(ylim)
	if(scale2 < scale)
		rep(mean(ylim), 2.) + ((ylim - mean(ylim)) * scale)/scale2
	else ylim
}
