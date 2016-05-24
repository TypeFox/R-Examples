over2 <-
function(para)
{
	alp1=para[[1]][1]
	alp2=para[[1]][2]
	mu1=para[[2]][1]
	mu2=para[[2]][2]
	sig1=para[[3]][1]
	sig2=para[[3]][2]
	part1=over1(alp1,mu1,sig1,alp2,mu2,sig2)
	part2=over1(alp2,mu2,sig2,alp1,mu1,sig1)
	(part1+part2)/2
}
