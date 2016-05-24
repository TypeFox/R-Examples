formula3 <-
function(para,n)
{
	overlap=over3(para)
	tover1=log(overlap[1]/(1-overlap[1]))
	tover2=log(overlap[2]/(1-overlap[2]))
	an2=-1.602-0.240*tover1-0.240*tover2-130.394/n
	an=0.35*exp(an2)/(1+exp(an2))
	an
}
