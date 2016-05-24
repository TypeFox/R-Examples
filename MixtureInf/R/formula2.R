formula2 <-
function(para,n)
{
	overlap=over2(para)
	tover=log(overlap/(1-overlap))
	an=0.35*exp(-1.859    -0.577*tover -60.453/n  )/(1+exp(-1.859    -0.577*tover -60.453/n))
	an
}
