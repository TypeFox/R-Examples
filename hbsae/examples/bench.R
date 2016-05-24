d <- generateFakeData()

# compute small area estimates
sae <- fSAE(y0 ~ x + area2, data=d$sam, area="area", popdata=d$Xpop)

# calibrate to overall population total
sae.c <- bench(sae, rhs=sum(d$mY0*sae$Narea))
plot(sae, sae.c)
