set.seed(10)

############ test whether predict.geolmStd results are correct

# check accuracy of method (compare to geoR results)
results.diff.uk <- matrix(99, ncol = 2, nrow = 16)

# generate possible parameter values for covariance models
sigmasq <- 1/rgamma(8, 2)
phi <- 1/rgamma(8, 2)
error.var <- 1/rgamma(8, 3)
micro <- c(rep(0, 4), runif(4))
kappa <- runif(8, .25, 3)
cm <- sample(c("exponential", "spherical", "matern"), 8, replace = TRUE)
n = rpois(8, 150)
np = rpois(8, 300)
ntimes = 1

# store universal kriging results
uk.results = matrix(99, nrow = 8, ncol = 3)

for(i in 1:8)
{
  # generate data
	fields <- grf(100, cov.pars = c(sigmasq[i], phi[i]), cov.model = cm[[i]], 
		kappa = kappa[i], messages = FALSE)
	# extract response
	y <- fields$data
	# get observed coordinates
	coords = fields$coords
	# generate unsampled coordinates
	pcoords = matrix(runif(np[i] * 2), ncol = 2)
	# combine both to check properties at samples sites
	acoords = rbind(coords, pcoords)

	# create geoR kriging model
	modeli <- krige.control(type = "ok", trend.d = "1st", trend.l = "1st",
	cov.model = cm[i], cov.pars = c(sigmasq[i], phi[i]), kappa = kappa[i], 
	nugget = (error.var[i] + micro[i]), micro.scale = micro[i])

	# create data frames needed for gear 
	data = data.frame(x1 = coords[,1], x2 = coords[,2], y = y)
	newdata = data.frame(x1 = acoords[,1], x2 = acoords[,2])
	
	# decide whether signal/filtered or unfiltered model
	if((i%%2) == 0) # unfiltered
	{
		output = output.control(signal = FALSE, messages = FALSE)
		# create appropriate covariance model for gear
		cmod = cmod.std(model = cm[i], psill = sigmasq[i], r = phi[i], 
		par3 = kappa[i], evar = 0, fvar = (micro[i] + error.var[i]))
	}else #filtered
	{
		output = output.control(signal = TRUE, messages = FALSE)
		# create appropriate covariance model for gear
		cmod = cmod.std(model = cm[i], psill = sigmasq[i], 
		                r = phi[i], 
		              par3 = kappa[i], evar = error.var[i], 
		              fvar = micro[i])
	}
	# gear geostatistic model
	gearmod = geolm(y ~ x1 + x2, data = data,
	              coordnames = c("x1", "x2"),
	              cmod = cmod)
	georout = krige.conv(fields, loc = acoords, krige = modeli, output = output)

	gearout = predict.geolmStd(gearmod,  newdata = newdata, sp = FALSE)
	
	uk.results[i, 1] = max(abs(range(georout$predict - gearout$pred)))
	uk.results[i, 2] = max(abs(range(georout$krige.var - gearout$mspe)))
	uk.results[i, 3] = max(abs(range(georout$beta.est - gearmod$coeff)))
}

test_that("all predict.geolmStd uk calculations are correct", {
  expect_true(max(uk.results) < 1e-10)
})

ok.results = matrix(99, nrow = 8, ncol = 3)

for(i in 1:8)
{
	fields <- grf(100, cov.pars = c(sigmasq[i], phi[i]), cov.model = cm[[i]], 
		kappa = kappa[i], messages = FALSE)
	y <- fields$data
	coords = fields$coords
	pcoords = matrix(runif(np[i] * 2), ncol = 2)
	acoords = rbind(coords, pcoords)

	x = matrix(1, nrow(coords))
	newx = matrix(1, nrow(acoords))	

	modeli <- krige.control(type = "ok", trend.d = "cte", trend.l = "cte",
	cov.model = cm[i], cov.pars = c(sigmasq[i], phi[i]), kappa = kappa[i], 
	nugget = (error.var[i] + micro[i]), micro.scale = micro[i])

	# create df for geolm
	data = data.frame(x1 = coords[,1], x2 = coords[,2], y = y)
	newdata = data.frame(x1 = acoords[,1], x2 = acoords[,2])
	
	# decide whether signal model
	if((i%%2) == 0)
	{
	  output = output.control(signal = FALSE, messages = FALSE)
	  cmod = cmod.std(model = cm[i], psill = sigmasq[i], r = phi[i], 
	                  par3 = kappa[i], evar = 0, 
	                  fvar = (micro[i] + error.var[i]))
	}else
	{
	  output = output.control(signal = TRUE, messages = FALSE)
	  cmod = cmod.std(model = cm[i], psill = sigmasq[i], 
	                  r = phi[i], 
	                  par3 = kappa[i], evar = error.var[i], 
	                  fvar = micro[i])
	}
	
	georout = krige.conv(fields, loc = acoords, krige = modeli, output = output)

	# create geolm for gear package
	gearmod = geolm(y ~ 1, data = data,
	              coordnames = c("x1", "x2"),
	              cmod = cmod)
	
	gearout = predict.geolmStd(gearmod,  newdata = newdata, sp = FALSE)
	
	ok.results[i, 1] = max(abs(range(georout$predict - gearout$pred)))
	ok.results[i, 2] = max(abs(range(georout$krige.var - gearout$mspe)))
	ok.results[i, 3] = max(abs(range(georout$beta.est - gearmod$coeff)))
}

test_that("all predict.geolmStd ok calculations are correct", {
  expect_true(max(ok.results) < 1e-10)
})

sk.results = matrix(99, nrow = 8, ncol = 2)

for(i in 1:8)
{
	fields <- grf(100, cov.pars = c(sigmasq[i], phi[i]), cov.model = cm[[i]], 
		kappa = kappa[i], messages = FALSE)
	y <- fields$data
	coords = fields$coords
	pcoords = matrix(runif(np[i] * 2), ncol = 2)
	acoords = rbind(coords, pcoords)
	mus = rnorm(8, 0, sd = 25)

	modeli <- krige.control(type = "sk", trend.d = "cte", trend.l = "cte",
	cov.model = cm[i], cov.pars = c(sigmasq[i], phi[i]), kappa = kappa[i], 
	nugget = (error.var[i] + micro[i]), micro.scale = micro[i], beta = mus[i])

	# decide whether signal model
	if((i%%2) == 0)
	{
	  output = output.control(signal = FALSE, messages = FALSE)
	  cmod = cmod.std(model = cm[i], psill = sigmasq[i], r = phi[i], 
	                  par3 = kappa[i], evar = 0, 
	                  fvar = (micro[i] + error.var[i]))
	}else
	{
	  output = output.control(signal = TRUE, messages = FALSE)
	  cmod = cmod.std(model = cm[i], psill = sigmasq[i], 
	                  r = phi[i], 
	                  par3 = kappa[i], evar = error.var[i], 
	                  fvar = micro[i])
	}
	
	georout = krige.conv(fields, loc = acoords, krige = modeli, output = output)

	data = data.frame(x1 = coords[,1], x2 = coords[,2], y = y)
	newdata = data.frame(x1 = acoords[,1], x2 = acoords[,2])
	
	# create geolm for gear package
	gearmod = geolm(y ~ 0, data = data,
	                 coordnames = c("x1", "x2"),
	                 cmod = cmod, mu = mus[i])
	
	gearout = predict.geolmStd(gearmod,  newdata = newdata, sp = FALSE)

	sk.results[i, 1] = max(abs(range(georout$predict - gearout$pred)))
	sk.results[i, 2] = max(abs(range(georout$krige.var - gearout$mspe)))
}

test_that("all fit.std sk calculations are correct", {
  expect_true(max(sk.results) < 1e-10)
})