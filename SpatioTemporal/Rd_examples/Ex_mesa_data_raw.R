##load the data
data(mesa.data.raw)

##extract matrix of observations (missing marked by NA)
obs.mat <- mesa.data.raw$obs
head(obs.mat)

##optionally observations can be given as a data.frame
obs <- data.frame(obs=c(obs.mat),
                  date=rep(rownames(obs.mat), dim(obs.mat)[2]),
                  ID=rep(colnames(obs.mat), each=dim(obs.mat)[1]))
##force date-format
obs$date <- as.Date(obs$date)

##drop unobserved
obs <- obs[!is.na(obs$obs),,drop=FALSE]

##create a 3D-array for the spatio-temporal covariate
ST <- array(mesa.data.raw$lax.conc.1500, dim =
            c(dim(mesa.data.raw$lax.conc.1500),1))
dimnames(ST) <- list(rownames(mesa.data.raw$lax.conc),
                     colnames(mesa.data.raw$lax.conc),
                     "lax.conc.1500")
##or use a list of matrices
ST.list <- list(lax.conc.1500=mesa.data.raw$lax.conc.1500)

###########################
## create STdata object ##
###########################
##Create the data-object
mesa.data <- createSTdata(obs.mat, mesa.data.raw$X, n.basis=2,
                          SpatioTemporal=ST)
mesa.data.2 <- createSTdata(obs, mesa.data.raw$X, n.basis=2,
                            SpatioTemporal=ST.list)

##This should yield equal structures,
##which are also the same as data(mesa.data)
all.equal(mesa.data, mesa.data.2)

###########################
## create STmodel object ##
###########################
##define land-use covariates, for intercept and trends
LUR <- list(~log10.m.to.a1+s2000.pop.div.10000+km.to.coast,
  ~km.to.coast, ~km.to.coast)
##and covariance model
cov.beta <- list(covf="exp", nugget=FALSE)
cov.nu <- list(covf="exp", nugget=~type, random.effect=FALSE)
##which locations to use
locations <- list(coords=c("x","y"), long.lat=c("long","lat"), others="type")
##create object
mesa.model <- createSTmodel(mesa.data, LUR=LUR, ST="lax.conc.1500",
                            cov.beta=cov.beta, cov.nu=cov.nu,
                            locations=locations)

##This should be the same as the data in data(mesa.model)
