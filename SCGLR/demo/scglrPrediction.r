library(SCGLR)

# load sample data
data(genus)

# get variable names from dataset
n <- names(genus)
ny <- n[grep("^gen",n)]    # Y <- names that begins with "gen"
nx <- n[-grep("^gen",n)]   # X <- remaining names

# remove "geology" and "surface" from nx
# as surface is offset and we want to use geology as additional covariate
nx <-nx[!nx%in%c("geology","surface")]

# build multivariate formula
# we also add "lat*lon" as computed covariate
form <- multivariateFormula(ny,c(nx,"I(lat*lon)"),c("geology"))

# split genus dataset
sub <- sample(1:nrow(genus),100,replace=FALSE)
sub_fit <- (1:nrow(genus))[-sub]

# define family 
fam <- rep("poisson",length(ny))

# fit the model
genus.scglr <- scglr(formula=form, data=genus, family=fam, K=4, offset=genus$surface, subset=sub_fit)

#xnew, the design matrix associated to sub-sample used for prediction
# note rhs parameter is introduced to take into account that the covariate part of
# the formula is composed of two differents sets
xnew <- model.matrix(form, data=genus[sub,], rhs=1:2)[,-1]

# prediction based on the scglr approch
pred.scglr <- multivariatePredictGlm(xnew,family=fam,beta=genus.scglr$beta,offset=genus$surface[sub])
cor.scglr <-diag(cor(pred.scglr,genus[sub,ny])) 
plot(cor.scglr, col="red",ylim=c(-1,1))

# prediction based on classical poisson glm
genus.glm <- multivariateGlm(formula=form, data=genus, family=fam, offset=genus$surface, subset=sub_fit)
coefs <- sapply(genus.glm,coef)

# note rhs parameter is introduced to take into account that the covariate part of
# the formula is composed of two differents sets
pred.glm <- multivariatePredictGlm(xnew,family=fam,beta=coefs,offset=genus$surface[sub])
cor.glm <- diag(cor(pred.glm,genus[sub,ny]))
points(cor.glm, col="blue")
