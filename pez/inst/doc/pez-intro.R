## ----include=FALSE--------------------
require(pez)
options(width=40)

## ----tidy=TRUE, size="small"----------
library(pez)
data(laja)
data <- comparative.comm(invert.tree, river.sites, invert.traits, river.env)

## ----tidy=TRUE, size="small"----------
site.subset <- data[1:5,]
spp.subset <- data[,1:3]

## ----tidy=TRUE, size="small"----------
species(data)[1:2]
species(data)[1:2] <- c("new", "names")
sites(data)[1:2] <- c("newer", "names")
data <- data[, colSums(data$comm) > 5]
traits(data)$new.trait <- rep("nonsense", nrow(traits(data)))
traits(data)$new.trait <- NULL

## ----tidy=TRUE, size="small"----------
shape.output <- pez.shape(data)
dim(shape.output)
shape.output[1:3,1:3]

## ----tidy=TRUE, size="small"----------
sqrt <- pez.shape(data, sqrt.phy=TRUE)
traits <- pez.shape(data, traitgram=1) #traits alone
traits <- pez.shape(data, traitgram=c(0,0.5))#phylogeny and both
traits <- pez.shape(data, ext.dist=as.dist(cophenetic(phy(data))))

## ----tidy=TRUE, fig.width=6, fig.height=4, size="small"----
dist <- pez.dissimilarity(data, "phylosor")
plot(hclust(dist$phylosor))

## ----tidy=TRUE, size="small"----------
metrics <- generic.metrics(data, c(.mpd,.pse,.ses.mpd))
#null.comparisons <- generic.null(data, c(.mpd,.pse))
metrics <- generic.metrics(data, c(.mpd,.mntd), dist=as.dist(cophenetic(phy(data))))

## ----tidy=TRUE, size="small"----------
phy <- eco.phy.regression(data, permute=10)
trait <- eco.trait.regression(data, permute=10, method="quantile", tau=c(0.25,0.5,0.7))
trait <- eco.trait.regression(data, altogether=FALSE)

## ----warning=FALSE, tidy=TRUE, size="small"----
model <- fingerprint.regression(data, eco.permute=10)

## ----traitgram, warning=FALSE, fig.width=4, fig.height=4.5, dev.args=list(pointsize=8)----
assemblage <- c("Nerophilus", "Hydroptila", "Psorophora",
                "Simuliidae", "Psychodidae", "Ceratopogon",
                "Nectopsyche", "Pedomoecus", "Ceratopsyche")
dataAssemblage <- data[, species(data) %in% assemblage]
traitgram.cc(dataAssemblage, "length")

## ----FPDist, size="small"-------------
fpd.data <- funct.phylo.dist(data, phyloWeight = 0.5, p = 2)

## ----sesFPD, size="small"-------------
ses.mfpd.data <- .ses.mpd(data, dist=fpd.data)
head(ses.mfpd.data)[,c("ntaxa", "mpd.obs", "mpd.obs.p")]

## ----pglmmSim, tidy=TRUE, size="small"----
# Basic parameters
nspp <- 15; nsite <- 10
# Fixed effects
beta0 <- beta1 <- 0
# Random effects' magnitudes
sd.B0 <- sd.B1 <- 1

# Generate environmental site variable
X <- matrix(1:nsite, nrow=1, ncol=nsite)
X <- (X - mean(X))/sd(X)

# Simulate phylogeny
phy <- compute.brlen(rtree(nspp), method = "Grafen", power = 0.5)
# Standardise phy. covariance matrix
Vphy <- vcv(phy); Vphy <- Vphy/(det(Vphy)^(1/nspp))

# Generate phylogenetic signal in parameters
iD <- t(chol(Vphy))
b0 <- beta0 + iD %*% rnorm(nspp, sd = sd.B0)
b1 <- beta1 + iD %*% rnorm(nspp, sd = sd.B1)

#Simulate presences
y <- matrix(outer(b0, array(1, dim = c(1, nsite))), nrow = nspp, ncol = nsite) + matrix(outer(b1, X), nrow = nspp, ncol = nsite)
e <- rnorm(nspp * nsite, sd = 0)
y <- y + matrix(e, nrow = nspp, ncol = nsite)
y <- matrix(y, nrow = nspp * nsite, ncol = 1)    
Y <- rbinom(n = length(y), size = 1, prob = exp(y)/(1 + exp(y)))
Y <- matrix(Y, nrow = nspp, ncol = nsite)

#Neat up the data to show structure
rownames(Y) <- 1:nspp; colnames(Y) <- 1:nsite

## ----pglmmModel, tidy=TRUE, size="small"----
# Transform data into 'long' format
# - Occurrence data
YY <- matrix(Y, nrow = nspp * nsite, ncol = 1)
# - Environmental (site) data
XX <- matrix(kronecker(X, matrix(1, nrow = nspp, ncol = 1)), nrow = nspp * nsite, ncol = 1)
site <- matrix(kronecker(1:nsite, matrix(1, nrow = nspp, ncol = 1)), nrow = nspp * nsite, ncol = 1)
sp <- matrix(kronecker(matrix(1, nrow = nsite, ncol = 1), 1:nspp), nrow = nspp * nsite, ncol = 1)
# - Make data.frame with all data
dat <- data.frame(Y = YY, X = XX, site = as.factor(site), sp = as.factor(sp))

# Setup random effects
# - 1: random intercept - species independent
re.1 <- list(1, sp = dat$sp, covar = diag(nspp))
# - 2: random intercept - species phylogenetically covary
re.2 <- list(1, sp = dat$sp, covar = Vphy)
# - 3: random slope - species independent
re.3 <- list(dat$X, sp = dat$sp, covar = diag(nspp))
# - 4: random intercept - species covary
re.4 <- list(dat$X, sp = dat$sp, covar = Vphy)   
# (Random effect for site)
re.site <- list(1, site = dat$site, covar = diag(nsite))

# Fit model!
model <- communityPGLMM(Y ~ X, data = dat, family = "binomial", sp = dat$sp, site = dat$site, random.effects = list(re.1, re.2, re.3, re.4), REML = TRUE, verbose = FALSE)

