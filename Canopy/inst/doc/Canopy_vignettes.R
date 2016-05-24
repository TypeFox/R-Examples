### R code from vignette source 'Canopy_vignettes.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: Installation (eval = FALSE)
###################################################
## install.packages('ape')
## install.packages('fields')
## install.packages('Canopy_0.99.0.tar.gz', repos = NULL, type="source")


###################################################
### code chunk number 2: Input
###################################################
library(Canopy)
data("MDA231")

projectname = MDA231$projectname ## name of project
R = MDA231$R; R ## mutant allele read depth (for SNAs)
X = MDA231$X; X ## total depth (for SNAs)
WM = MDA231$WM; WM ## observed major copy number (for CNA regions)
Wm = MDA231$Wm; Wm ## observed minor copy number (for CNA regions)
epsilonM = MDA231$epsilonM ## standard deviation of WM, pre-fixed here
epsilonm = MDA231$epsilonm ## standard deviation of Wm, pre-fixed here
## Matrix C specifices whether CNA regions harbor specific CNAs 
## only needed if overlapping CNAs are observed
C = MDA231$C; C
Y = MDA231$Y; Y ## whether SNAs are affected by CNAs


###################################################
### code chunk number 3: Sampling1 (eval = FALSE)
###################################################
## K = 3:6 # number of subclones
## numchain = 20 # number of chains with random initiations
## sampchain = canopy.sample(R = R, X = X, WM = WM, Wm = Wm, epsilonM = epsilonM, 
##             epsilonm = epsilonm, C = C, Y = Y, K = K, numchain = numchain, 
##             simrun = 50000, writeskip = 200, projectname = projectname,
##             cell.line = TRUE, plot.likelihood = TRUE)
## save.image(file = paste(projectname, '_postmcmc_image.rda',sep=''),
##            compress = 'xz')


###################################################
### code chunk number 4: Sampling2
###################################################
data("MDA231_sampchain")
sampchain = MDA231_sampchain
k = 3
K = 3:6
sampchaink = MDA231_sampchain[[which(K==k)]]


###################################################
### code chunk number 5: Sampling3
###################################################
length(sampchain) ## number of subtree spaces (K=3:6)
length(sampchain[[which(K==4)]]) ## number of chains for subtree space with 4 subclones
length(sampchain[[which(K==4)]][[1]]) ## number of posterior trees in each chain


###################################################
### code chunk number 6: BIC
###################################################
burnin = 100
thin = 10
bic = canopy.BIC(sampchain = sampchain, projectname = projectname, K = K,
               numchain = numchain, burnin = burnin, thin = thin, pdf = FALSE)
optK = K[which.max(bic)]


###################################################
### code chunk number 7: fig1
###################################################
par(mfrow=c(1,2))
projectname = 'MDA231'
numchain = 20
clikelihood = matrix(nrow = numchain, ncol = length(sampchaink[[1]]), data = NA)
for(numi in 1:numchain){
  for(i in 1:ncol(clikelihood)){
    clikelihood[numi,i] = sampchaink[[numi]][[i]]$likelihood
  }
}
plot(1:ncol(clikelihood), clikelihood[1,], type='l', xlab = 'Iteration',
     ylab = 'Log-likelihood', col = 1, ylim = c(min(clikelihood), 
                                                max(clikelihood)))
for(numi in 2:numchain){
  points(1:ncol(clikelihood), clikelihood[numi,], type = 'l', col = numi)
}
title(main=paste('Posterior likelihood', k, 'clones', numchain,
            'chains'),cex=0.6)
plot(K, bic, xlab = 'Number of subclones', ylab = 'BIC', type = 'b', xaxt = "n")
axis(1, at = K)
abline(v = (3:6)[which.max(bic)], lty = 2)
title('BIC for model selection')


###################################################
### code chunk number 8: Post
###################################################
post = canopy.post(sampchain = sampchain, projectname = projectname, K = K,
                 numchain = numchain, burnin = burnin, thin = thin, optK = optK,
                 C = C, post.config.cutoff = 0.05)
samptreethin = post[[1]]   # list of all post-burnin and thinning trees
samptreethin.lik = post[[2]]   # likelihoods of trees in samptree
config = post[[3]] # configuration for each posterior tree
config.summary = post[[4]] # configuration summary
print(config.summary)
# first column: tree configuration
# second column: posterior configuration probability in the entire tree space
# third column: posterior configuration likelihood in the subtree space


###################################################
### code chunk number 9: Plot
###################################################
config.i = config.summary[which.max(config.summary[,3]),1]
cat('Configuration', config.i, 'has the highest posterior likelihood!\n')
output.tree = canopy.output(post, config.i, C)
canopy.plottree(output.tree, pdf = FALSE)


###################################################
### code chunk number 10: fig2
###################################################
canopy.plottree(output.tree, pdf = FALSE)


###################################################
### code chunk number 11: Try (eval = FALSE)
###################################################
## data(toy)
## projectname = 'toy'
## R = toy$R; X = toy$X; WM = toy$WM; Wm = toy$Wm
## epsilonM = toy$epsilonM; epsilonm = toy$epsilonm; Y = toy$Y
## 
## K = 3:6; numchain = 10
## sampchain = canopy.sample(R = R, X = X, WM = WM, Wm = Wm, epsilonM = epsilonM, 
##                           epsilonm = epsilonm, C = NULL, Y = Y, K = K, 
##                           numchain = numchain, simrun = 50000, writeskip = 200,
##                           projectname = projectname, cell.line = FALSE,
##                           plot.likelihood = TRUE)


###################################################
### code chunk number 12: fig3
###################################################
data(toy)
canopy.plottree(toy$besttree, pdf = FALSE)


###################################################
### code chunk number 13: sessionInfo
###################################################
toLatex(sessionInfo())


