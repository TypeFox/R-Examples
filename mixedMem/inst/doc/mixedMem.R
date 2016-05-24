## ----echo = FALSE, warning=FALSE-----------------------------------------
library(xtable, quietly = TRUE)

## ------------------------------------------------------------------------
library(mixedMem)
data(ANES)
# Dimensions of the data set: 279 individuals with 19 responses each
dim(ANES)
# The 19 variables and their categories
# The specific statements for each variable can be found using help(ANES)
# Variables titled EQ are about Equality
# Variables titled IND are about Econonic Individualism
# Variables titled ENT are about Free Enterprise
colnames(ANES)
# Distribution of responses
table(unlist(ANES))

## ------------------------------------------------------------------------
# Sample Size
Total <- 279
# Number of variables
J <- 19 
# we only have one replicate for each of the variables
Rj <- rep(1, J)
# Nijr indicates the number of ranking levels for each variable.
# Since all our data is multinomial it should be an array of all 1s
Nijr <- array(1, dim = c(Total, J, max(Rj)))
# Number of sub-populations
K <- 3
# There are 3 choices for each of the variables ranging from 0 to 2.
Vj <- rep(3, J)
# we initialize alpha to .2
alpha <- rep(.2, K)
# All variables are multinomial
dist <- rep("multinomial", J)
# obs are the observed responses. it is a 4-d array indexed by i,j,r,n
# note that obs ranges from 0 to 2 for each response
obs <- array(0, dim = c(Total, J, max(Rj), max(Nijr)))
obs[, , 1, 1] <- as.matrix(ANES)

# Initialize theta randomly with Dirichlet distributions
set.seed(123)
theta <- array(0, dim = c(J, K, max(Vj)))
for (j in 1:J) {
    theta[j, , ] <- gtools::rdirichlet(K, rep(.8, Vj[j]))
}

# Create the mixedMemModel
# This object encodes the initialization points for the variational EM algorithim
# and also encodes the observed parameters and responses
initial <- mixedMemModel(Total = Total, J = J, Rj = Rj,
                         Nijr = Nijr, K = K, Vj = Vj, alpha = alpha,
                         theta = theta, dist = dist, obs = obs)

## ------------------------------------------------------------------------
computeELBO(initial)
st = proc.time()
#printStatus 1 indicates that status updates will be printed
# printMod 25 indicates that status updates will only be printed ever 25th step
out <- mmVarFit(initial, printStatus = 1, printMod = 25)
end <- proc.time()
computeELBO(out)
time <- end - st

## ------------------------------------------------------------------------
summary(out)

## ------------------------------------------------------------------------
#read in GM estimates
# Note that Gross and Manrique-Vallier do not report estimates for the third group, so they have been set to 0 in the data
data(gmv_theta)
#
optimal.perm <- findLabels(out, gmv_theta)
# display the permutation as well as the weighted squared error loss
optimal.perm
# save object with permuted labels
out.permute <- permuteLabels(out, optimal.perm$perm)

## ------------------------------------------------------------------------
# calculate the expected value
model.exp <- matrix(0, nrow = J, ncol = max(Vj))
for(j in 1:J){
  model.exp[j, ] <- (out.permute$alpha / sum(out.permute$alpha)) %*% out.permute$theta[j, , ]
}
model.exp <- model.exp*out.permute$Total

# aggregate the observed counts
obs.counts <- t(apply(obs[, , 1, 1] + 1, MARGIN = 2, tabulate))
chi.squared <- sum( (obs.counts - model.exp)^2 / model.exp)

## ----results = 'asis', echo = FALSE--------------------------------------
tab <- cbind(model.exp[, 1], obs.counts[, 1],
             model.exp[, 2], obs.counts[, 2],
             model.exp[, 3], obs.counts[, 3])
rownames(tab) <- colnames(ANES)
colnames(tab) <- c("Exp Agree", "Obs Agree",  "Exp Can't Decide", "Obs Can't Decide", "Exp Disagree", "Obs Disagree")
xtable(tab, digits = 1, caption = "Expected responses based on fitted parameters and Observed responses")


## ------------------------------------------------------------------------
# The chi squared statistic
chi.squared
# p-value
pchisq(chi.squared, df = (J-1) * max(Vj))

## ------------------------------------------------------------------------
# The estimated quantities can be exctracted from the output model
names(out)
out$alpha

## ----theta, fig.height = 7, fig.width = 7, fig.align = 'center', fig.cap = "Individuals were asked to indicate their level of agreement with 19 opinion based statements. The fitted multinomial response probabilities for each ideology bloc to each of the first 10 statements are displayed. On the horizontal axis, 0 indicates agree, 1 indicates can't decide, and 2 indicates disagree. The estimates from our variational analysis are denoted by the dots, the Gross and Manrique-Vallier results using MCMC are shown by the red X's. Since they do not report estimates for Group 3, the Can't Decide group, those estimates are shown as 0's. The estimates for the first 10 questions are shown."----
plot(out.permute, type = "theta", compare = gmv_theta, varNames = colnames(ANES),
         groupNames = c("Conservative", "Liberal", "Undecided"), nrow = 10, indices = c(1:10))

## ----consprop, fig.width = 4, fig.height = 4, fig.align='center', fig.cap = "Propensity to agree with each opinion-based statement for the conservative bloc"----
pop1VarOrder <- colnames(ANES)[order(out.permute$theta[, 1, 1], decreasing = T)]
pop1VarAgree <- sort(out.permute$theta[, 1, 1], decreasing = T)
barplot(height = pop1VarAgree, names.arg = pop1VarOrder,
        main = "Propensity to Agree",
        cex.names = .7, las = 2, xlab = "Value Statements",
        ylab = expression(paste(theta["j,1,1"])),
        col = ifelse(pop1VarAgree > .5, "forestgreen", "darkred"))

## ----libprop, fig.width = 4, fig.height = 4,fig.align='center',  fig.cap = "Propensity to agree with each opinion-based statement for the liberal bloc"----
pop2VarOrder <- colnames(ANES)[order(out.permute$theta[, 2, 1], decreasing = T)]
pop2VarAgree <- sort(out.permute$theta[, 2, 1], decreasing = T)
barplot(height = pop2VarAgree,
        names.arg = pop2VarOrder, main = "Propensity to Agree",
        cex.names = .7, las = 2, xlab = "Value Statements",
        ylab = expression(paste(theta["j,2,1"])),
        col = ifelse(pop2VarAgree > .5, "forestgreen", "darkred"))

## ------------------------------------------------------------------------
# Point estimates for lambda
lambda.point <- out.permute$phi / rowSums(out.permute$phi)
# number of individuals which exhibit more than .3 degree of membership
# in the undecided group
sum(lambda.point[, 3] >= .3)
# number of can not decide responses from those with high membership in undecided group
sum(ANES[which(lambda.point[, 3] >= .3),] == 1)

## ------------------------------------------------------------------------
relativeFrequency <- out.permute$alpha / sum(out.permute$alpha)

## ----echo = FALSE, results = 'asis'--------------------------------------
table.alpha = rbind(out.permute$alpha, relativeFrequency)
rownames(table.alpha) = c("Estimated Alpha", "Estimate Relative Frequency")
colnames(table.alpha) = c("Conservatives", "Liberals", "Undecided")
xtable(table.alpha, caption = "Variational Estimates of Alpha", digits = 3)

## ----polarizing, fig.width = 4, fig.height = 4, fig.align = 'center', fig.cap = "Using Hellinger distance indicates that the most polarizing statements are IND 5, IND 3, and EQ 7. These were all issues involving opportunity for advancement"----
hellingerDist <- (1/sqrt(2)) * sqrt(rowSums( (sqrt(out.permute$theta[, 1, ]) 
                                          - sqrt(out.permute$theta[, 2, ]))^2))
barplot(sort(hellingerDist, decreasing = T),
        names.arg = colnames(ANES)[order(hellingerDist, decreasing = T)],
        main = "Hellinger Distance",
        cex.names = .7, las = 2, ylab = "Hellinger Distance",
        ylim = c(0, 1))
mtext("Between Conservatives and Liberals")
colnames(ANES)[order(hellingerDist, decreasing = T)][1:3]

## ------------------------------------------------------------------------
# number of individuals with at least 40% membership in
# both conservative and liberal blocs
sum( (lambda.point[, 1] > .4) & (lambda.point[, 2] > .4))

## ----mems, fig.width = 4, fig.height = 4, fig.align = 'center', fig.cap = "Estimated memberships for first 42 individuals"----
plot(out.permute, type = "membership", indices = c(1:20), nrow = 5, ncol = 4, main = "Estimated Memberships",
     groupNames = c("Cons", "Lib", "Cant Decide"), fitNames = c("mixedMem"))

## ----posteriorMem, fig.cap = "We observe a relatively high rate of intra-individual mixing", fig.width = 5, fig.align = 'center', fig.height = 5----
index <- order(lambda.point[, 1])

# variance of posterior membership in conservative bloc
var.Mem <- out.permute$phi[, 1] * (rowSums(out.permute$phi) -
                                 out.permute$phi[, 1]) /
  (rowSums(out.permute$phi)^2 * (rowSums(out.permute$phi) + 1))
# plot posterior means
plot(sort(lambda.point[, 1]), pch = 19,
     main = "Posterior Membership in Conservative Bloc",
     ylab = "Posterior Membership", xlab = "Individual",
     cex = .8, ylim = c(0, 1))

# marginal distirbution of Dirichlet, is Beta distribution, so we can get posterior CI
# plot posterior 90% CI
ci_up <- qbeta(.975, out.permute$phi[index, 1], rowSums(out.permute$phi[index, c(2:3)]))
ci_low <- qbeta(.025, out.permute$phi[index, 1], rowSums(out.permute$phi[index, c(2:3)]))
segments(x0 = c(1:out$Total), y0 = ci_up, y1 = ci_low, col = "gray", lwd = .3, lty = 1)
legend("bottomright", legend = c("posterior mean", "95% CI"), pch = c(19, NA),
       lty = c(NA, 1), col = c("black", "gray"))

