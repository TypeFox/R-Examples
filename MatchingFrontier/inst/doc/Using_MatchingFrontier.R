### R code from vignette source 'Using_MatchingFrontier.Rnw'

###################################################
### code chunk number 1: Using_MatchingFrontier.Rnw:198-200 (eval = FALSE)
###################################################
## library(devtools) 
## install_github('ChristopherLucas/MatchingFrontier')


###################################################
### code chunk number 2: Using_MatchingFrontier.Rnw:206-210
###################################################
cat("","$ curl -OL https://github.com/ChristopherLucas/MatchingFrontier/archive/master.zip\n",
"$ unzip master.zip\n",
"$ cd MatchingFrontier-master\n",
"$ R CMD INSTALL package\n")


###################################################
### code chunk number 3: Using_MatchingFrontier.Rnw:323-338
###################################################
# Load the package and the data
library(MatchingFrontier)
data('lalonde')

# Create a vector of column names to indicate which variables we 
# want to match on. We will match on everything except the treatment
# and the outcome.
match.on <- colnames(lalonde)[!(colnames(lalonde) %in% c('re78', 'treat'))]
match.on # Print variables in match.on
# Make the mahalanobis frontier
mahal.frontier <- makeFrontier(dataset = lalonde, 
                            treatment = 'treat', 
                            outcome = 're78', 
                            match.on = match.on)
mahal.frontier


###################################################
### code chunk number 4: Using_MatchingFrontier.Rnw:358-367
###################################################
# Make the L1 frontier
L1.frontier <- makeFrontier(dataset = lalonde, 
                            treatment = 'treat', 
                            outcome = 're78', 
                            match.on = match.on,
                            QOI = 'SATT',
                            metric = 'L1',
                            ratio = 'fixed')
L1.frontier


###################################################
### code chunk number 5: Using_MatchingFrontier.Rnw:381-408
###################################################
# Set base form
my.form <- as.formula(re78 ~ treat + age + black + education + hispanic +
                          married + nodegree + re74 + re75)

# Estimate effects for the mahalanobis frontier
mahal.estimates <- estimateEffects(mahal.frontier, 
                                   're78 ~ treat',
                                   mod.dependence.formula = my.form,
                                   continuous.vars = c('age', 
                                       'education', 
                                       're74', 
                                       're75'),
                                   prop.estimated = .1,
                                   means.as.cutpoints = TRUE
                                   )

# Estimate effects for the L1 frontier
L1.estimates <- estimateEffects(L1.frontier, 
                                're78 ~ treat',
                                mod.dependence.formula = my.form,
                                continuous.vars = c('age', 
                                    'education', 
                                    're74', 
                                    're75'),
                                prop.estimated = .1,
                                means.as.cutpoints = TRUE
                                )


###################################################
### code chunk number 6: l1_frontier_plain
###################################################
# Plot frontier
plotFrontier(L1.frontier)


###################################################
### code chunk number 7: l1_frontier_pretty
###################################################
# Plot frontier
plotFrontier(L1.frontier,
             cex.lab = 1.4,
             cex.axis = 1.4,
             type = 'l',
             panel.first = 
                grid(NULL, 
                     NULL, 
                     lwd = 2)
             )


###################################################
### code chunk number 8: l1_frontier_est
###################################################
# Plot estimates
plotEstimates(L1.estimates, 
              ylim = 
                  c(-10000, 
                    3000),
              cex.lab = 1.4,
              cex.axis = 1.4,
              panel.first = 
                 grid(NULL,
                      NULL,
                      lwd = 2,
                      )
              )


###################################################
### code chunk number 9: l1_means
###################################################
# Plot estimates
plotMeans(L1.frontier)


###################################################
### code chunk number 10: l1_parplot
###################################################
# Make parallel plot
parallelPlot(L1.frontier,
             N = 200,
             variables = c('age',
             're74',
             're75',
             'black'),
             treated.col = 'gray',
             control.col = 'blue'
             )


###################################################
### code chunk number 11: Using_MatchingFrontier.Rnw:579-581
###################################################
n <- 200 # Identify the point from which to select the data
matched.data <- generateDataset(L1.frontier, N = n)


###################################################
### code chunk number 12: tab1
###################################################
library(stargazer)

mod1 <- lm(re78 ~ treat, data = matched.data)
mod2 <- lm(paste('re78 ~ treat +', paste(match.on, collapse = ' + ')), 
           data = matched.data)

stargazer(mod1, mod2, dep.var.labels = c('re78', 're78'))


