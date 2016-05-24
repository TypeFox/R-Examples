## ----inst1,eval=FALSE,tidy=FALSE-----------------------------------------
#  install.packages("MPTinR")

## ----load,eval=TRUE,tidy=FALSE-------------------------------------------
library("MPTinR")

## ----mod1,eval=FALSE,tidy=FALSE------------------------------------------
#  # Tree for old items: First 'yes', then 'no'
#  Do + (1 - Do) * g
#  (1-Do)*(1-g)
#  
#  #Tree for new items: First 'yes', then 'no'
#  (1-Dn) * g
#  Dn + (1-Dn) * (1 - g)

## ----rest1,eval=FALSE,tidy=FALSE-----------------------------------------
#  " # quotes need to be removed
#  W1 < W2 < W3
#  X4 = X3
#  Y1 = Y3 = 0.5
#  Z = 0 #Restrictions may also contain comments.
#  "

## ----mod_check-----------------------------------------------------------

mod_2htm_1 <- "
# Tree for old items: First 'yes', then 'no'
Do + (1 - Do) * g
(1-Do)*(1-g)

#Tree for new items: First 'yes', then 'no'
(1-Dn) * g
Dn + (1-Dn) * (1 - g)
"
check.mpt(textConnection(mod_2htm_1))


## ------------------------------------------------------------------------
d.broeder.agg <- c(145, 95, 170, 1990, 402, 198, 211, 1589, 868, 332, 
                   275, 925, 1490, 310, 194, 406, 1861, 299, 94, 146)

## ----mod_check_2,eval=TRUE,tidy=FALSE------------------------------------
mod_2htm_2 <- "
# Tree for old items (10%): First 'yes', then 'no' 
Do + (1 - Do) * g1
(1-Do)*(1-g1)

#Tree for new items (90%): First 'yes', then 'no'
(1-Dn) * g1
Dn + (1-Dn) * (1 - g1)
    
# Tree for old items (25%): First 'yes', then 'no'
Do + (1 - Do) * g2
(1-Do)*(1-g2)

#Tree for new items  (75%): First 'yes', then 'no'
(1-Dn) * g2
Dn + (1-Dn) * (1 - g2)

# Tree for old items (50%): First 'yes', then 'no'
Do + (1 - Do) * g3
(1-Do)*(1-g3)

#Tree for new items  (50%): First 'yes', then 'no'
(1-Dn) * g3
Dn + (1-Dn) * (1 - g3)

# Tree for old items (75%): First 'yes', then 'no'
Do + (1 - Do) * g4
(1-Do)*(1-g4)

#Tree for new items  (25%): First 'yes', then 'no'
(1-Dn) * g4
Dn + (1-Dn) * (1 - g4)

# Tree for old items (90%): First 'yes', then 'no'
Do + (1 - Do) * g5
(1-Do)*(1-g5)

#Tree for new items  (10%): First 'yes', then 'no'
(1-Dn) * g5
Dn + (1-Dn) * (1 - g5)
"



## ----check3--------------------------------------------------------------
check.mpt(textConnection(mod_2htm_2))


## ----dat1----------------------------------------------------------------
data("d.broeder")
head(d.broeder)


## ----fit_br_1------------------------------------------------------------
br.2htm <- fit.mpt(d.broeder, textConnection(mod_2htm_2), fia = 25000)

## ----fit_br_1b-----------------------------------------------------------
br.2htm.ineq <- fit.mpt(d.broeder, textConnection(mod_2htm_2), 
                        list("g1 < g2 < g3 < g4 < g5"), fia = 25000)

## ------------------------------------------------------------------------
d.broeder[c(2, 6, 7),]

## ------------------------------------------------------------------------
str(br.2htm, 1)

## ------------------------------------------------------------------------
br.2htm$model.info

## ----br_res1-------------------------------------------------------------
br.2htm[["goodness.of.fit"]][["aggregated"]]

br.2htm[["parameters"]][["aggregated"]]

## ------------------------------------------------------------------------
br.2htm[["parameters"]][["aggregated"]][,"estimates"]

br.2htm.ineq[["parameters"]][["aggregated"]][,"estimates"]

## ------------------------------------------------------------------------
br.2htm[["goodness.of.fit"]][2:3]

br.2htm.ineq[["goodness.of.fit"]][2:3]

## ------------------------------------------------------------------------
br.2htm$information.criteria[2:3]

br.2htm.ineq$information.criteria[2:3]

## ------------------------------------------------------------------------
select.mpt(list(br.2htm, br.2htm.ineq), output = "full")

## ------------------------------------------------------------------------
select.mpt(list(br.2htm, br.2htm.ineq), output = "full", dataset = (1:40)[-c(2, 7)])

## ----fit_br_2------------------------------------------------------------
br.2htmr <- fit.mpt(d.broeder, textConnection(mod_2htm_2), list("Do = Dn"), fia = 25000)

## ----fit_br_2b-----------------------------------------------------------
br.2htmr.ineq <- fit.mpt(d.broeder, textConnection(mod_2htm_2), 
                        list("g1 < g2 < g3 < g4 < g5", "Do = Dn"), fia = 25000)

## ------------------------------------------------------------------------
select.mpt(list(br.2htm, br.2htm.ineq, br.2htmr, br.2htmr.ineq), output = "full")

## ------------------------------------------------------------------------
str(br.2htmr$parameters$individual)

## ------------------------------------------------------------------------

apply(br.2htm$parameters$individual[,1,], 1, mean)
apply(br.2htm$parameters$individual[,1,], 1, median)

apply(br.2htmr$parameters$individual[,1,], 1, mean)
apply(br.2htmr$parameters$individual[,1,], 1, median)


## ------------------------------------------------------------------------
br.2htm$parameters$mean

## ----fit_br_3------------------------------------------------------------
br.1htm <- fit.mpt(d.broeder, textConnection(mod_2htm_2), list("Dn = 0"), fia = 25000)

## ----fit_br_3b-----------------------------------------------------------
br.1htm.ineq <- fit.mpt(d.broeder, textConnection(mod_2htm_2), 
                        list("g1 < g2 < g3 < g4 < g5", "Dn = 0"), fia = 25000)

## ------------------------------------------------------------------------
select.mpt(list(br.2htm, br.2htm.ineq, br.2htmr, br.2htmr.ineq, br.1htm, br.1htm.ineq), 
           output = "full")[,1:16]

## ------------------------------------------------------------------------
br.2htm.2 <- fit.mpt(colSums(d.broeder), textConnection(mod_2htm_2))

## ------------------------------------------------------------------------
t(br.2htm.2[["parameters"]])

bs.data <- gen.data(br.2htm.2[["parameters"]][,1], 200, 
                    textConnection(mod_2htm_2), data = colSums(d.broeder))

br.2htm.bs <- fit.mpt(bs.data, textConnection(mod_2htm_2), fit.aggregated = FALSE)

## ------------------------------------------------------------------------
apply(br.2htm.bs[["parameters"]][["individual"]][,1,],
      1, quantile, probs = c(0.025, 0.975))

## ----misfit-plots, fig.width=8.5, fig.height=5, out.width='0.95\\textwidth'----
br.2htm.2 <- fit.mpt(colSums(d.broeder), textConnection(mod_2htm_2), show.messages = FALSE)

axis.labels <- c("10%-O", "90%-N", "25%-O", "75%-N", "50%-O", "50%-N", 
                 "75%-O", "25%-N", "90%-O", "10%-N")
prediction.plot(br.2htm.2, textConnection(mod_2htm_2), axis.labels = axis.labels,
                args.plot = list(main = "Absolute deviations (frequency scale)"))


