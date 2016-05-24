### R code from vignette source 'raschtree.Rnw'
### Encoding: ASCII

###################################################
### code chunk number 1: prep_install (eval = FALSE)
###################################################
## install.packages("psychotree")


###################################################
### code chunk number 2: prep_library
###################################################
library("psychotree")


###################################################
### code chunk number 3: prep_data
###################################################
data("SPISA", package = "psychotree")


###################################################
### code chunk number 4: recode (eval = FALSE)
###################################################
## mydata$resp <- as.matrix(mydata[ , 1:5])


###################################################
### code chunk number 5: recode2 (eval = FALSE)
###################################################
## mydata <- mydata[ , -(1:5)]


###################################################
### code chunk number 6: delete_cases (eval = FALSE)
###################################################
## mydata <- subset(mydata, rowMeans(resp, na.rm = TRUE) > 0 &
##   rowMeans(resp, na.rm = TRUE) < 1)


###################################################
### code chunk number 7: fit_raschtree (eval = FALSE)
###################################################
## my_first_raschtree <- raschtree(spisa ~ age + gender +
##   semester + elite + spon, data = SPISA)


###################################################
### code chunk number 8: fit_raschtree_cache
###################################################
if(file.exists("raschtree-spisa.rda")) load("raschtree-spisa.rda") else {
my_first_raschtree <- raschtree(spisa ~ age + gender +
  semester + elite + spon, data = SPISA)
save(my_first_raschtree, file = "raschtree-spisa.rda")
}
file.remove("raschtree-spisa.rda")


###################################################
### code chunk number 9: plot_raschtree
###################################################
plot(my_first_raschtree)


###################################################
### code chunk number 10: plot_raschtree_col
###################################################
library("colorspace")
plot(my_first_raschtree, 
      col = rep(rainbow_hcl(5, c = 65, l = 65), each = 9))


###################################################
### code chunk number 11: coef_raschtree
###################################################
coef(my_first_raschtree, node = 4)


###################################################
### code chunk number 12: itempar_raschtree
###################################################
itempar(my_first_raschtree, node = 4)


