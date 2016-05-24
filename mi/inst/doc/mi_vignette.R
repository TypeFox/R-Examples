## ----step0-----------------------------------------------------
options(width = 65)
suppressMessages(library(mi))
data(nlsyV, package = "mi")

## ----step1-----------------------------------------------------
mdf <- missing_data.frame(nlsyV)

## ----step1.5---------------------------------------------------
show(mdf) # momrace is guessed to be ordered

## ----, step2---------------------------------------------------
mdf <- change(mdf, y = c("income", "momrace"), what = "type", 
                     to = c("non", "un"))
show(mdf)

## ----, step3---------------------------------------------------
summary(mdf)
image(mdf)
hist(mdf)

## ----, step4---------------------------------------------------
rm(nlsyV)       # good to remove large unnecessary objects to save RAM
options(mc.cores = 2)
imputations <- mi(mdf, n.iter = 30, n.chains = 4, max.minutes = 20)
show(imputations)

## ----, step5A--------------------------------------------------
round(mipply(imputations, mean, to.matrix = TRUE), 3)
Rhats(imputations)

## ----, step5B--------------------------------------------------
imputations <- mi(imputations, n.iter = 5)

## ----, step6---------------------------------------------------
plot(imputations)
plot(imputations, y = c("ppvtr.36", "momrace"))
hist(imputations)
image(imputations)
summary(imputations)

## ----, step7---------------------------------------------------
analysis <- pool(ppvtr.36 ~ first + b.marr + income + momage + momed + momrace, 
                 data = imputations, m = 5)
display(analysis)

## ----, step8---------------------------------------------------
dfs <- complete(imputations, m = 2)

