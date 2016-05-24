## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(  #default code chunk options
    dev = "CairoPNG",      #nicer PNG figures
    fig.width = 6,
    fig.height = 4
)           
pander::panderOptions("table.split.table", Inf)     #don't split wide tables in output
pander::panderOptions("table.style", "rmarkdown")   #table style that's supported by github

## ----message=FALSE-------------------------------------------------------
library(dplyr)      #%>%
library(lsmeans)    #lsmeans
library(DescTools)  #EtaSq
library(car)        #sigmaHat
library(ARTool)     #art, artlm
library(ggplot2)    #ggplot, stat_..., geom_..., etc

## ------------------------------------------------------------------------
data(InteractionTestData)
df = InteractionTestData    #save some typing

## ------------------------------------------------------------------------
#we'll be doing type 3 tests, so we want sum-to-zero contrasts
options(contrasts = c("contr.sum", "contr.poly"))
m.linear = lm(Y ~ X1*X2, data=df)

## ------------------------------------------------------------------------
m.art = art(Y ~ X1*X2, data=df)

## ------------------------------------------------------------------------
m.art.anova = anova(m.art)
print(m.art.anova, verbose=TRUE)

## ------------------------------------------------------------------------
m.art.anova$eta.sq.part = with(m.art.anova, `Sum Sq`/(`Sum Sq` + `Sum Sq.res`))
m.art.anova

## ------------------------------------------------------------------------
EtaSq(m.linear, type=3)

## ------------------------------------------------------------------------
x2.contrasts = summary(pairs(lsmeans(m.linear, ~ X2)))

## ------------------------------------------------------------------------
x2.contrasts$d = x2.contrasts$estimate / sigmaHat(m.linear)
x2.contrasts

## ------------------------------------------------------------------------
m.art.x2 = artlm(m.art, "X2")
x2.contrasts.art = summary(pairs(lsmeans(m.art.x2, ~ X2)))
x2.contrasts.art$d = x2.contrasts.art$estimate / sigmaHat(m.art.x2)
x2.contrasts.art

## ------------------------------------------------------------------------
x2.contrasts.ci = confint(pairs(lsmeans(m.linear, ~ X2)))
x2.contrasts.ci = within(x2.contrasts.ci, {
    d.upper.CL = upper.CL / sigmaHat(m.linear)
    d.lower.CL = lower.CL / sigmaHat(m.linear)
    d = estimate / sigmaHat(m.linear)
})
x2.contrasts.ci

## ------------------------------------------------------------------------
x2.contrasts.art.ci = confint(pairs(lsmeans(m.art.x2, ~ X2)))
x2.contrasts.art.ci = within(x2.contrasts.art.ci, {
    d.upper.CL = upper.CL / sigmaHat(m.art.x2)
    d.lower.CL = lower.CL / sigmaHat(m.art.x2)
    d = estimate / sigmaHat(m.art.x2)
})
x2.contrasts.art.ci

## ----cohens-d-comparison-------------------------------------------------
rbind(
        cbind(x2.contrasts.ci, model="linear"), 
        cbind(x2.contrasts.art.ci, model="ART")
    ) %>%
    ggplot(aes(x=model, y=d, ymin=d.lower.CL, ymax=d.upper.CL)) +
    geom_pointrange() +
    facet_grid(contrast ~ .) + 
    coord_flip()

