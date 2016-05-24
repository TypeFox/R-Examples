## Mental health data: mosaics for glm() and gnm() models
library(gnm)
library(vcdExtra)
data(Mental)

# display the frequency table
(Mental.tab <- xtabs(Freq ~ mental+ses, data=Mental))

# fit independence model
# Residual deviance: 47.418 on 15 degrees of freedom
indep <- glm(Freq ~ mental+ses,
                family = poisson, data = Mental)
deviance(indep)


long.labels <- list(set_varnames = c(mental="Mental Health Status", ses="Parent SES"))
mosaic(indep,residuals_type="rstandard", labeling_args = long.labels, labeling=labeling_residuals,
       main="Mental health data: Independence")
# as a sieve diagram
mosaic(indep, labeling_args = long.labels, panel=sieve, gp=shading_Friendly,
       main="Mental health data: Independence")
 
# fit linear x linear (uniform) association.  Use integer scores for rows/cols 
Cscore <- as.numeric(Mental$ses)
Rscore <- as.numeric(Mental$mental)

# column effects model (ses)
coleff <- glm(Freq ~ mental + ses + Rscore:ses,
                family = poisson, data = Mental)
mosaic(coleff,residuals_type="rstandard", 
 labeling_args = long.labels, labeling=labeling_residuals, suppress=1, gp=shading_Friendly,
 main="Mental health data: Col effects (ses)")

# row effects model (mental)
roweff <- glm(Freq ~ mental + ses + mental:Cscore,
                family = poisson, data = Mental)
mosaic(roweff,residuals_type="rstandard", 
 labeling_args = long.labels, labeling=labeling_residuals, suppress=1, gp=shading_Friendly,
 main="Mental health data: Row effects (mental)")
               
linlin <- glm(Freq ~ mental + ses + Rscore:Cscore,
                family = poisson, data = Mental)

# compare models
anova(indep, roweff, coleff, linlin)
AIC(indep, roweff, coleff, linlin)
            
mosaic(linlin,residuals_type="rstandard", 
 labeling_args = long.labels, labeling=labeling_residuals, suppress=1, gp=shading_Friendly,
 main="Mental health data: Linear x Linear")


##  Goodman Row-Column association model fits well (deviance 3.57, df 8)
Mental$mental <- C(Mental$mental, treatment)
Mental$ses <- C(Mental$ses, treatment)
RC1model <- gnm(Freq ~ mental + ses + Mult(mental, ses),
                family = poisson, data = Mental)

mosaic(RC1model,residuals_type="rstandard", 
 labeling_args = long.labels, labeling=labeling_residuals, suppress=1, gp=shading_Friendly,
 main="Mental health data: RC1 model")
