# Occupational status data from Goodman (1979) and Duncan (1979)
# Fit a variety of models.  Compare mosaic using expected= to mosaic.glm

library(gnm)
library(vcdExtra)
data(occupationalStatus)
occupationalStatus

# define long labels for use in mosaics
long.labels <- list(set_varnames = c(origin="origin: Father's status", 
                                     destination="destination: Son's status"))

mosaic(occupationalStatus, shade=TRUE, main="Occupational status: Independence model", 
    labeling_args = long.labels, legend=FALSE)


indep.glm <- glm(Freq ~ origin + destination, family=poisson, data=occupationalStatus)
mosaic(indep.glm, main="Independence model",
  labeling_args = long.labels, legend=FALSE, gp=shading_Friendly)


quasi.independence <-  gnm(Freq ~ origin + destination + Diag(origin,destination), family=poisson, data=occupationalStatus)
str(quasi.independence$data)
anova(quasi.independence)

## BUGLET (vcd): the diagonal cells should be considered 0 here--- outlined in black
mosaic(occupationalStatus, expected=fitted(quasi.independence), main="Quasi-independence model",
    labeling_args = long.labels, legend=FALSE, gp=shading_Friendly)


## using mosaic.gnm
 mosaic(quasi.independence, main="Quasi-independence model",
  labeling_args = long.labels, legend=FALSE, gp=shading_Friendly, labeling=labeling_residuals)


symmetry <- glm(Freq ~ Symm(origin, destination), family=poisson, data=occupationalStatus)
mosaic(occupationalStatus, expected=fitted(symmetry), main="Symmetry model",
	gp=shading_Friendly, labeling=labeling_residuals, labeling_args = long.labels )

# using mosaic.glm --- OK
mosaic(symmetry, main="Symmetry model",
	gp=shading_Friendly, labeling=labeling_residuals, labeling_args = long.labels )

quasi.symm <- glm(Freq ~ origin + destination + Symm(origin, destination), family=poisson, data=occupationalStatus)
anova(quasi.symm)
mosaic(occupationalStatus, expected=fitted(quasi.symm), main="Quasi-symmetry model")

# model comparisons
anova(independence, quasi.independence, quasi.symm, test="Chisq")

# compare symmetry to quasi summetry
anova(symmetry, quasi.symm, test="Chisq")

# association models
# uniform association, aka linear x linear association
Rscore <- as.vector(row(occupationalStatus))
Cscore <- as.vector(col(occupationalStatus))
uniform <- gnm(Freq ~ origin + destination + Rscore:Cscore, 
	family=poisson, data=occupationalStatus)
mosaic(uniform, main="Uniform association model", labeling_args = long.labels, legend=FALSE,
	gp=shading_Friendly, labeling=labeling_residuals )

RChomog <- gnm(Freq ~ origin + destination + Diag(origin, destination) + MultHomog(origin, destination),
	family=poisson, data=occupationalStatus)
mosaic(RChomog, main="RC homogeneous model", labeling_args = long.labels, legend=FALSE,
	gp=shading_Friendly, labeling=labeling_residuals)

# RC1 - heterogeneous association
RC1 <- gnm(Freq ~ origin + destination + Diag(origin, destination) + Mult(origin, destination),
	family=poisson, data=occupationalStatus)
mosaic(RC1, main="RC heterogeneous model", labeling_args = long.labels, legend=FALSE,
	gp=shading_Friendly, labeling=labeling_residuals)
