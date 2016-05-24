####### Example 1: Estimate PCSEs for a balanced dataset. #########

# Load the data.
data(agl)

# OLS Estimation for a model of growth in OECD countries

agl.lm <- lm(growth ~ lagg1 + opengdp + openex + openimp + central + leftc +
             inter + as.factor(year), data=agl)
summary(agl.lm)

# Estimate Panel-Corrected Standard-Errors and summarize the results.

agl.pcse <- pcse(agl.lm, groupN=agl$country, groupT=agl$year)
summary(agl.pcse)


####### Example 2: Estimate PCSEs for an unbalanced dataset. #########

# Load the data.
data(aglUn)

# Note: this is the orignal agl dataset with 10 observations randomly deleted.

# OLS Estimation for a model of growth in OECD countries.

aglUn.lm <- lm(growth ~ lagg1 + opengdp + openex + openimp + central + leftc +
             inter + as.factor(year), data=aglUn)
summary(aglUn.lm)

# Estimate Panel-Corrected Standard-Errors with Pairwise Selection
# and summarize the results.

aglUn.pcse1 <- pcse(aglUn.lm, groupN=aglUn$country, groupT=aglUn$year,
                   pairwise=TRUE)
summary(aglUn.pcse1)

# Estimate Panel-Corrected Standard-Errors with Casewise Selection
# and summarize the results.

aglUn.pcse2 <- pcse(aglUn.lm, groupN=aglUn$country, groupT=aglUn$year,
                   pairwise=FALSE)
summary(aglUn.pcse2)


# Extract Panel-Corrected Variance Covariance Matrix and provide to coeftest().
cv <- vcovPC(aglUn.lm, groupN=aglUn$country, groupT=aglUn$year)
# both calls below provide the same results.
coeftest(aglUn.lm, vcov.=cv)
coeftest(aglUn.lm, vcov.=vcovPC(aglUn.lm, aglUn$country, aglUn$year))



