## ------------------------------------------------------------------------
library(rUnemploymentData)

data(df_state_unemployment)
?df_state_unemployment

head(df_state_unemployment)

## ------------------------------------------------------------------------
?boxplot
boxplot(df_state_unemployment[, -1],
        main="USA State Unemployment Data",
        xlab="Year",
        ylab="Percent Unemployment")

## ------------------------------------------------------------------------
?state_unemployment_choropleth
state_unemployment_choropleth(year=2013)

## ------------------------------------------------------------------------
?animated_state_unemployment_choropleth
# animated_state_unemployment_choropleth()

