## ------------------------------------------------------------------------
library(rUnemploymentData)

data(df_county_unemployment)
?df_county_unemployment

head(df_county_unemployment)

## ------------------------------------------------------------------------
?boxplot
boxplot(df_county_unemployment[, c(-1, -2, -3)],
        main="USA County Unemployment Data",
        xlab="Year",
        ylab="Percent Unemployment")

## ------------------------------------------------------------------------
?county_unemployment_choropleth
county_unemployment_choropleth(year=2013)

## ------------------------------------------------------------------------
?animated_county_unemployment_choropleth
# animated_county_unemployment_choropleth()

