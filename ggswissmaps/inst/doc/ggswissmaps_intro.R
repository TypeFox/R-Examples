## ----eval=FALSE----------------------------------------------------------
#  # install.packages("devtools")
#  devtools::install_github("gibonet/ggswissmaps")

## ------------------------------------------------------------------------
library(ggswissmaps)

## ------------------------------------------------------------------------
ls("package:ggswissmaps")

## ----eval = FALSE--------------------------------------------------------
#  utils::help(package = "ggswissmaps")

## ------------------------------------------------------------------------
data(package = "ggswissmaps")

## ------------------------------------------------------------------------
data("shp_df")
class(shp_df)
length(shp_df)
names(shp_df)

# Data description
?shp_df

## ------------------------------------------------------------------------
names(maps2)

# By name
maps2[["g1k15"]]

# By index
maps2[[5]]

## ------------------------------------------------------------------------
ggplot(shp_df[["g1k15"]], aes(x = long, y = lat, group = group)) +
  geom_path() +
  coord_equal() +
  theme_white_f()

