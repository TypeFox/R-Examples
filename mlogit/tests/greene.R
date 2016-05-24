library("mlogit")
data("TravelMode", package = "AER")
TravelMode$choice2 <- TravelMode$choice == "yes"
TravelMode$incair <- with(TravelMode, income * (mode == "air"))
tm_cl <- mlogit(choice2 ~ gcost + wait + incair, data = TravelMode, shape = "long", choice = "choice", alt.var = "mode", reflevel = "car")
# Greene table 21.15 first column
tm_hl <- mlogit(choice2 ~ gcost + wait + incair, data = TravelMode, shape = "long", choice = "choice", alt.var = "mode", reflevel = "car", heterosc = TRUE)
# Greene table 21.15 second column
tm_nl <- mlogit(choice2 ~ gcost + wait + incair, data = TravelMode, shape = "long", choice = "choice", alt.var = "mode", reflevel = "car", nests = list(fly = "air", ground = c("bus", "car", "train")), unscaled = TRUE)

## tm_ml <- mlogit(choice2 ~ gcost + wait + incair, data = TravelMode,
## shape = "long", choice = "choice", alt.var = "mode", reflevel =
## "car", rpar = c(altair = "n", alttrain = "n", altbus = "n", gcost =
## "n", wait = "n", incair = "n"), R = 400)

