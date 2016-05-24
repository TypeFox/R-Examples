## ----include=FALSE--------------------------------------------------
require(knitr)
opts_chunk$set(engine='R', tidy = FALSE, cache = TRUE)
options(width=70)

## ----FOCUS_2006_Z_data, echo=TRUE, eval=TRUE------------------------
require(mkin)
LOD = 0.5
FOCUS_2006_Z = data.frame(
  t = c(0, 0.04, 0.125, 0.29, 0.54, 1, 2, 3, 4, 7, 10, 14, 21, 
        42, 61, 96, 124),
  Z0 = c(100, 81.7, 70.4, 51.1, 41.2, 6.6, 4.6, 3.9, 4.6, 4.3, 6.8, 
         2.9, 3.5, 5.3, 4.4, 1.2, 0.7),
  Z1 = c(0, 18.3, 29.6, 46.3, 55.1, 65.7, 39.1, 36, 15.3, 5.6, 1.1, 
         1.6, 0.6, 0.5 * LOD, NA, NA, NA),
  Z2 = c(0, NA, 0.5 * LOD, 2.6, 3.8, 15.3, 37.2, 31.7, 35.6, 14.5, 
         0.8, 2.1, 1.9, 0.5 * LOD, NA, NA, NA),
  Z3 = c(0, NA, NA, NA, NA, 0.5 * LOD, 9.2, 13.1, 22.3, 28.4, 32.5, 
         25.2, 17.2, 4.8, 4.5, 2.8, 4.4))

FOCUS_2006_Z_mkin <- mkin_wide_to_long(FOCUS_2006_Z)

## ----FOCUS_2006_Z_fits_1, echo=TRUE, fig.height=4-------------------
Z.2a <- mkinmod(Z0 = list(type = "SFO", to = "Z1"),
                Z1 = list(type = "SFO"))
m.Z.2a <- mkinfit(Z.2a, FOCUS_2006_Z_mkin, quiet = TRUE)
plot(m.Z.2a)
summary(m.Z.2a, data = FALSE)

## ----FOCUS_2006_Z_fits_2, echo=TRUE, fig.height=4-------------------
Z.2a.ff <- mkinmod(Z0 = list(type = "SFO", to = "Z1"),
                   Z1 = list(type = "SFO"),
                   use_of_ff = "max")

m.Z.2a.ff <- mkinfit(Z.2a.ff, FOCUS_2006_Z_mkin, quiet = TRUE)
plot(m.Z.2a.ff)
summary(m.Z.2a.ff, data = FALSE)

## ----FOCUS_2006_Z_fits_3, echo=TRUE, fig.height=4-------------------
Z.3 <- mkinmod(Z0 = list(type = "SFO", to = "Z1", sink = FALSE),
               Z1 = list(type = "SFO"), use_of_ff = "max")
m.Z.3 <- mkinfit(Z.3, FOCUS_2006_Z_mkin, quiet = TRUE)
plot(m.Z.3)
summary(m.Z.3, data = FALSE)

## ----FOCUS_2006_Z_fits_5, echo=TRUE, fig.height=4-------------------
Z.5 <- mkinmod(Z0 = list(type = "SFO", to = "Z1", sink = FALSE),
               Z1 = list(type = "SFO", to = "Z2", sink = FALSE),
               Z2 = list(type = "SFO"))
m.Z.5 <- mkinfit(Z.5, FOCUS_2006_Z_mkin, quiet = TRUE)
plot(m.Z.5)
summary(m.Z.5, data = FALSE)

## ----FOCUS_2006_Z_fits_6, echo=TRUE, fig.height=4-------------------
Z.FOCUS <- mkinmod(Z0 = list(type = "SFO", to = "Z1", sink = FALSE),
                   Z1 = list(type = "SFO", to = "Z2", sink = FALSE),
                   Z2 = list(type = "SFO", to = "Z3"),
                   Z3 = list(type = "SFO"))
m.Z.FOCUS <- mkinfit(Z.FOCUS, FOCUS_2006_Z_mkin, 
                     quiet = TRUE)
plot(m.Z.FOCUS)
summary(m.Z.FOCUS, data = FALSE)

## ----FOCUS_2006_Z_residuals_6, echo=TRUE----------------------------
par(mfrow = c(2, 2))
mkinresplot(m.Z.FOCUS, "Z0", lpos = "bottomright")
mkinresplot(m.Z.FOCUS, "Z1", lpos = "bottomright")
mkinresplot(m.Z.FOCUS, "Z2", lpos = "bottomright")
mkinresplot(m.Z.FOCUS, "Z3", lpos = "bottomright")

## ----FOCUS_2006_Z_fits_6_ff, echo=TRUE, fig.height=4----------------
Z.FOCUS.ff <- mkinmod(Z0 = list(type = "SFO", to = "Z1", sink = FALSE),
                   Z1 = list(type = "SFO", to = "Z2", sink = FALSE),
                   Z2 = list(type = "SFO", to = "Z3"),
                   Z3 = list(type = "SFO"), 
                   use_of_ff = "max")
m.Z.FOCUS.ff <- mkinfit(Z.FOCUS.ff, FOCUS_2006_Z_mkin, quiet = TRUE)
plot(m.Z.FOCUS.ff)
summary(m.Z.FOCUS.ff, data = FALSE)

## ----FOCUS_2006_Z_fits_7, echo=TRUE, fig.height=4-------------------
Z.mkin.1 <- mkinmod(Z0 = list(type = "SFO", to = "Z1", sink = FALSE),
                    Z1 = list(type = "SFO", to = "Z2", sink = FALSE),
                    Z2 = list(type = "SFO", to = "Z3"),
                    Z3 = list(type = "SFORB"))
m.Z.mkin.1 <- mkinfit(Z.mkin.1, FOCUS_2006_Z_mkin, quiet = TRUE)
plot(m.Z.mkin.1)
summary(m.Z.mkin.1, data = FALSE)

## ----FOCUS_2006_Z_fits_8, echo=TRUE, fig.height=4-------------------
Z.mkin.2 <- mkinmod(Z0 = list(type = "SFORB", to = "Z1", sink = FALSE),
                    Z1 = list(type = "SFO"))
m.Z.mkin.2 <- mkinfit(Z.mkin.2, FOCUS_2006_Z_mkin, quiet = TRUE)
plot(m.Z.mkin.2)
summary(m.Z.mkin.2, data = FALSE)

## ----FOCUS_2006_Z_fits_9, echo=TRUE, fig.height=4-------------------
Z.mkin.3 <- mkinmod(Z0 = list(type = "SFORB", to = "Z1", sink = FALSE),
                    Z1 = list(type = "SFO", to = "Z2", sink = FALSE),
                    Z2 = list(type = "SFO"))
m.Z.mkin.3 <- mkinfit(Z.mkin.3, FOCUS_2006_Z_mkin, quiet = TRUE)
plot(m.Z.mkin.3)
summary(m.Z.mkin.3, data = FALSE)

## ----FOCUS_2006_Z_fits_10, echo=TRUE, fig.height=4------------------
Z.mkin.4 <- mkinmod(Z0 = list(type = "SFORB", to = "Z1", sink = FALSE),
                    Z1 = list(type = "SFO", to = "Z2", sink = FALSE),
                    Z2 = list(type = "SFO", to = "Z3"),
                    Z3 = list(type = "SFO"))
m.Z.mkin.4 <- mkinfit(Z.mkin.4, FOCUS_2006_Z_mkin, 
                      quiet = TRUE)
plot(m.Z.mkin.4)
summary(m.Z.mkin.4, data = FALSE)

## ----FOCUS_2006_Z_fits_11, echo=TRUE, fig.height=4------------------
Z.mkin.5 <- mkinmod(Z0 = list(type = "SFORB", to = "Z1", sink = FALSE),
                    Z1 = list(type = "SFO", to = "Z2", sink = FALSE),
                    Z2 = list(type = "SFO", to = "Z3"),
                    Z3 = list(type = "SFORB"))
m.Z.mkin.5 <- mkinfit(Z.mkin.5, FOCUS_2006_Z_mkin, quiet = TRUE)
plot(m.Z.mkin.5)
summary(m.Z.mkin.5, data = FALSE)$bpar

## ----FOCUS_2006_Z_fits_11a, echo=TRUE-------------------------------
m.Z.mkin.5a <- mkinfit(Z.mkin.5, FOCUS_2006_Z_mkin, 
                       parms.ini = c(k_Z3_bound_free = 0),
                       fixed_parms = "k_Z3_bound_free",
                       quiet = TRUE)
summary(m.Z.mkin.5a, data = FALSE)$bpar

## ----FOCUS_2006_Z_fits_11b, echo=TRUE-------------------------------
mkinparplot(m.Z.mkin.5a)

## ----FOCUS_2006_Z_fits_11b_endpoints, echo=TRUE---------------------
endpoints(m.Z.mkin.5a)

## ----FOCUS_2006_Z_residuals_11--------------------------------------
par(mfrow = c(2, 2))
mkinresplot(m.Z.mkin.5, "Z0", lpos = "bottomright")
mkinresplot(m.Z.mkin.5, "Z1", lpos = "bottomright")
mkinresplot(m.Z.mkin.5, "Z2", lpos = "bottomright")
mkinresplot(m.Z.mkin.5, "Z3", lpos = "bottomright")

