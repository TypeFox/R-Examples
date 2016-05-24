library(testthat)
library(sValues)
context("Growth Regressions Example")

test_that("All variables - plot, coef, print",
          {
            data(economic_growth)
            eg_sv <- sValues(GR6096 ~ ., data = economic_growth)
            eg_sv <- sValues(economic_growth)
            economic_growth <- as.matrix(economic_growth)
            eg_sv <- sValues(economic_growth)
            print <- capture.output(print(eg_sv))
            print_test <- c("Data: economic_growth,    Formula: GR6096 ~ .", 
                            "R2 bounds: 0.1 - 0.5 - 1", "", "abs(S-value) > 1:", "  R2 (0.1, 1): None ", 
                            "  R2 (0.1, 0.5): None ", "  R2 (0.5, 1): BUDDHA CONFUC EAST IPRICE1 P60 RERD ", 
                            "", "abs(t-value) > 2:", "  Bayesian (R2 = 0.1): EAST ", "  Bayesian (R2 = 0.5): IPRICE1 ", 
                            "  Bayesian (R2 = 1): IPRICE1 ", "  Unconstrained OLS: IPRICE1 ")
            expect_equal(print, print_test)
            expect_warning(plot(eg_sv))
            plot(eg_sv, R2_bounds = c(0.5, 1))
            plot(eg_sv, type = "beta_plot", variable = "P60", error_bar = TRUE)
            plot(eg_sv, type = "beta_plot", variable = "P60", error_bar = TRUE, ext_bounds_shades = TRUE)
            coefs_eg <- coef(eg_sv) 
            coefs <- coefs_eg[c(1, 10, 20, 30), ]
            coefs_test <- structure(list(ols_simple = c(-0.444307543185112, 0.317848904703609, 
                                                        0.427446392923573, -0.340132299393341), b_bayes_0.1 = c(-0.0626078340226298, 
                                                                                                                -0.0103341238291921, 0.0477623111391932, -0.0326432001494661), 
                                         b_bayes_0.5 = c(-0.127764488394334, -0.0728465969961728, 
                                                         0.0501290429384938, -0.0417722974044236), b_bayes_1 = c(-0.159631480538356, 
                                                                                                                 -0.120982541843825, 0.0395733146012488, -0.0404052157558836
                                                         ), ols_all = c(-0.38276153699654, -0.329368422529548, 0.0665969050474676, 
                                                                        -3.40107574865638), t_ols_simple = c(-4.57241940143088, 3.09070163354801, 
                                                                                                             4.35916329577475, -3.33468691333559), t_bayes_0.1 = c(-1.83232400528223, 
                                                                                                                                                                   -0.28100808305331, 1.37447642448477, -0.929994575874222), 
                                         t_bayes_0.5 = c(-2.13816098407061, -0.949821797291483, 0.758257386089319, 
                                                         -0.594947428125908), t_bayes_1 = c(-2.23433504898702, -1.1811457086151, 
                                                                                            0.465735987308541, -0.426945470443487), t_ols_all = c(-2.48322648139655, 
                                                                                                                                                  -0.780337787396919, 0.192295802036046, -0.993373271117278
                                                                                            ), s_R2_0.1_1 = c(-0.444077526028489, -0.172218950280348, 
                                                                                                              0.141216331847595, -0.104206425890907), s_R2_0.1_0.5 = c(-0.550901861715267, 
                                                                                                                                                                       -0.175324866474424, 0.246912222956557, -0.17360274975807), 
                                         s_R2_0.5_1 = c(-1.96654615169665, -0.765067673951897, 0.450348250081866, 
                                                        -0.347043877176559)), .Names = c("ols_simple", "b_bayes_0.1", 
                                                                                         "b_bayes_0.5", "b_bayes_1", "ols_all", "t_ols_simple", "t_bayes_0.1", 
                                                                                         "t_bayes_0.5", "t_bayes_1", "t_ols_all", "s_R2_0.1_1", "s_R2_0.1_0.5", 
                                                                                         "s_R2_0.5_1"), row.names = c("IPRICE1", "GDPCH60L", "DENS65C", 
                                                                                                                      "GOVSH61"), class = "data.frame")
            expect_equal(coefs, coefs_test)
          }
)
          


test_that("14 favorites - plot, coef, print",
          {
            eg_sv_14 <-  sValues(GR6096 ~GDPCH60L + OTHFRAC + ABSLATIT +
                                   LT100CR + BRIT + GOVNOM1 + WARTIME +
                                   SCOUT + P60 + PRIEXP70 + OIL +
                                   H60 + POP1560 + POP6560, data = economic_growth)
            eg_sv_14
            coefs_eg_14 <- coef(eg_sv_14)
            coefs <- coefs_eg_14[,c(1,5,10, 12)]
            coefs_test <- structure(list(ols_simple = c(0.574220973453999, -0.490737489826083, 
                                                        0.404251494995227, 0.257511431206419, 0.393683750693167, 0.0758899001467707, 
                                                        0.317848904703609, -0.0726445801554096, -0.0193709519195873, 
                                                        -0.134823341581747, -0.227699055555507, 0.234415512922383, 0.307620101946417, 
                                                        0.0301227650909381), ols_all = c(0.524787655926683, -0.341624751605586, 
                                                                                         0.155041992954349, 0.145140940016124, 0.314045514456665, 0.107167468191506, 
                                                                                         -0.461969619343032, -0.0788272512782842, 0.101786060297976, -0.050200615638352, 
                                                                                         0.122043501316606, -0.101676034332394, 0.0823645761222889, -0.0318076162530259
                                                        ), t_ols_all = c(4.13969315321054, -2.60519266758056, 1.57433600400239, 
                                                                         1.49650957749142, 2.17613019611532, 1.17019634601861, -2.60798962825166, 
                                                                         -0.898798602886412, 1.14239775766554, -0.567284533392752, 0.531882448067402, 
                                                                         -0.440455571635058, 0.622380430145112, -0.35826965443253), s_R2_0.1_0.5 = c(1.61749743576061, 
                                                                                                                                                     -1.19314465523358, 0.962734936343899, 0.847562745977143, 0.797361733689061, 
                                                                                                                                                     0.699033642812357, -0.47916415858698, -0.521734105004576, 0.3019333578508, 
                                                                                                                                                     -0.287450127188205, 0.21300718533884, -0.218519949777515, 0.190790659990426, 
                                                                                                                                                     0.0059482297917709)), .Names = c("ols_simple", "ols_all", "t_ols_all", 
                                                                                                                                                                                      "s_R2_0.1_0.5"), class = "data.frame", row.names = c("P60", "PRIEXP70", 
                                                                                                                                                                                                                                           "LT100CR", "OTHFRAC", "ABSLATIT", "BRIT", "GDPCH60L", "GOVNOM1", 
                                                                                                                                                                                                                                           "OIL", "WARTIME", "POP1560", "POP6560", "H60", "SCOUT"))
            expect_equal(coefs, coefs_test)
          }
)



test_that("favorites among all",
          {
            favorites <- c("GDPCH60L", "OTHFRAC", "ABSLATIT", "LT100CR",
                           "BRIT", "GOVNOM1", "WARTIME", "SCOUT",
                           "P60", "PRIEXP70", "OIL", "H60",
                           "POP1560", "POP6560")
            eg_sv_fav <- sValues(GR6096 ~ ., data = economic_growth, R2_bounds = c(0.5, 1),
                                 favorites = favorites, R2_favorites = c(0.4, 0.8))
            print <- capture.output(print(eg_sv_fav))
            print_test <- c("Data: economic_growth,    Formula: GR6096 ~ .", "R2 bounds: 0.5 - 1", 
                            "Favorites: GDPCH60L OTHFRAC ABSLATIT LT100CR BRIT GOVNOM1 and 8 more.", 
                            "R2 favorites: 0.4 - 0.8", "", "abs(S-value) > 1:", "  R2 (0.5, 1): EAST GDPCH60L IPRICE1 OTHFRAC P60 PRIEXP70 ", 
                            "", "abs(t-value) > 2:", "  Bayesian (R2 = 0.5): P60 ", "  Bayesian (R2 = 1): P60 ", 
                            "  Unconstrained OLS: IPRICE1 ")
            expect_equal(print, print_test)
            
            plot(eg_sv_fav, R2_bounds = c(0.5, 1))
            plot(eg_sv_fav, type = "beta_plot", variable = "P60", error_bar = TRUE)
            coefs_eg_fav <- coef(eg_sv_fav)
            coefs <- coefs_eg_fav[c(1,10, 20, 30), ]
            coefs_test <- structure(list(ols_simple = c(0.574220973453999, 0.0758899001467707, 
                                                        -0.0193709519195873, -0.428949288564192), b_bayes_0.5 = c(0.309573188237554, 
                                                                                                                  0.0731104456720827, 0.0398917806422729, -0.0251870028129465), 
                                         b_bayes_1 = c(0.337842968637128, 0.0707743296948595, 0.0436644657560457, 
                                                       -0.0331150219493617), ols_all = c(0.457434290347926, 0.00902127560990199, 
                                                                                         0.0550214316388249, -0.516636953631072), t_ols_simple = c(6.46641442303491, 
                                                                                                                                                   0.701693850612851, -0.178624868589281, -4.37793861637264), 
                                         t_bayes_0.5 = c(3.35201526618303, 0.983115219116318, 0.539243828579907, 
                                                         -0.604578074647588), t_bayes_1 = c(3.16523799587686, 0.8382045304001, 
                                                                                            0.517687141216749, -0.572515740411568), t_ols_all = c(1.46799326697437, 
                                                                                                                                                  0.0429235545175516, 0.254732809421546, -1.05769127075003), 
                                         s_R2_0.5_1 = c(2.82217085352512, 0.839271061401995, 0.481037939723666, 
                                                        -0.338913067247792)), .Names = c("ols_simple", "b_bayes_0.5", 
                                                                                         "b_bayes_1", "ols_all", "t_ols_simple", "t_bayes_0.5", "t_bayes_1", 
                                                                                         "t_ols_all", "s_R2_0.5_1"), row.names = c("P60", "BRIT", "OIL", 
                                                                                                                                   "TROPICAR"), class = "data.frame")
            expect_equal(coefs, coefs_test)
          }
)

test_that("Factors", {
  data("mtcars")
  mtcars$cyl <- factor(mtcars$cyl)
  sValues(mtcars) # this should not throw an error
  
  data("economic_growth")
  e1 <- sValues(economic_growth)
  s1 <- s_values(e1)
  economic_growth$BRIT <- as.factor(economic_growth$BRIT)
  e2 <- sValues(economic_growth)
  s2 <- s_values(e2)
  row.names(s2) <- row.names(s1)
  expect_equal(s1, s2)
})