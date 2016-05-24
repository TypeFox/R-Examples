
 library(fuzzyRankTests)

 x <- c(-1.2, -0.7, -0.3, 0.1, 0.2, 0.3, 0.4, 0.9, 0.9, 1.0, 1.0,
     1.1, 1.5, 1.7, 1.9, 3.5, 5.1)

 fuzzy.sign.ci(x)

 x <- c(-1.2, -0.7, 0.2, 0.2, 0.2, 0.3, 0.4, 0.9, 0.9, 1.0, 1.0,
     1.1, 1.5, 1.7, 1.9, 3.5, 5.1)

 fuzzy.sign.ci(x)

 x <- c(-1.2, -0.7, 0.2, 0.2, 0.2, 0.2, 0.4, 0.9, 0.9, 1.0, 1.0,
     1.1, 1.5, 1.7, 1.9, 3.5, 5.1)

 fuzzy.sign.ci(x)

 x <- c(-1.2, -0.7, -0.3, 0.1, 0.2, 0.3, 0.3, 0.3, 0.9, 1.0, 1.0,
     1.1, 1.5, 1.7, 1.9, 3.5, 5.1)

 fuzzy.sign.ci(x)

 x <- seq(1, 6)

 fuzzy.sign.ci(x)

 x <- seq(1, 4)

 fuzzy.sign.ci(x)

 fuzzy.sign.ci(x, conf.level = 0.5)

 fuzzy.sign.ci(x, conf.level = 0.25)

 x <- seq(1, 5)

 fuzzy.sign.ci(x)

 fuzzy.sign.ci(x, conf.level = 0.5)

 fuzzy.sign.ci(x, conf.level = 0.25)

 fuzzy.sign.ci(x, conf.level = 0.1)

 fuzzy.sign.ci(x, conf.level = 0.01)

 fuzzy.sign.ci(x, conf.level = 0.0)

 fuzzy.sign.ci(x, conf.level = 1.0)

 ##### one-sided #####

 x <- c(-1.2, -0.7, -0.3, 0.1, 0.2, 0.3, 0.4, 0.9, 0.9, 1.0, 1.0,
     1.1, 1.5, 1.7, 1.9, 3.5, 5.1)

 fuzzy.sign.ci(x, alt = "less")

 fuzzy.sign.ci(x, alt = "great")

 fuzzy.sign.ci(x, alt = "great", conf = 0.5)

 fuzzy.sign.ci(x, alt = "great", conf = 0.25)

 fuzzy.sign.ci(x, alt = "great", conf = 0)

 fuzzy.sign.ci(x, alt = "great", conf = 1)

 ##### example from paper #####

 x <- c(3.13, 3.48, 3.50, 4.70, 4.76, 4.82, 5.28, 5.67, 5.82, 8.67)

 fuzzy.sign.ci(x)

