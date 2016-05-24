
 library(fuzzyRankTests)

 x <- c(-1.2, -0.7, -0.3, 0.1, 0.2, 0.3, 0.4, 0.9, 0.9, 1.0, 1.0,
     1.1, 1.5, 1.7, 1.9, 3.5, 5.1)

 fuzzy.sign.test(x, alternative = "greater")

 x <- c(-1.2, -0.7, 0.0, 0.0, 0.0, 0.3, 0.4, 0.9, 0.9, 1.0, 1.0,
     1.1, 1.5, 1.7, 1.9, 3.5, 5.1)

 print(fuzzy.sign.test(x, alternative = "greater"), digits = 3)

 print(fuzzy.sign.test(- x, alternative = "less"), digits = 3)

 print(fuzzy.sign.test(x, alternative = "two.sided"), digits = 3)

 x2 <- c(-1.2, -0.7, 0.0, 0.0, 0.0, 0.3, 0.4, 0.9, 0.9, 1.0)

 print(fuzzy.sign.test(x2, alternative = "two.sided"), digits = 3)

 print(fuzzy.sign.test(- x2, alternative = "two.sided"), digits = 3)

 x <- c(-3.5, -2.3, -1.2, -0.7, 0.0, 0.0, 0.0, 0.0, 0.4, 0.9, 0.9, 1.0, 1.0,
     1.1, 1.9, 3.5, 5.1)

 print(fuzzy.sign.test(x, alternative = "two.sided"), digits = 3)

 x <- c(-4.1, -3.5, -2.3, -1.2, -0.7, 0.0, 0.0, 0.0, 0.0, 0.4, 0.9, 0.9, 1.0,
     1.1, 1.9, 3.5)

 print(fuzzy.sign.test(x, alternative = "two.sided"), digits = 3)

 x <- x[- length(x)]

 print(fuzzy.sign.test(x, alternative = "two.sided"), digits = 3)

 x <- seq(-2, 2)

 print(fuzzy.sign.test(x, alternative = "two.sided"), digits = 3)

 x <- x[x != 0]

 print(fuzzy.sign.test(x, alternative = "two.sided"), digits = 3)

 ##### now check with alpha #####

 print(fuzzy.sign.test(x, alternative = "two.sided", alpha = 0.75),
     digits = 3)

 x <- c(-1.2, -0.7, 0.0, 0.0, 0.0, 0.3, 0.4, 0.9, 0.9, 1.0, 1.0,
     1.1, 1.5, 1.7, 1.9, 3.5, 5.1)

 print(fuzzy.sign.test(x, alternative = "greater"), digits = 3)

 print(fuzzy.sign.test(x, alternative = "greater", alpha = 0.10),
     digits = 3)

 print(fuzzy.sign.test(x, alternative = "greater", alpha = 0.05),
     digits = 3)

 print(fuzzy.sign.test(x, alternative = "greater", alpha = 0.01),
     digits = 3)

 print(fuzzy.sign.test(x, alternative = "greater", alpha = 0.001),
     digits = 3)

 print(fuzzy.sign.test(x, alternative = "greater", alpha = 0.0001),
     digits = 3)

