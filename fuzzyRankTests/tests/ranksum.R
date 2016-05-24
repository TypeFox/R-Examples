
 library(fuzzyRankTests)

 x <- c(1, 2, 3, 4, 4, 4, 5, 6, 7)
 y <- c(4, 5, 7, 7, 8, 9, 10, 11)

 fuzzy.ranksum.test(x, y)

 wilcox.test(x, y)

 fuzzy.ranksum.test(x, y, alt = "less")

 wilcox.test(x, y, alt = "less")

 fuzzy.ranksum.ci(x, y)

 wilcox.test(x, y, conf.int = TRUE)

 fuzzy.ranksum.ci(x, y, alt = "less")

 wilcox.test(x, y, conf.int = TRUE, alt = "less")

 fuzzy.ranksum.ci(x, y, alt = "great")

 wilcox.test(x, y, conf.int = TRUE, alt = "great")

