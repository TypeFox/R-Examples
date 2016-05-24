
### Regression tests for the transformation functions

set.seed(290875)
library("coin")
isequal <- coin:::isequal
options(useFancyQuotes = FALSE)


### NA handling: continuous
x <- c(1L, 2L, NA, 3L, 3L, NA, 4L, 5L, NA)
cc <- complete.cases(x)

id_trafo(x)
id_trafo(x[cc])

rank_trafo(x)
rank_trafo(x[cc])
rank_trafo(x, ties.method = "random")
rank_trafo(x[cc], ties.method = "random")

normal_trafo(x)
normal_trafo(x[cc])
normal_trafo(x, ties.method = "average-scores")
normal_trafo(x[cc], ties.method = "average-scores")

median_trafo(x)
median_trafo(x[cc])
median_trafo(x, mid.score = "0.5")
median_trafo(x[cc], mid.score = "0.5")
median_trafo(x, mid.score = "1")
median_trafo(x[cc], mid.score = "1")

savage_trafo(x)
savage_trafo(x[cc])
savage_trafo(x, ties.method = "average-scores")
savage_trafo(x[cc], ties.method = "average-scores")

consal_trafo(x)
consal_trafo(x[cc])
consal_trafo(x, a = c(2, 5))
consal_trafo(x[cc], a = c(2, 5))
consal_trafo(x, ties.method = "average-scores")
consal_trafo(x[cc], ties.method = "average-scores")
consal_trafo(x, ties.method = "average-scores", a = c(2, 5))
consal_trafo(x[cc], ties.method = "average-scores", a = c(2, 5))

koziol_trafo(x)
koziol_trafo(x[cc])
koziol_trafo(x, j = 2)
koziol_trafo(x[cc], j = 2)
koziol_trafo(x, ties.method = "average-scores")
koziol_trafo(x[cc], ties.method = "average-scores")
koziol_trafo(x, ties.method = "average-scores", j = 2)
koziol_trafo(x[cc], ties.method = "average-scores", j = 2)

klotz_trafo(x)
klotz_trafo(x[cc])
klotz_trafo(x, ties.method = "average-scores")
klotz_trafo(x[cc], ties.method = "average-scores")

mood_trafo(x)
mood_trafo(x[cc])
mood_trafo(x, ties.method = "average-scores")
mood_trafo(x[cc], ties.method = "average-scores")

ansari_trafo(x)
ansari_trafo(x[cc])
ansari_trafo(x, ties.method = "average-scores")
ansari_trafo(x[cc], ties.method = "average-scores")

fligner_trafo(x)
fligner_trafo(x[cc])
fligner_trafo(x, ties.method = "average-scores")
fligner_trafo(x[cc], ties.method = "average-scores")

maxstat_trafo(x)
maxstat_trafo(x[cc])
maxstat_trafo(x, minprob = 0.3, maxprob = 0.51)
maxstat_trafo(x[cc], minprob = 0.3, maxprob = 0.51)


### NA handling: survival
x <- c(1, 2, NA, 3, 3, NA, 4, 5, NA)
cc <- complete.cases(x)

logrank_trafo(Surv(x))
logrank_trafo(Surv(x)[cc])
logrank_trafo(Surv(x), ties.method = "Hothorn-Lausen")
logrank_trafo(Surv(x)[cc], ties.method = "Hothorn-Lausen")
logrank_trafo(Surv(x), ties.method = "average-scores")
logrank_trafo(Surv(x)[cc], ties.method = "average-scores")

x <- c(1, 2, 3, 3, 3, 4, 4, 5, 5)
e <- rep(c(0, NA, 1, 1), length.out = 9)
cc <- complete.cases(x, e)

logrank_trafo(Surv(x, e))
logrank_trafo(Surv(x, e)[cc])
logrank_trafo(Surv(x, e), ties.method = "Hothorn-Lausen")
logrank_trafo(Surv(x, e)[cc], ties.method = "Hothorn-Lausen")
logrank_trafo(Surv(x, e), ties.method = "average-scores")
logrank_trafo(Surv(x, e)[cc], ties.method = "average-scores")

x <- c(1, 2, NA, 3, 3, NA, 4, 5, NA)
e <- rep(c(0, NA, 1, 1), length.out = 9)
cc <- complete.cases(x, e)

logrank_trafo(Surv(x, e))
logrank_trafo(Surv(x, e)[cc])
logrank_trafo(Surv(x, e), ties.method = "Hothorn-Lausen")
logrank_trafo(Surv(x, e)[cc], ties.method = "Hothorn-Lausen")
logrank_trafo(Surv(x, e), ties.method = "average-scores")
logrank_trafo(Surv(x, e)[cc], ties.method = "average-scores")


### NA handling: factor
x <- factor(c(1, 1, NA, 2, NA, 3, 3, NA, 4), labels = as.roman(1:4))
ox <- ordered(x)
cc <- complete.cases(x)

f_trafo(x)
f_trafo(x[cc])

of_trafo(x)
of_trafo(x[cc])
of_trafo(x, scores = 5:8)
of_trafo(x[cc], scores = 5:8)
of_trafo(x, scores = list(s1 = 5:8, s2 = 9:12))
of_trafo(x[cc], scores = list(s1 = 5:8, s2 = 9:12))

of_trafo(ox)
of_trafo(ox[cc])
of_trafo(ox, scores = 5:8)
of_trafo(ox[cc], scores = 5:8)
of_trafo(ox, scores = list(s1 = 5:8, s2 = 9:12))
of_trafo(ox[cc], scores = list(s1 = 5:8, s2 = 9:12))

fmaxstat_trafo(x)
fmaxstat_trafo(x[cc])
fmaxstat_trafo(x, minprob = 0.49)
fmaxstat_trafo(x[cc], minprob = 0.49)

ofmaxstat_trafo(x)
ofmaxstat_trafo(x[cc])
ofmaxstat_trafo(x, minprob = 0.49)
ofmaxstat_trafo(x[cc], minprob = 0.49)

mcp_trafo(x = "Tukey")(data.frame(x))
mcp_trafo(x = "Tukey")(data.frame(x = x[cc]))

x[9] <- NA
ox[9] <- NA
cc <- complete.cases(x)

f_trafo(x)
f_trafo(x[cc])

of_trafo(x)
of_trafo(x[cc])
of_trafo(x, scores = 5:8)
of_trafo(x[cc], scores = 5:8)
of_trafo(x, scores = list(s1 = 5:8, s2 = 9:12))
of_trafo(x[cc], scores = list(s1 = 5:8, s2 = 9:12))

of_trafo(ox)
of_trafo(ox[cc])
of_trafo(ox, scores = 5:8)
of_trafo(ox[cc], scores = 5:8)
of_trafo(ox, scores = list(s1 = 5:8, s2 = 9:12))
of_trafo(ox[cc], scores = list(s1 = 5:8, s2 = 9:12))

fmaxstat_trafo(x)
fmaxstat_trafo(x[cc])
fmaxstat_trafo(x, minprob = 0.4, maxprob = 0.51)
fmaxstat_trafo(x[cc], minprob = 0.4, maxprob = 0.51)

ofmaxstat_trafo(x)
ofmaxstat_trafo(x[cc])
ofmaxstat_trafo(x, minprob = 0.4, maxprob = 0.51)
ofmaxstat_trafo(x[cc], minprob = 0.4, maxprob = 0.51)

mcp_trafo(x = "Tukey")(data.frame(x))
mcp_trafo(x = "Tukey")(data.frame(x = x[cc]))


### Weighted logrank scores
x <- c(1, 2, 3, 3, 3, 6, 6, 6, 9, 10)
e <- c(1, 0, 1, 0, 1, 1, 0, 1, 0, 1)

logrank_trafo(Surv(x, e))
logrank_trafo(Surv(x, e), ties.method = "Hothorn-Lausen")
logrank_trafo(Surv(x, e), ties.method = "average-scores")

logrank_trafo(Surv(x, e),
              type = "Gehan-Breslow")
logrank_trafo(Surv(x, e), ties.method = "Hothorn-Lausen",
              type = "Gehan-Breslow")
logrank_trafo(Surv(x, e), ties.method = "average-scores",
              type = "Gehan-Breslow")

logrank_trafo(Surv(x, e),
              type = "Tarone-Ware")
logrank_trafo(Surv(x, e), ties.method = "Hothorn-Lausen",
              type = "Tarone-Ware")
logrank_trafo(Surv(x, e), ties.method = "average-scores",
              type = "Tarone-Ware")

logrank_trafo(Surv(x, e),
              type = "Prentice")
logrank_trafo(Surv(x, e), ties.method = "Hothorn-Lausen",
              type = "Prentice")
logrank_trafo(Surv(x, e), ties.method = "average-scores",
              type = "Prentice")

logrank_trafo(Surv(x, e),
              type = "Prentice-Marek")
logrank_trafo(Surv(x, e), ties.method = "Hothorn-Lausen",
              type = "Prentice-Marek")
logrank_trafo(Surv(x, e), ties.method = "average-scores",
              type = "Prentice-Marek")

logrank_trafo(Surv(x, e),
              type = "Andersen-Borgan-Gill-Keiding")
logrank_trafo(Surv(x, e), ties.method = "Hothorn-Lausen",
              type = "Andersen-Borgan-Gill-Keiding")
logrank_trafo(Surv(x, e), ties.method = "average-scores",
              type = "Andersen-Borgan-Gill-Keiding")

logrank_trafo(Surv(x, e),
              type = "Fleming-Harrington")
logrank_trafo(Surv(x, e), ties.method = "Hothorn-Lausen",
              type = "Fleming-Harrington")
logrank_trafo(Surv(x, e), ties.method = "average-scores",
              type = "Fleming-Harrington")

logrank_trafo(Surv(x, e),
              type = "Self")
logrank_trafo(Surv(x, e), ties.method = "Hothorn-Lausen",
              type = "Self")
logrank_trafo(Surv(x, e), ties.method = "average-scores",
              type = "Self")
