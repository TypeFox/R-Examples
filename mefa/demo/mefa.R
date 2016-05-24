wait <- function(vign=0) {
  if (vign!=0) {
    if (vign==1) {
        ANSWER <- readline("Do you want to open the vignette? (y/n) ")
        if (substr(tolower(ANSWER), 1, 1) == "y")
        vignette("mefa", package = "mefa")}
    if (vign==2) {
        ANSWER <- readline("Do you want to check out mefa website? (y/n) ")
        if (substr(tolower(ANSWER), 1, 1) == "y")
        mefaweb()}
  } else {
    ANSWER <- readline("Please press ENTER to continue ... ")
  }
}

cat("## Demo for the 'mefa' package")

wait(1)

cat("## Load the package and the example data set")

library(mefa)
data(dol.count, dol.samp, dol.taxa)
str(dol.count)
str(dol.samp)
str(dol.taxa)

wait()

cat("## Object classes")

cat("## 'stcs'")

x1 <- stcs(dol.count)
str(x1)
unique(x1$count)

wait()

x2 <- stcs(dol.count, expand = TRUE)
str(x2)
sum(x2$count)
unique(x2$count)

wait()

x3 <- stcs(dol.count, drop.zero = TRUE)
str(x3)
unique(x3$count)

wait()

cat("## 'mefa'")

m1 <- mefa(x1)
m1
m1$xtab["LT1", ]

wait()

str(m1$xtab)
str(m1$segm)

wait()

mefa(x1, nested = TRUE)

wait()

mefa(x1, segment = FALSE)

wait()

mefa(x1, dol.samp)
mefa(x1, , dol.taxa)

wait()

m2 <- mefa(x1, dol.samp, dol.taxa)
m2
str(m2$xtab)
str(m2$samp)
str(m2$taxa)

wait()

m2.sub <- mefa(x1, dol.samp[-c(1:5), ], dol.taxa[-c(1:80), ], xtab.fixed = FALSE)
m2.sub
str(m2.sub$xtab)
str(m2.sub$samp)
str(m2.sub$taxa)

wait()

mefalogo()

wait()

cat("## S3 methods")

dim(m2)
dimnames(m2)

wait()

summary(m2)

wait()

names(summary(m2))
summary(m2)$s.rich
summary(m2)$mfill

wait()

opar <- par(mfrow = c(1, 2))
plot(m2, 1, main="A")
plot(m2, 4, type="rank", trafo="log", main="B")
par(opar)

wait()

opar <- par(mfrow = c(1, 2))
boxplot(m2, 2, main = "A")
boxplot(m2, 3, main = "B")
par(opar)

wait()

molten <- melt(m2, "method")
str(molten)

wait()

m3 <- mefa(molten, dol.samp, dol.taxa)
m3

wait()

opar <- par(mfrow = c(1, 3))
image(m3, trafo = "log", sub = "all segments", main="A")
for (i in 1:2)
    image(m3, segm = i, trafo = "log",
    sub = dimnames(m3)$segm[i], main = LETTERS[i + 1])
par(opar)

wait()

ex1 <- m2[1:20, 11:15, "fresh"]
dim(ex1)
dim(ex1$samp)
dim(ex1$taxa)


wait()

ex2 <- m2[m2$samp$method == "time"]
levels(ex2$samp$method)
ex3 <- m2[m2$samp$method == "time", drop = TRUE]
levels(ex3$samp$method)

wait()

size.5 <- as.factor(is.na(m3$taxa$size) | m3$taxa$size < 5)
levels(size.5) <- c("large", "small")
m4 <- aggregate(m3, "microhab", size.5)
t(m4$xtab)
lapply(m4$segm, t)

wait()

cat("## Writing reports")

#set.seed(1234)
#m5 <- m2[ , sample(1:dim(m2)[2], 10)]
#report(m5, "report.tex", tex = TRUE, 
#    segment = TRUE, taxa.name = 1, author.name = 2,
#    drop.redundant = 1)

cat("## see mefadocs(\"SampleReport\") also")

wait()

cat("## Data analysis")

mod.amin <- glm(m2$xtab[, "amin"] ~ .,
    data = m2$samp, family = poisson)
summary(mod.amin)

wait()

library(MASS)
mod.abu <- glm.nb(summary(m2)$s.abu ~ .^2,
    data = m2$samp)
summary(mod.abu)

wait()

prop.fr <- cbind(summary(m2[ , , "fresh"])$s.abu, summary(m2)$s.abu)
mod.fr <- glm(prop.fr ~ .^2,
    data = m2$samp, family = binomial)
summary(mod.fr)

wait()

if (require(vegan)) {
m6 <- m2[summary(m2)$s.abu != 0, , ]
m6.ado <- adonis(m6$xtab ~ .^2,
    data = m6$samp, permutations = 100)
m6.ado
}

wait()

if (require(vegan)) {
m2.cca <- cca(m2$segm[["fresh"]] ~ ., data=m2$samp, 
    subset=rowSums(m2$segm[["fresh"]]) > 0)
plot(m2.cca)
}

wait()

m.list <- list()
n1 <- rep(c("time", "quadrat"), each = 2)
n2 <- rep(c("fresh", "broken"), 2)
n3 <- paste(n1, n2, sep=".")
for (i in 1:4) {
    m.list[[n3[i]]] <-
    aggregate(m2[m2$samp$method == n1[i], , n2[i]], "microhab")
}

wait()

opar <- par(mfrow=c(2, 2))
for (i in 1:4) {
    tmp <- hclust(dist(m.list[[i]]$xtab), "ward")
    plot(tmp, main = LETTERS[i], sub = names(m.list)[i], xlab = "")
}
par(opar)

cat("## End of mefa demo\n")

wait(2)

