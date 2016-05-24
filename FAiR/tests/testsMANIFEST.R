notests <- FALSE
if(notests) q(save = "no")
stopifnot(require(FAiR))

bs <- 2 # number of bootstraps

man <- make_manifest(USJudgeRatings[,-c(1:5)], how = "unbiased")
man <- make_manifest(USJudgeRatings[,-c(1:5)], how = "mle")
man <- make_manifest(USJudgeRatings[,-c(1:5)], how = "ranks")
man <- make_manifest(USJudgeRatings[,-c(1:5)], how = "lambda")
man <- make_manifest(USJudgeRatings[,-c(1:5)], how = "mcd")

man <- make_manifest(USJudgeRatings[,-c(1:5)], how = "unbiased", shrink = TRUE)
man <- make_manifest(USJudgeRatings[,-c(1:5)], how = "mle", shrink = TRUE)
man <- make_manifest(USJudgeRatings[,-c(1:5)], how = "ranks", shrink = TRUE)
man <- make_manifest(USJudgeRatings[,-c(1:5)], how = "lambda", shrink = TRUE)
man <- make_manifest(USJudgeRatings[,-c(1:5)], how = "mcd", shrink = TRUE)

man <- make_manifest(USJudgeRatings[,-c(1:5)], how = "unbiased", bootstrap = bs)
man <- make_manifest(USJudgeRatings[,-c(1:5)], how = "mle", bootstrap = bs)
man <- make_manifest(USJudgeRatings[,-c(1:5)], how = "ranks", bootstrap = bs)
man <- make_manifest(USJudgeRatings[,-c(1:5)], how = "lambda", bootstrap = bs)
man <- make_manifest(USJudgeRatings[,-c(1:5)], how = "mcd", bootstrap = bs)

man <- make_manifest(USJudgeRatings[,-c(1:5)], how = "unbiased", bootstrap = bs, 
			shrink = TRUE)
man <- make_manifest(USJudgeRatings[,-c(1:5)], how = "mle", bootstrap = bs, 
			shrink = TRUE)
man <- make_manifest(USJudgeRatings[,-c(1:5)], how = "ranks", bootstrap = bs,
			shrink = TRUE)
man <- make_manifest(USJudgeRatings[,-c(1:5)], how = "lambda", bootstrap = bs,
			shrink = TRUE)
man <- make_manifest(USJudgeRatings[,-c(1:5)], how = "mcd", bootstrap = bs,
			shrink = TRUE)

stopifnot(require(mvnmle))
data(missvals)
man <- make_manifest(missvals)
man <- make_manifest(missvals, bootstrap = bs) # many warnings

stopifnot(require(polycor))
example(hetcor)
man <- make_manifest(data)
man <- make_manifest(data, shrink = TRUE)
man <- make_manifest(data, bootstrap = bs)
man <- make_manifest(data, bootstrap = bs, shrink = TRUE)

man <- make_manifest(data, how = "ranks")
# man <- make_manifest(data, how = "ranks", shrink = TRUE)
man <- make_manifest(data, how = "ranks", bootstrap = bs)
# man <- make_manifest(data, how = "ranks", bootstrap = bs, shrink = TRUE)

## test constructors
man <- make_manifest(covmat = Harman23.cor)
man <- make_manifest(covmat = cov2cor(ability.cov$cov), n.obs = ability.cov$n.obs,
			sds = sqrt(diag(ability.cov$cov)))
hc  <- hetcor(x1, x2, y1, y2, std.err = TRUE)
man <- make_manifest(covmat = hc)
example(CovMcd)
man <- make_manifest(covmat = c1)
example(factanal)
df <- as.data.frame(cbind(v1 = v1, v2 = v2, v3 = v3, v4 = v4, v5 = v5, v6 = v6))
man <- make_manifest(~ v1 + v2 + v3 + v4 + v5 + v6, data = df, how = "unbiased")