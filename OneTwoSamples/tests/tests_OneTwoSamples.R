library("OneTwoSamples")

## This x is used in examples 1 to 4.
set.seed(1)
x=rnorm(10, mean = 1, sd = 0.2); x

## Example 1: one normal sample, mu is known, sigma is known
one_sample(x, mu = 1, sigma = 0.2)
one_two_sample(x, mu = 1, sigma = 0.2)

interval_estimate4(x, sigma = 0.2)
mean_test1(x, mu = 1, sigma = 0.2)
interval_var3(x, mu = 1)
var_test1(x, sigma2 = 0.2^2, mu = 1)

## Example 2: one normal sample, mu is unknown, sigma is known
one_sample(x, sigma = 0.2)
one_two_sample(x, sigma = 0.2)

interval_estimate4(x, sigma = 0.2)
mean_test1(x, sigma=0.2)
interval_var3(x)
var_test1(x, sigma2 = 0.2^2)

## Example 3: one normal sample, mu is known, sigma is unknown
one_sample(x, mu = 1)
one_two_sample(x, mu = 1)

t.test(x, mu = 1)
interval_var3(x, mu = 1)
var_test1(x, mu = 1)

## Example 4: one normal sample, mu is unknown, sigma is unknown
one_sample(x)
one_two_sample(x)

t.test(x)
interval_var3(x)
var_test1(x)

## This x is used in examples 5 to 6.
set.seed(1)
x=rexp(10, rate = 1); x

## Example 5: one non-normal sample, sigma is known
one_sample(x, sigma = 1)
one_two_sample(x, sigma = 1)
interval_estimate3(x, sigma = 1)

## Example 6: one non-normal sample, sigma is unknown
one_sample(x)
one_two_sample(x)
interval_estimate3(x)

## The x, y, y2 defined below are used in examples 7-12.
set.seed(1)
x=rnorm(10, mean = 1, sd = 0.2); x
y=rnorm(20, mean = 2, sd = 0.3); y
y2=rnorm(20, mean = 2, sd = 0.2); y2

## The common characteristics in examples 7-12: two normal samples x and y, n1!=n2
## Example 7: sigma1, sigma2 are known; mu1, mu2 are known
one_two_sample(x, y, sigma = c(0.2, 0.3), mu = c(1, 2))

interval_estimate5(x, y, sigma = c(0.2, 0.3))
mean_test2(x, y, sigma = c(0.2, 0.3))
interval_var4(x, y, mu = c(1, 2))
var_test2(x, y, mu = c(1, 2))
ks.test(x, y)
wilcox.test(x, y)

## Example 8: sigma1 = sigma2 are unknown; mu1, mu2 are known
one_two_sample(x, y2, var.equal = TRUE, mu = c(1, 2))

t.test(x, y2, var.equal = TRUE)
interval_var4(x, y2, mu = c(1, 2))
var_test2(x, y2, mu = c(1, 2))
ks.test(x, y)
wilcox.test(x, y)

## Example 9: sigma1 != sigma2 are unknown; mu1, mu2 are known
one_two_sample(x, y, mu = c(1, 2))

t.test(x, y)
interval_var4(x, y, mu = c(1, 2))
var_test2(x, y, mu = c(1, 2))
ks.test(x, y)
wilcox.test(x, y)

## Example 10: sigma1, sigma2 are known; mu1, mu2 are unknown
one_two_sample(x, y, sigma = c(0.2, 0.3))

interval_estimate5(x, y, sigma = c(0.2, 0.3))
mean_test2(x, y, sigma = c(0.2, 0.3))
var.test(x, y)
ks.test(x, y)
wilcox.test(x, y)

## Example 11: sigma1 = sigma2 are unknown; mu1, mu2 are unknown
one_two_sample(x, y2, var.equal = TRUE)

t.test(x, y2, var.equal = TRUE)
var.test(x, y2)
ks.test(x, y)
wilcox.test(x, y)

## Example 12: sigma1 != sigma2 are unknown; mu1, mu2 are unknown
one_two_sample(x, y)

t.test(x, y)
var.test(x, y)
ks.test(x, y)
wilcox.test(x, y)

## The "women" data in package datasets is used in the paper.
## women$height and women$weight: they are normals!
women
x = women$height; x
y = women$weight; y
one_two_sample(x, y)

t.test(x, y)
var.test(x, y)
ks.test(x, y)
binom.test(sum(x<y), length(x))
wilcox.test(x, y, paired = TRUE)
cor.test(x, y, method = "pearson")
cor.test(x, y, method = "kendall")
cor.test(x, y, method = "spearman")

