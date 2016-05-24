library(ThreeArmedTrials)

# Test for non-inferiority test. lambda_P=8, lambda_R = 4, lambda_E = 5, and phi = 1
# Delta = (lambda_P-lambda_E)/(lambda_P-lambda_R)
xExp <-rnbinom(60, mu=5, size=1)
xRef <-rnbinom(40, mu=4, size=1)
xPla <-rnbinom(40, mu=8, size=1)
Delta <- (8-5)/(8-4)
taNegbin.test(xExp, xRef, xPla, Delta, method = 'RML')
taNegbin.test(xExp, xRef, xPla, Delta, method = 'ML')
taNegbin.test(xExp, xRef, xPla, Delta, method = 'SampleVariance')

# Test for superiority test. lambda_P=8, lambda_R = 5, lambda_E = 4, and phi = 1
# Delta = (lambda_P-lambda_E)/(lambda_P-lambda_R)
xExp <-rnbinom(60, mu=5, size=1)
xRef <-rnbinom(40, mu=4, size=1)
xPla <-rnbinom(40, mu=8, size=1)
Delta <- (8-5)/(8-4)
taNegbin.test(xExp, xRef, xPla, Delta, method = 'RML')
taNegbin.test(xExp, xRef, xPla, Delta, method = 'ML')
taNegbin.test(xExp, xRef, xPla, Delta, method = 'SampleVariance')


# Test for type = 'unrestricted' & equal sample size allocation: calculation of n, power, and sig.level. Expect 1038, 0.8, 0.025
power.taNegbin.test(rateExp = 2, rateRef = 2, ratePla = 4, shape = 0.5, Delta = 0.8, sig.level = 0.025, power = 0.8, type = 'unrestricted', allocation = c(1/3, 1/3, 1/3))$n
power.taNegbin.test(rateExp = 2, rateRef = 2, ratePla = 4, shape = 0.5, Delta = 0.8, sig.level = 0.025, n = 1038, type = 'unrestricted', allocation = c(1/3, 1/3, 1/3))$power
power.taNegbin.test(rateExp = 2, rateRef = 2, ratePla = 4, shape = 0.5, Delta = 0.8, power = 0.8007362, n = 1038, type = 'unrestricted', allocation = c(1/3, 1/3, 1/3))$sig.level

# Test for type = 'restricted' & equal sample size allocation: calculation of n, power, and sig.level. Expect 1092, 0.8, 0.025
power.taNegbin.test(rateExp = 2, rateRef = 2, ratePla = 4, shape = 0.5, Delta = 0.8, sig.level = 0.025, power = 0.8, type = c('restricted'), allocation = c(1/3, 1/3, 1/3))$n
power.taNegbin.test(rateExp = 2, rateRef = 2, ratePla = 4, shape = 0.5, Delta = 0.8, sig.level = 0.025, n = 1092, type = c('restricted'), allocation = c(1/3, 1/3, 1/3))$power
power.taNegbin.test(rateExp = 2, rateRef = 2, ratePla = 4, shape = 0.5, Delta = 0.8, n = 1092, power = 0.8000113, type = c('restricted'), allocation = c(1/3, 1/3, 1/3))$sig.level

#Test for type = 'unrestricted' & sample size allocation 2:2:1 : calculation of n, power, and sig.level. Expect 923, 0.8, 0.025
power.taNegbin.test(rateExp = 2, rateRef = 2, ratePla = 4, shape = 0.5, Delta = 0.8, sig.level = 0.025, power = 0.8, type = c('unrestricted'), allocation = c(2/5, 2/5, 1/5))$n
power.taNegbin.test(rateExp = 2, rateRef = 2, ratePla = 4, shape = 0.5, Delta = 0.8, sig.level = 0.025, n = 923, type = c('unrestricted'), allocation = c(2/5, 2/5, 1/5))$power
power.taNegbin.test(rateExp = 2, rateRef = 2, ratePla = 4, shape = 0.5, Delta = 0.8, power = 0.8011693, n = 923, type = c('unrestricted'), allocation = c(2/5, 2/5, 1/5))$sig.level

# Test for type = 'restricted' & sample size allocation 2:2:1 : calculation of n, power, and sig.level. Expect 1048, 0.8, 0.025
power.taNegbin.test(rateExp = 2, rateRef = 2, ratePla = 4, shape = 0.5, Delta = 0.8, sig.level = 0.025, power = 0.8, type = c('restricted'), allocation = c(2/5, 2/5, 1/5))$n
power.taNegbin.test(rateExp = 2, rateRef = 2, ratePla = 4, shape = 0.5, Delta = 0.8, sig.level = 0.025, n = 1048, type = c('restricted'), allocation = c(2/5, 2/5, 1/5))$power
power.taNegbin.test(rateExp = 2, rateRef = 2, ratePla = 4, shape = 0.5, Delta = 0.8, n = 1048, power = 0.8012539, type = c('restricted'), allocation = c(2/5, 2/5, 1/5))$sig.level
power.taNegbin.test(rateExp = 2, rateRef = 2, ratePla = 4, shape = 0.5, Delta = 0.8, n = 1048, power = 0.8012539, type = c('restricted'), allocation = c(2/5, 2/5, 1/5))

# 
# taNegBin.OptAllocation(rateExp = 2, rateRef = 2, ratePla = 4, shape = 0.5, Delta = 0.8, type = 'unrestricted', n = NULL)
# taNegBin.OptAllocation(rateExp = 2, rateRef = 2, ratePla = 4, shape = 0.5, Delta = 0.8, type = 'unrestricted', n = 300)
# taNegBin.OptAllocation(rateExp = 2, rateRef = 2, ratePla = 4, shape = 0.5, Delta = 0.8, type = 'restricted', n = 500, sig.level = 0.025)
# 
# 
taNegBin.OptAllocation(rateExp = 2, rateRef = 2, ratePla = 4, shape = 0.5, Delta = 0.8, type = 'unrestricted', n = 1048, sig.level = 0.025)
taNegBin.OptAllocation(rateExp = 2, rateRef = 2, ratePla = 4, shape = 0.5, Delta = 0.8, type = 'unrestricted')
power.taNegbin.test(rateExp = 2, rateRef = 2, ratePla = 4, shape = 0.5, Delta = 0.8, n = 1048, power = 0.8546295, type = c('unrestricted'), allocation = c(488, 391, 169)/1048)$sig.level
# 
# 
#taNegBin.OptAllocation(rateExp = 2, rateRef = 2, ratePla = 4, shape = 0.5, Delta = 0.8, type = 'restricted', n = 500, sig.level = 0.025)
# power.taNegbin.test(rateExp = 2, rateRef = 2, ratePla = 4, shape = 0.5, Delta = 0.8, n = 500, sig.level = 0.025, type = 'restricted', allocation = c(0.530, 0.324, 0.146))

# # Expect recalculation of 'allocation' and 'n'
# power.taNegbin.test(rateExp = 2, rateRef = 2, ratePla = 4, shape = 0.5, Delta = 0.8, n = 1001, power = 0.8, allocation = c(0.25, 0.5, 0.25))
# 
# # Expect error since n, power, and sig.level are NULL
# power.taNegbin.test(rateExp = 2, rateRef = 2, ratePla = 4, shape = 0.5, Delta = 0.8, type = c('unrestricted'), allocation = c(1/3, 1/3, 1/3))
# 
# # Expect error since power is larger than 1
# power.taNegbin.test(rateExp = 2, rateRef = 2, ratePla = 4, shape = 0.5, Delta = 0.8, sig.level = 0.025, power = 1.8, type = c('unrestricted'), allocation = c(1/3, 1/3, 1/3))
# 
# # Expect error since Delta is not larger than 0
# power.taNegbin.test(rateExp = 2, rateRef = 2, ratePla = 4, shape = 0.5, Delta = 0, sig.level = 0.025, power = 0.8, type = c('unrestricted'), allocation = c(1/3, 1/3, 1/3))
# 
# # Expect error since shape is not larger than 0
# power.taNegbin.test(rateExp = 2, rateRef = 2, ratePla = 4, shape = -0.5, Delta = 0.8, sig.level = 0.025, power = 0.8, type = c('unrestricted'), allocation = c(1/3, 1/3, 1/3))
# 
# # Expect error since rateExp is missing
# power.taNegbin.test(rateRef = 2, ratePla = 4, shape = 0.5, Delta = 0.8, sig.level = 0.025, power = 0.8, type = c('unrestricted'), allocation = c(1/3, 1/3, 1/3))
