data("war1800")

## Balance of power model
f1 <- esc + war ~ balanc + s_wt_re1 | 0 | balanc | balanc + s_wt_re1
m1 <- egame12(f1, data = war1800, subset = !is.na(regime1) & !is.na(regime2))

## Regime type model
f2 <- esc + war ~ regime1 | 0 | regime1 + regime2 | regime1 + regime2
m2 <- egame12(f2, data = war1800)

## Comparing two strategic models
vuong(model1 = m1, model2 = m2)
clarke(model1 = m1, model2 = m2)

## Comparing strategic model to logit - must specify `outcome1` appropriately
logit1 <- glm(war ~ balanc + s_wt_re1, data = m1$model, family=binomial)
vuong(model1 = m1, outcome1 = 3, model2 = logit1)
clarke(model1 = m1, outcome1 = 3, model2 = logit1)

logit2 <- glm(sq ~ regime1 + regime2, data = war1800, family=binomial)
vuong(model1 = m2, outcome1 = 1, model2 = logit2)
clarke(model1 = m2, outcome1 = 1, model2 = logit2)

## Ultimatum model
data(data_ult)
f3 <- offer + accept ~ w1 + w2 + x1 + x2 | w1 + w2 + z1 + z2
m3 <- ultimatum(f3, maxOffer = 15, data = data_ult)
ols1 <- lm(offer ~ w1 + w2 + x1 + x2 + z1 + z2, data = data_ult)
vuong(model1 = m3, outcome1 = "offer", model2 = ols1)
clarke(model1 = m3, outcome1 = "offer", model2 = ols1)
