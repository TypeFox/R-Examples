## ---- message = FALSE, fig.width = 6-------------------------------------
library(survey)
library(ggplot2)
library(dplyr)

data(api)

out <- apistrat %>%
  mutate(hs_grad_pct = cut(hsg, c(0, 20, 100), include.lowest = TRUE,
                           labels = c("<20%", "20+%"))) %>%
  group_by(stype, hs_grad_pct) %>%
  summarize(api_diff = weighted.mean(api00 - api99, pw),
            n = n())

ggplot(data = out, aes(x = stype, y = api_diff, group = hs_grad_pct, fill = hs_grad_pct)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(y = 0, ymax = api_diff, label = n), position = position_dodge(width = 0.9), vjust = -1)

## ---- message = FALSE----------------------------------------------------
library(srvyr)

# simple random sample
srs_design_srvyr <- apisrs %>% as_survey_design(ids = 1, fpc = fpc)

srs_design_survey <- svydesign(ids = ~1, fpc = ~fpc, data = apisrs)

## ---- message = FALSE----------------------------------------------------
# selecting variables to keep in the survey object (stratified example)
strat_design_srvyr <- apistrat %>%
  as_survey_design(1, strata = stype, fpc = fpc, weight = pw,
                variables = c(stype, starts_with("api")))

strat_design_survey <- svydesign(~1, strata = ~stype, fpc = ~fpc,
                                 variables = ~stype + api99 + api00 + api.stu,
                                 weight = ~pw, data = apistrat)

## ---- message = FALSE----------------------------------------------------
strat_design_srvyr <- strat_design_srvyr %>%
  mutate(api_diff = api00 - api99) %>%
  rename(api_students = api.stu)

strat_design_survey$variables$api_diff <- strat_design_survey$variables$api00 -
  strat_design_survey$variables$api99
names(strat_design_survey$variables)[names(strat_design_survey$variables) == "api.stu"] <- "api_students"

## ---- message = FALSE----------------------------------------------------
# Using srvyr
out <- strat_design_srvyr %>%
  summarize(api_diff = survey_mean(api_diff, vartype = "ci"))

out

# Using survey
out <- svymean(~api_diff, strat_design_survey)

out
confint(out)

## ---- message = FALSE----------------------------------------------------
# Using srvyr
strat_design_srvyr %>%
  group_by(stype) %>%
  summarize(api_increase = survey_total(api_diff >= 0),
            api_decrease = survey_total(api_diff < 0))

# Using survey
svyby(~api_diff >= 0, ~stype, strat_design_survey, svytotal)

## ---- message = FALSE----------------------------------------------------
# Using srvyr
strat_design_srvyr %>%
  group_by(stype) %>%
  summarize(n = unweighted(n()))

# Using survey
svyby(~api99, ~stype, strat_design_survey, unwtd.count)

## ---- message = FALSE, fig.width = 6-------------------------------------
strat_design <- apistrat %>%
  as_survey_design(strata = stype, fpc = fpc, weight  = pw)

out <- strat_design %>%
  mutate(hs_grad_pct = cut(hsg, c(0, 20, 100), include.lowest = TRUE,
                           labels = c("<20%", "20+%"))) %>%
  group_by(stype, hs_grad_pct) %>%
  summarize(api_diff = survey_mean(api00 - api99, vartype = "ci"),
            n = unweighted(n()))

ggplot(data = out, aes(x = stype, y = api_diff, group = hs_grad_pct, fill = hs_grad_pct,
                       ymax = api_diff_upp, ymin = api_diff_low)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(position = position_dodge(width = 0.9), width = 0.1) +
  geom_text(aes(y = 0, label = n), position = position_dodge(width = 0.9), vjust = -1)



## ---- message = FALSE----------------------------------------------------
glm <- svyglm(api00 ~ ell + meals + mobility, design = strat_design)
summary(glm)

## ---- message = FALSE----------------------------------------------------
srs_design_srvyr <- apisrs %>% as_survey_design_(ids = "1", fpc = "fpc")

strat_design_srvyr %>%
  group_by_("stype") %>%
  summarize_(api_increase = "survey_total(api_diff >= 0)",
             api_decrease = "survey_total(api_diff < 0)")



