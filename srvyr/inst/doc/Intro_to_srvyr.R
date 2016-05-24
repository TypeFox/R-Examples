## ----message=FALSE-------------------------------------------------------
library(srvyr)
library(survey)
library(pander)
data(api)

# simple random sample
srs_design <- apisrs %>% as_survey_design(ids = 1, fpc = fpc)
srs_design

srs_design %>% 
  summarise(Total = survey_total(enroll)) %>% 
  pander()

srs_design %>% 
  summarise(Mean = survey_mean(enroll)) %>%
  pander()

# weighted sample
nofpc <- apisrs %>% as_survey_design(weights = pw)
nofpc

nofpc %>% 
  summarise(Total = survey_total(enroll)) %>% 
  pander()

nofpc %>% 
  summarise(Mean = survey_mean(enroll)) %>% 
  pander()

# stratified sample
strat_design <- apistrat %>% as_survey_design(strata = stype, fpc = fpc)
strat_design

strat_design %>% 
  summarise(Total = survey_total(enroll)) %>% 
  pander()

strat_design %>% 
  summarise(Mean=survey_mean(enroll)) %>% 
  pander()

# try with mutate
srs_design %>% 
  mutate(apidiff = api00 - api99) %>% 
  summarise(Mean = survey_mean(apidiff)) %>%
  pander()

srs_design %>% 
  mutate(apidiffpercent = (api00 - api99) / api99) %>% 
  summarise(Mean = survey_mean(apidiffpercent)) %>% 
  pander()

# try with group_by
strat_design %>% 
  group_by(stype) %>% 
  summarise(Totals = survey_total(enroll)) %>% 
  pander()


