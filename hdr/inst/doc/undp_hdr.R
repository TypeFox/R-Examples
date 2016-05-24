## ------------------------------------------------------------------------
library(hdr)

# Get a data frame with id and indicator names
head(hdr_indicators)

## ------------------------------------------------------------------------
hdi <- get_data(indicator = 137506, country = "DEU", year = 2013)
head(hdi)

## ------------------------------------------------------------------------
df <- get_data(103606)
head(df)

## ---- fig.width = 6, fig.height = 4--------------------------------------
br <- get_data(c(24806, 36806), year = 2010:2013)

library(dplyr)
library(tidyr)

br <- br %>% 
  group_by(id, iso3c) %>% 
  summarise(mean_val = mean(value, na.rm = TRUE)) %>% 
  spread(id, mean_val) %>% 
  .[complete.cases(.), ] %>% 
  setNames(c("iso3c", "fem_ed", "birth_rate"))

library(ggplot2)

ggplot(br, aes(x = fem_ed, y = birth_rate, label = iso3c)) +
  geom_text(size = 3, alpha = 0.75) +
  geom_smooth() +
  theme_light(8) +
  labs(y = "\nAdolescent birth rate (women aged 15-19 years)",
       x = "Population with at least secondary education, female/male ratio\n",
       title = "Relationship between female education and adolescent birth rate (2010-2013 means)")

