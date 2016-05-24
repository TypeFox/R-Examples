## ----head----------------------------------------------------------------
library(dplyr)
library(ameco)
head(ameco)

## ----find_var------------------------------------------------------------
ameco %>% 
  filter(sub.chapter == "01 Population") %>% 
  .$title %>% 
  unique()

## ----example, fig.width = 7, fig.height = 4------------------------------
library(ggplot2)

ameco %>% 
  dplyr::filter(title == "Total population",
         year == 2015,
         cntry %in% c("USA", "JPN", "DEU", "FRA", "ESP", "ITA")) %>% 
  ggplot(aes(x = reorder(country, -value), y = value / 1000)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(x = NULL, y = "Population (millions)", title = "Total population")

