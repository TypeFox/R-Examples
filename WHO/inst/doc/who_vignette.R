## ----options, echo=FALSE-------------------------------------------------
knitr::opts_chunk$set(cache = FALSE, warning = FALSE, error = FALSE, 
                      message = FALSE, fig.path = "")
library(WHO)

## ----install, eval=FALSE-------------------------------------------------
#  # From CRAN
#  install.packages("WHO")
#  
#  # From Github
#  library(devtools)
#  install_github("expersso/WHO")
#  
#  library(WHO)

## ----get_codes-----------------------------------------------------------
library(dplyr)

codes <- get_codes()
glimpse(codes)

## ----find_series---------------------------------------------------------
codes[grepl("[Ll]ife expectancy", codes$display), ]

## ----example_1, fig.width=5, fig.height=3--------------------------------
library(ggplot2)

df <- get_data("WHOSIS_000001")

head(df)

df %>% 
  filter(sex == "Both sexes") %>% 
  group_by(region, year) %>%
  summarise(value = mean(value)) %>% 
  ggplot(aes(x = year, y = value, color = region, linetype = region)) +
  geom_line(size = 1) +
  theme_light(9) +
  labs(x = NULL, y = "Life expectancy at birth (years)\n", 
       linetype = NULL, color = NULL,
       title = "Evolution of life expectancy (by region)\n")

