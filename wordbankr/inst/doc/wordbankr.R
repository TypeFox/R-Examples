## ---- include=FALSE------------------------------------------------------
library(wordbankr)
library(knitr)
library(dplyr)
library(ggplot2)
opts_chunk$set(message = FALSE, warning = FALSE, cache = FALSE)
theme_set(theme_bw())

## ------------------------------------------------------------------------
english_ws_admins <- get_administration_data("English", "WS")
head(english_ws_admins)
all_admins <- get_administration_data()
head(all_admins)

## ------------------------------------------------------------------------
spanish_wg_items <- get_item_data("Spanish", "WG")
head(spanish_wg_items)
all_items <- get_item_data()
head(all_items)

## ------------------------------------------------------------------------
eng_ws_canines <- get_instrument_data(instrument_language = "English",
                                      instrument_form = "WS",
                                      items = c("item_26", "item_46"))
head(eng_ws_canines)

## ---- fig.width=6, fig.height=4------------------------------------------
animals <- get_item_data("English", "WS") %>%
  filter(category == "animals")

## ------------------------------------------------------------------------
animal_data <- get_instrument_data(instrument_language = "English",
                                   instrument_form = "WS",
                                   items = animals$item_id,
                                   administrations = english_ws_admins)

## ------------------------------------------------------------------------
animal_summary <- animal_data %>%
  mutate(produces = value == "produces") %>%
  group_by(age, data_id) %>%
  summarise(num_animals = sum(produces, na.rm = TRUE)) %>%
  group_by(age) %>%
  summarise(median_num_animals = median(num_animals, na.rm = TRUE))
  
ggplot(animal_summary, aes(x = age, y = median_num_animals)) +
  geom_point()

