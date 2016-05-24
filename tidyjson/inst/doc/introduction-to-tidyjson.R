## ----, echo = FALSE------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
options(dplyr.print_min = 4L, dplyr.print_max = 4L)

## ----, message = FALSE---------------------------------------------------
library(tidyjson)   # this library
library(dplyr)      # for %>% and other dplyr functions

# Define a simple JSON array of people
people <- '
[
  {
    "name": "bob",
    "age": 32
  }, 
  {
    "name": "susan", 
    "age": 54
  }
]'

# Structure the data
people %>%                  # %>% is the magrittr pipeline operator 
  gather_array %>%          # gather (stack) the array by index
  spread_values(            # spread (widen) values to widen the data.frame
    name = jstring("name"), # value of "name" becomes a character column
    age = jnumber("age")    # value of "age" becomes a numeric column
  )

## ----, message = FALSE---------------------------------------------------
library(jsonlite)
jsonlite::fromJSON(people, simplifyDataFrame = TRUE)

## ------------------------------------------------------------------------
purch_json <- '
[
  {
    "name": "bob", 
    "purchases": [
      {
        "date": "2014/09/13",
        "items": [
          {"name": "shoes", "price": 187},
          {"name": "belt", "price": 35}
        ]
      }
    ]
  },
  {
    "name": "susan", 
    "purchases": [
      {
        "date": "2014/10/01",
        "items": [
          {"name": "dress", "price": 58},
          {"name": "bag", "price": 118}
        ]
      },
      {
        "date": "2015/01/03",
        "items": [
          {"name": "shoes", "price": 115}
        ]
      }
    ]
  }
]'

## ------------------------------------------------------------------------
library(jsonlite)
# Parse the JSON into a data.frame
purch_df <- jsonlite::fromJSON(purch_json, simplifyDataFrame = TRUE)
# Examine results
purch_df

## ------------------------------------------------------------------------
str(purch_df)

## ------------------------------------------------------------------------
items <- lapply(purch_df$purchases, `[[`, "items")
prices <- lapply(items, lapply, `[[`, "price")
vapply(lapply(prices, unlist), sum, integer(1))

## ------------------------------------------------------------------------
purch_df %>% group_by(name) %>% do({
  .$purchases[[1]] %>% rowwise %>% do({
    .$items[, "price", drop = FALSE]
    })
  }) %>% summarize(price = sum(price))

## ------------------------------------------------------------------------
purch_items <- purch_json %>%
  gather_array %>%                                     # stack the users 
  spread_values(person = jstring("name")) %>%          # extract the user name
  enter_object("purchases") %>% gather_array %>%       # stack the purchases
  spread_values(purchase.date = jstring("date")) %>%   # extract the purchase date
  enter_object("items") %>% gather_array %>%           # stack the items
  spread_values(                                       # extract item name and price
    item.name = jstring("name"),
    item.price = jnumber("price")
  ) %>%
  select(person, purchase.date, item.name, item.price) # select only what is needed

## ------------------------------------------------------------------------
purch_items

## ------------------------------------------------------------------------
purch_items %>% group_by(person) %>% summarize(spend = sum(item.price))

## ------------------------------------------------------------------------
# Using a single character string
x <- '{"key": "value"}' %>% as.tbl_json
x
attr(x, "JSON")

## ------------------------------------------------------------------------
# Using a vector of JSON strings
y <- c('{"key1": "value1"}', '{"key2": "value2"}') %>% as.tbl_json
y

## ------------------------------------------------------------------------
attr(y, "JSON")

## ------------------------------------------------------------------------
df <- data.frame(
  x = 1:2,
  JSON = c('{"key1": "value1"}', '{"key2": "value2"}'),
  stringsAsFactors = FALSE
) 
z <- df %>% as.tbl_json(json.column = "JSON")
z
attr(z, "JSON")

## ------------------------------------------------------------------------
c('{"a": 1}', '[1, 2]', '"a"', '1', 'true', 'null') %>% json_types

## ------------------------------------------------------------------------
'[1, "a", {"k": "v"}]' %>% gather_array %>% json_types

## ------------------------------------------------------------------------
'{"name": "bob", "age": 32}' %>% gather_keys %>% json_types

## ------------------------------------------------------------------------
'{"name": {"first": "bob", "last": "jones"}, "age": 32}' %>%
  spread_values(
    first.name = jstring("name", "first"), 
    age = jnumber("age")
  )

## ------------------------------------------------------------------------
'{"first": "bob", "last": "jones"}' %>% 
  gather_keys() %>%
  append_values_string()

## ------------------------------------------------------------------------
c('{"name": "bob", "children": ["sally", "george"]}', '{"name": "anne"}') %>% 
  spread_values(parent.name = jstring("name")) %>%
  enter_object("children") %>% 
  gather_array %>% 
  append_values_string("children")

## ------------------------------------------------------------------------
c('[1, 2, 3]', '{"k1": 1, "k2": 2}', '1', {}) %>% json_lengths

## ------------------------------------------------------------------------
'{"key": "value", "array": [1, 2, 3]}' %>% prettify

## ------------------------------------------------------------------------
'{"key": "value", "array": [1, 2, 3]}' %>% 
  gather_keys %>% json_types %>% json_lengths

## ------------------------------------------------------------------------
library(jsonlite)
worldbank[1] %>% prettify

## ------------------------------------------------------------------------
amts <- worldbank %>%
  spread_values(
    total = jnumber("totalamt")
  ) %>% 
  enter_object("majorsector_percent") %>% gather_array %>%
  spread_values(
    sector = jstring("Name"),
    pct = jnumber("Percent")
  ) %>%
  mutate(total.m = total / 10^6) %>%
  select(document.id, sector, total.m, pct) %>%
  tbl_df 
amts

## ------------------------------------------------------------------------
amts %>% 
  group_by(document.id) %>%
  summarize(pct.total = sum(pct)) %>%
  group_by(pct.total) %>%
  tally

## ------------------------------------------------------------------------
summary(amts$total.m)

## ------------------------------------------------------------------------
amts %>%
  group_by(sector) %>%
  summarize(
    spend.portion = sum(total.m * pct / 100)
  ) %>%
  ungroup %>%
  mutate(spend.dist = spend.portion / sum(spend.portion)) %>%
  arrange(desc(spend.dist))

## ----, fig.width = 7, fig.height = 6-------------------------------------
library(ggplot2)
key_stats <- companies %>% 
  gather_keys %>% json_types %>% group_by(key, type) %>% tally
key_stats
ggplot(key_stats, aes(key, n, fill = type)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip()

## ----, fig.width = 7, fig.height = 2-------------------------------------
companies %>%
  enter_object("funding_rounds") %>%
  gather_array %>% 
  gather_keys %>% json_types %>% group_by(key, type) %>% tally %>%
  ggplot(aes(key, n, fill = type)) +
    geom_bar(stat = "identity", position = "stack") +
    coord_flip()

## ------------------------------------------------------------------------
rounds <- companies %>%
  spread_values(
    id = jstring("_id", "$oid"),
    name = jstring("name"),
    category = jstring("category_code")
  ) %>%
  enter_object("funding_rounds") %>%
  gather_array %>%
  spread_values(
    round = jstring("round_code"),
    raised = jnumber("raised_amount")
  )
rounds %>% glimpse

## ----, fig.width = 7, fig.height = 2-------------------------------------
rounds %>%
  filter(
    !is.na(raised),
    round %in% c('a', 'b', 'c'),
    category %in% c('enterprise', 'software', 'web')
  ) %>%
  group_by(category, round) %>%
  summarize(raised = mean(raised)) %>%
  ggplot(aes(round, raised / 10^6, fill = round)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(y = "Raised (m)") +
    facet_grid(. ~ category)

