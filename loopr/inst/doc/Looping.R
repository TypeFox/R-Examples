## ---- message=FALSE------------------------------------------------------
library(loopr)
library(dplyr)
library(magrittr)
library(knitr)

## ------------------------------------------------------------------------
loop = loopClass$new()

## ------------------------------------------------------------------------
id = c(1, 2, 3, 4)
toFix = c(0, 0, 1, 1)
group = c(1, 1, 1, 0)
example = data_frame(id, toFix, group)
kable(example)

## ------------------------------------------------------------------------
stack = stackClass$new()

## ------------------------------------------------------------------------
stack$push(1, name = "first")
stack$push(2, name = "second")
stack$push(3, name = "third")

## ------------------------------------------------------------------------
stack$peek

## ------------------------------------------------------------------------
stack$stack %>%
  as.data.frame %>%
  kable

## ------------------------------------------------------------------------
stack$height

## ------------------------------------------------------------------------
stack$pop
stack$pop
stack$pop

## ------------------------------------------------------------------------
stack$stack

## ------------------------------------------------------------------------
"first" %>%
loop$begin()

## ------------------------------------------------------------------------
"second" %>%
  loop$end(paste)

## ------------------------------------------------------------------------
"first" %>%
  loop$begin()

"second" %>%
  loop$cross(paste)

## ---- eval=FALSE---------------------------------------------------------
#  end(endData, FUN, ...) = FUN(stack$pop, endData, ...)
#  
#  cross(crossData, FUN, ...) = FUN(crossData, stack$pop, ...)

## ------------------------------------------------------------------------
insertData =
  example %>%
  filter(toFix == 0) %>%
  mutate(toFix = 1) %>%
  select(-group)

kable(insertData)

## ------------------------------------------------------------------------
insert(example, insertData, by = "id") %>%
  kable

## ------------------------------------------------------------------------
oldColumn1 = c(0, 0);
newColumn1 = c(1, NA)
oldColumn2 = c(0, 0);
newColumn2 = c(NA, 1)
columnData = data_frame(oldColumn1, newColumn1, oldColumn2, newColumn2)
kable(columnData)

## ------------------------------------------------------------------------
columnData %>%
  amendColumns(
    c("oldColumn1", "oldColumn2"), 
    c("newColumn1", "newColumn2")) %>%
  kable

## ------------------------------------------------------------------------
oldColumn = c(0, 0)
newColumn = c(1, NA)
columnData %>%
  fillColumns(c("newColumn1", "newColumn2"),
              c("oldColumn1", "oldColumn2")) %>%
  kable

## ------------------------------------------------------------------------
amendData = insertData

example %>%
  amend(amendData, by = "id") %>%
  kable

## ------------------------------------------------------------------------
example %>% 
  group_by(id) %>%
  amend(amendData) %>%
  kable

## ------------------------------------------------------------------------
kable(example)

## ------------------------------------------------------------------------
example %>%
  ungroup %>%
  loop$begin() %>%
    filter(group == 0) %>%
    mutate(toFix = 0) %>%
  loop$end(insert, by = "id") %>%
  kable

## ------------------------------------------------------------------------
example %>%
  group_by(group) %>%
  loop$begin() %>%
    summarize(toFix = mean(toFix)) %>%
    mutate(group = rev(group)) %>%
  loop$end(amend) %>%
  kable

## ------------------------------------------------------------------------
example %>%
  mutate(group = group + 1) %>%
  loop$begin() %>%
    names %>%
    paste0("Suffix") %>%
  loop$end(setNames) %>%
  kable

## ------------------------------------------------------------------------
example %>%
  mutate(replication = 1) %>%
  loop$begin() %>%
    mutate(replication = 2) %>%
  loop$end(bind_rows) %>%
  kable

## ------------------------------------------------------------------------
example %>%
  loop$begin(name = "original") %>%
    filter(group == 1) %>%
    loop$begin(name = "filtered") %>%
       names %>%
       paste0("Extra") %>%
    loop$end(setNames) %>%
    rename(id = idExtra) %>%
  loop$end(amend, by = "id") %>%
  kable

