# TODO: add new weight columns to BTS demo
# load library and dataset
library(panelaggregation)
data(btsdemo)
head(btsdemo)
# adapt the levels to positive, equal and negative
# in order to suit the naming defaults. other levels work too, 
# but you'd need to specify multipliers in computeBalance then
levels(btsdemo$question_1) <- c("pos","eq","neg")

# compute the weighted shares and display store in wide format 
# to get a basis for further steps
level1 <- computeShares(btsdemo,"question_1","weight", 
                        by = c("date_qtrly","group", "altGroup", "sClass"))

# compute balance, don't have to do much here, because
# (pos, eq, neg) is the default for the possible answers
level1_wbalance <- computeBalance(level1)

# Select a particular grouping combination and a timeseries that 
# should be extracted from the level 1 aggregation.
ts1 <- extractTimeSeries(level1_wbalance,
                         "date_qtrly",
                         list(group = "C", altGroup = "a", sClass = "S"),
                         freq = 4,
                         item = "balance",
                         variable = "question_1")
ts1
# Plot a standard R ts using the plot method for ts
plot(ts1, main = attributes(ts1)$ts_key)

# Add weight column to the aggregated results
# In order to join the tables, we need to know what weight to assign to each row.
# This is done by having via a common key, for example c('group', 'altGroup').
# In this example we would assign a different weight for each 
#   c('group', 'altGroup') combination (e.g. c('A', 'a')).
btsweight1 <- btsdemo[, list(weight = sum(weight)), by = 'group']
btsagg1 <- joinDataTables(level1_wbalance, btsweight1, 'group')

# Compute second level aggregation, this time on fewer columns and using a different set of weights.
level2_balance <- computeWeightedMeans(btsagg1, c('item_pos', 'item_eq', 'item_neg', 'balance'), 
                                       'weight', c("date_qtrly","group", "sClass"))

# Select a particular grouping combination and a timeseries that 
# should be extracted from the level 2 aggregation.
ts2 <- extractTimeSeries(level2_balance,
                         "date_qtrly",
                         list(group = "C", sClass = "S"),
                         freq = 4,
                         item = "balance",
                         variable = "question_1")
ts2
# Plot a standard R ts using the plot method for ts
plot(ts2, main = attributes(ts2)$ts_key)

# Add weight column to the aggregated results
# In order to join the tables, we need to know what weight to assign to each row.
# This is done by having via a common key, for example c('group', 'altGroup').
# In this example we would assign a different weight for each 
#   c('group', 'altGroup') combination (e.g. c('A', 'a')).
btsweight2 <- btsdemo[, list(weight = sum(weight)), by = 'sClass']
btsagg2 <- joinDataTables(level2_balance, btsweight2, 'sClass')

# Compute third level of aggregation, on the whole sector, using yet another set of weights.
level3_balance <- computeWeightedMeans(btsagg2, 'balance', 'weight', c("date_qtrly", "sClass"))

# Select a particular grouping combination and a timeseries that 
# should be extracted from the level 2 aggregation.
ts3 <- extractTimeSeries(level3_balance,
                         "date_qtrly",
                         list(sClass = "S"),
                         freq = 4,
                         item = "balance",
                         variable = "question_1")
ts3
# Plot a standard R ts using the plot method for ts
plot(ts3, main = attributes(ts3)$ts_key)
