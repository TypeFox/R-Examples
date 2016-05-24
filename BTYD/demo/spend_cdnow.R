# First, load the appropriate data. For spend models, we need
# to know how many transactions customers made and how much
# they spent on average.

data(cdnowSummary)
ave.spend <- cdnowSummary$m.x
tot.trans <- cdnowSummary$cbs[,"x"]

# Now we can estimate model parameters. spend.LL, which is used
# by spend.EstimateParameters, will give you warning if you pass
# it a customer with zero transactions. To avoid this, we can
# remove customers with zero transactions (remember to remove them
# from both vectors, as spend.LL requires the average spend and frequency
# vectors to be of equal length)

ave.spend <- ave.spend[which(tot.trans > 0)]
tot.trans <- tot.trans[which(tot.trans > 0)]

# We will let the spend function use default starting parameters.

params <- spend.EstimateParameters(ave.spend, tot.trans)
params

# Now we can make estimation about individual customers
# The following will calculate the expected transaction
# value of a customer who spent an average of $40 over 2 transactions.
spend.expected.value(params, m.x=40, x=2)

# This would also work if we wanted to compare a vector of values:
spend.expected.value(params, m.x=30:40, x=2)

# Finally, we can plot the actual and expected average
# transaction value across customers. 
spend.plot.average.transaction.value(params, ave.spend, tot.trans)

