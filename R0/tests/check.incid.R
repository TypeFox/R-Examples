#Loading package
library(R0)

## Data is taken from the paper by Nishiura for key transmission parameters of an institutional
## outbreak during 1918 influenza pandemic in Germany
data(Germany.1918)
Germany.1918

## check.incid will extract names from the vector and coerce them as dates
check.incid(Germany.1918)

## Had Germany.1918 not have names() set, output would have been with index dates
## To force such an output, we here impose t=1:126. 
## Erasing names(Germany.1918) would have produced the same
## If so, then the epid$t vector returned will be replacement values.
check.incid(Germany.1918, t=1:126)

## You can also choose not to provide a complete date vector, but to only
## indicated the first day of the observation, and the number of days between each
## observation. In this example we will assume a time step of 7 days.
check.incid(Germany.1918, date.first.obs="1918-01-01", time.step=7)

## Finally, if no names() are available for the dataset and date.first.obs is not provided,
## setting time.step to any integer value will generate a t vector starting 
## from 1 and incrementing by the time.step parameter.
