#
# Test warning when normalization is not possible
#

library(usl)

d <- data.frame(clients=c(   2,    4,    8,    12,    16,    20,    24,    30),
                reqrate=c(22.7, 45.4, 76.9, 109.3, 100.0, 137.6, 143.2, 145.3))

try(usl(reqrate ~ clients, d))
