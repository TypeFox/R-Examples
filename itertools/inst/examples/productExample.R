library(itertools)
library(foreach)

it <- product(a=LETTERS[1:10], b=1, x=1:3)

success <- 
  foreach(a=LETTERS[1:10], .combine='c', .final=all) %:%
    foreach(b=1, .combine='c') %:%
      foreach(x=1:3, actual=it, .combine='c') %do%
        identical(list(a=a, b=b, x=x), actual)

print(success)
