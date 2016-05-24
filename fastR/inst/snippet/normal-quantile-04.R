x <- rbinom(40,50,0.4); x         # sample of size 40 from Binom(50,0.4)
myplot <- xqqmath(~x,fitline=TRUE)
