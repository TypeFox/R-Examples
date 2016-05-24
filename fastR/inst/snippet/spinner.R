probs <- dbinom(0:50,50,0.25)
sum(probs[probs <= dbinom(8,50,0.25)])    # sum the small probs
binom.test(8,50,0.25)                     # check with binom.test()
