numerator <-                 # P(K=2 & Q=2) 
    choose(4,2) *            # pick 2 Kings
    choose(4,2) *            # pick 2 Queens
    choose(52-4-4,1) /       # pick 1 other card
    choose(52,5)

denominator <-               # P(Q = 2)
    choose(4,2) *            # pick 2 Queens
    choose(52-4,3) /         # pick 3 other cards
    choose(52,5)

numerator/denominator
