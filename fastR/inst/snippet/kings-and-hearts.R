numerator <-                 # P(K=2 & H=2) 
(   1 *                      # pick King of Hearts
    choose(3,1) *            # pick one other King
    choose(12,1) *           # pick one other heart
    choose(52-1-12-3,2) +    # pick 2 other cards
#
    choose(3,2) *            # pick 2 Kings (not hearts)
    choose(12,2) *           # pick 2 hearts (not King)
    choose(52-1-12-3,1)      # pick 1 other card
) / choose(52,5)

denominator <-               # P(H = 2)
    choose(13,2) *           # pick 2 hearts
    choose(52-13,3) /        # pick 3 other cards
    choose(52,5)

numerator/denominator
