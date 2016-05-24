dbinom(10,10,0.8)       # make 10 straight

1 - pnbinom(4,10,0.80)  # at least 5 misses <-> at least 15 shots
pbinom(9,14,0.8)        # 9 or fewer makes in 14 tries <-> at least 15 shots
pbinom(10,15,0.8)       # 10 or fewer makes in 15 tries is INCORRECT

1 - pnbinom(4,10,0.70)  # at least 5 misses = at least 15 shots
pbinom(9,14,0.7)        # 9 or fewer makes in 14 tries -> at least 15 shots
pbinom(10,15,0.7)       # 10 or fewer makes in 15 tries is INCORRECT
