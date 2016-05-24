### using binomial dist
1- pbinom(1,3,0.6)              # win at least 2 of 3
1- pbinom(2,5,0.6)              # win at least 3 of 5
1- pbinom(3,7,0.6)              # win at least 4 of 7
### using neg binomial dist
pnbinom(1,2,0.6)                # lose <= 1 time  before 2 wins
pnbinom(2,3,0.6)                # lose <= 2 times before 3 wins
pnbinom(3,4,0.6)                # lose <= 3 times before 4 wins
