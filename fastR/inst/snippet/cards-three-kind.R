(choose(13,1) *               # a number to have three of
choose(4,3) *                 # three of that number
choose(12,2) *                # two other numbers
choose(4,1)*choose(4,1)) /    # one of each of those numbers
choose(52,5)
