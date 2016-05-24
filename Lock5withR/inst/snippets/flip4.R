Flips <- do(10000) * rflip(10)
tally(~heads, data = Flips)
tally(~heads, data = Flips, format = 'percent')
tally(~heads, data = Flips, format = 'proportion')


