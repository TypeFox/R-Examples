1 - pbinom(11,20,0.25)          # 11 or fewer correct fails
1 - pbinom(11,20,1/3)           # 11 or fewer correct fails
1 - pbinom(11,20, 0.5 + 0.4 * 1/3 + 0.1 * 1/4)
