(choose(13,2)*                # two numbers (the pairs)
choose(4,2)*choose(4,2)*      # two suits for each pair
choose(11,1)*choose(4,1)) /   # one more card of a different number
choose(52,5)
