8*5*4 / choose(17,3)  # 1 sock of each kind means no pairs
1 - (8*5*4 / choose(17,3))  # so this is prob of getting a pair
# or do it this way
( choose(8,2)*9 + choose(5,2) * 12 + choose(4,2) * 13 +
  choose(8,3) + choose(5,3) + choose(4,3) ) / choose(17,3)  
