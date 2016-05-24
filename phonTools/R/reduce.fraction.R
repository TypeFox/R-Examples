# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


reduce.fraction = function (ratio){
  if (length (ratio) != 2) stop ('Ratio must contain two numbers.')
  if (ratio[1] %% 1 > 0 | ratio[2] %% 1 > 0) stop ('Elements of ratio must be whole numbers.')

  numerator = ratio[1]
  denomenator = ratio[2]

  remainder = -1
  while (remainder != 0){
    remainder = numerator %% denomenator    
    numerator  = denomenator
    if (remainder != 0) denomenator = remainder
  }
  ratio = ratio / denomenator
  return (ratio)
}

