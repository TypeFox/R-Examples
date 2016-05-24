# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


combocalc = function (objects, choose, order = FALSE, repetition = TRUE){
  if (length(objects) > 1) stop ('Incorrect objects input.')
  if (!is.numeric (objects)) stop ('Incorrect objects input.')
  if (length(choose) > 1) stop ('Incorrect order input.')
  if (!is.numeric (choose)) stop ('Incorrect order input.')
  
  if (order == TRUE){
   if (repetition == TRUE) combos = objects^choose
   if (repetition == FALSE) combos = gamma (objects+1) / gamma (objects-choose+1)
  }
  if (order == FALSE){
    if (repetition == TRUE) combos = gamma (objects+choose) / (gamma(choose+1)*gamma(objects))
    if (repetition == FALSE) combos = gamma (objects+1) / (gamma(objects-choose+1)*gamma(choose+1))
  }
  combos
}
