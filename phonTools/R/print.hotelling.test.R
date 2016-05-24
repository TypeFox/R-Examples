# Copyright (c) 2015 Santiago Barreda
# All rights reserved.

print.hotelling.test <-
function (x, ...){

  if (x$samples == 1){
    cat ('\nHotelling\'s One-Sample T2-test\n')
    cat ('Alternative hypothesis: Means are not all equal to zero.\n\n')
  }
  if (x$samples == 2){
    cat ('\nHotelling\'s Two-Sample T2-test\n')
    cat ('Alternative hypothesis: Group means are not all equal.\n\n')
  }

  cat ('f = ', x$f.value, ', df1 = ', x$df1, ', df2 = ', x$df2,', p.value = ', x$p.value, '\n\n')
}
