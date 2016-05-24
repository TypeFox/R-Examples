# Copyright (c) 2015 Santiago Barreda
# All rights reserved.

multiplot = function (n = 2, type = 'r', sizes = rep (1/n, n), show = FALSE){
  if (n < 2) stop ('Must specify at least two panels.')
    if (length(sizes) != n) stop ('Incorrect number of panel sizes specified.')
  sizes = sizes / sum(sizes)
  if (type == 'r') layout (mat = matrix(1:n, n), heights = sizes)
  else if (type == 'c') layout (mat = matrix(1:n, 1,n), widths = sizes)
  else stop ('Invalid type. Enter r for by rows and c for by columns.')
  if (show == TRUE) layout.show(n)
}

