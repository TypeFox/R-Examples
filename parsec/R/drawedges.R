drawedges <-
function(C, vertices, ...) {
  n <- dim(C)[1]
  for(a in 1:n) for(b in 1:n) if(C[a, b]) lines(vertices$x[c(a, b)], vertices$y[c(a, b)], ...)
}
