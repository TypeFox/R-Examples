# returns -log determinant of r

`aspectDeterminant` <-
function(r) {
  list(f = -log(det(r)), g = -solve(r))
}

