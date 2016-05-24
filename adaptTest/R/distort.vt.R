`distort.vt` <-
function (f,p1,p2) function (t) pmax(0, pmin(1, f(t)+p2-f(p1)))

