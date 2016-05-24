anglePoint <-
function(A, angl, len)
{ # internally used by smoothArc();
  c(   A[1]-sin(angl)*len,
       A[2]-cos(angl)*len);
}

