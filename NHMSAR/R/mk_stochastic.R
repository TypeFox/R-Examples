mk_stochastic <-
function(T) {
if (is.vector(T)) {
	T = normalise(T)
  } else {
  n = length(dim(T));
  # Copy the normaliser plane for each i.
  normaliser = apply(T,1,sum);
  normaliser = matrix(normaliser,dim(T)[1],dim(T)[2]);
  # Set zeros to 1 before dividing
  # This is valid since normaliser(i) = 0 iff T(i) = 0
  normaliser = normaliser + (normaliser==0);
  T = T / normaliser;
  }
return(T)
}
