"convGausExp" <-
function (k, x, mu, tau) 
{
      m <- rep(0, length(x) * length(k))
      m <- as.matrix(.C("calcCirf", m = as.double(m), as.double(k), 
                        as.double(x), as.double(tau), as.double(mu),
                        as.integer(length(k)), 
                        as.integer(length(x)), PACKAGE="TIMP")$m)
      dim(m) <- c(length(x), length(k))
      m
}