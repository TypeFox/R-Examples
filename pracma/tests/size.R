##
##  s i z e . R  tests
##

size <- pracma::size
numel <- pracma::numel
ndims <- pracma::ndims
isempty <- pracma::isempty

identical(size(1:8), c(1, 8))
identical(size(1:8, 1), 1)
identical(size(1:8, 2), 8)
identical(size(1:8, 3), 1)
identical(size(matrix(1:12, 3, 4)), c(3L, 4L))

identical(numel(array(0, c(4,4,2))), 32)
identical(numel(1:100), 100)

identical(ndims(array(NA, c(4,4,2))), 3L)
identical(ndims(list(a=1:5)), 2L)

identical(isempty(numeric(0)), TRUE)
identical(isempty(matrix(0, 1, 0)), TRUE)
identical(isempty(matrix(0, 1, 1)), FALSE)
identical(isempty(array(NA, c(2,2,2))), FALSE)
identical(isempty(array(NA, c(2,0,2))), TRUE)
