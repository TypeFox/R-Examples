build.normalhull <-
function(inputframe) {
  ## Return a matrix containing the system [A|b], where
  ## Ax <= -b describes the convex hull
  TEMPORARY_LPFILE = "tmplp.dat"
  d = dim(inputframe)[2]
  hullinput = build.hull.list(inputframe)
  ## This is a hack, which tells quickhull to output the normals
  ## to a file.
  convhulln(hullinput, paste("n TO", TEMPORARY_LPFILE))
  input = scan(TEMPORARY_LPFILE, quiet=TRUE)
  matrix = input[3:length(input)]
  dim(matrix) = c(input[1], input[2])
  matrix = t(matrix)
  matrix
}
