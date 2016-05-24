#
# Convert input object 'obj' to supported matrix 'mat' if necessary.
# Returns:
#   0   - if conversion was not necessary.
#   1   - if conversion was necessary and failed.
#   mat - if conversion was necessary and succeded.
#
obtain_sparsematrix <- function(obj) {
  if (is(obj, "CsparseMatrix") | is(obj, "TsparseMatrix"))
    return(0)

  if (is(obj, "Matrix"))
    return(as(obj, "CsparseMatrix"))
  
  return(1)
}
