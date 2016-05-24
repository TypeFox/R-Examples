"flist.match" <-
  function (pattern = NULL)
{
  flist = flist.data()

  if (!is.null(pattern))
  {
    incvec = grepl(pattern, flist)
    flist = flist[incvec]
  }

  return(flist)
}
