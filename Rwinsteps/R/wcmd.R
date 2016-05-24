wcmd <- function(title = "R2Winsteps Run", data, item1, ni,
  name1, namelen = item1 - name1, codes = 0:1, csv = "y",
  hlines = "y", tfile = NULL, arglist = NULL, anchor = NULL,
  labels = NULL, extra = NULL) {

  out <- list(title = title, data = data, item1 = item1,
    ni = ni, name1 = name1, namelen = namelen, codes = codes,
    csv = csv, hlines = hlines, tfile = tfile,
    arglist = arglist, anchor = anchor, extra = extra,
    labels = labels)

  out <- as.wcmd(out)

  return(out)
}
