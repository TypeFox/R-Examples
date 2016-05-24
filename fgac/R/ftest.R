"ftest" <-
function(x,y, test = c("wilcox.test", "t.test")) {
  test <- match.arg(test)
  switch(test,
         wilcox.test = wilcox.test(x,y,exact=FALSE),
         t.test = t.test(x,y,exact=FALSE))
}

