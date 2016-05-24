consExcel <-
function(tst.set_sorc, tst.set_trgt, method = c ("gold", "aligns"), out1 = "gold.xlsx", out2 = "align.xlsx", nrec = -1, minlen = 5, maxlen = 40, ul_s = FALSE, ul_t = TRUE, removePt = TRUE, all = FALSE)
{
method = match.arg(method)

p1 = prepareData (tst.set_sorc, tst.set_trgt, nrec = nrec, minlen = minlen, maxlen = maxlen, ul_s = ul_s, ul_t = ul_t, removePt = removePt, all = all, word_align = FALSE)
len = p1 $ used
p1 = unlist (p1, recursive = FALSE)
p1 = sapply(3 : length(p1), function(x) c('null', p1[[x]]))

if ( method == "gold")
{
wb1 <- createWorkbook ("data")
for (j in 1 : len)
{
m1 = matrix (0, length (p1 [[j]]) + 1, length (p1 [[j + len]]) + 1)
m1 [2 : nrow (m1), 1] = p1 [[j]]; m1 [1, 2 : ncol (m1)] = p1 [[j + len]]; m1 [1, 1] = ''
addWorksheet (wb1, as.character(j))
writeData (wb1, sheet =j, m1)
saveWorkbook (wb1, out1, overwrite = TRUE)
}
print ("Now, please edit 'gold.xlsx' to enter Sure/Possible alignments (Sure=1, Possible=2)")
}

if ( method == "aligns")
{
wb2 <- createWorkbook ("data")
for(j in 1 : len)
{
m1 = matrix (0, length (p1 [[j]]) + 1, length (p1 [[j + len]]) + 1)
m1 [2 : nrow (m1), 1] = p1 [[j]]; m1 [1, 2 : ncol (m1)] = p1 [[j + len]]; m1 [1,1]=''
addWorksheet (wb2, as.character(j))
writeData (wb2, sheet = j, m1)
saveWorkbook (wb2, out2, overwrite = TRUE)
}
print("Now, please edit 'align.xlsx' to enter '3' for alignments.")
}
}
