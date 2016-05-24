fix.gold <-
function (tst.set_sorc, tst.set_trgt, nrec = -1, method = c("gold","aligns"), minlen = 5, maxlen = 40, ul_s = FALSE, ul_t = TRUE, removePt = TRUE, all = FALSE, num)
{
method = match.arg (method)

p1 = prepareData (tst.set_sorc, tst.set_trgt, nrec = nrec, minlen = minlen, maxlen = maxlen, ul_s = ul_s, ul_t = ul_t, removePt = removePt, all = all, word_align = FALSE)
len = p1 $ used
p1 = unlist (p1, recursive = FALSE)
p1 = sapply(3 : length(p1), function(x) c('null', p1[[x]]))

if (method == 'gold') readline("Now, press 'Enter' to continue and edit the matrix to enter Sure/Possible alignments (Sure=1,Possible=2).")
if (method == 'aligns') readline("Now, press 'Enter' to continue and edit the matrix to enter '3' for alignments.")

mm = sapply (1 : len, function (x) {m = matrix (0, length (p1 [[x]]) + 1, length (p1 [[x + len]]) + 1);
     m [2 : nrow (m), 1] = p1 [[x]]; m [1, 2 : ncol(m)] = p1 [[x+len]]; m [1, 1] = ''; m})
fg = mm [[num]]
fix (fg)
}
