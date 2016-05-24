mydictionary <-
function (file_train1,file_train2, nrec = -1, iter = 15, prob = 0.8, minlen=5, maxlen = 40, ul_s = FALSE, ul_t = TRUE, lang1 = 'Farsi', lang2 = 'English', removePt = TRUE, dtfile = NULL, f1 = 'fa', e1 = 'en')
{
date1 = as.POSIXlt (Sys.time(), "Iran")

e = f = c()

if(is.null (dtfile))
{
dd1 = word_alignIBM1 (file_train1, file_train2, minlen = minlen, maxlen = maxlen, iter = iter, nrec = nrec, ul_s = ul_s, ul_t = ul_t, removePt = removePt, input = TRUE)
save(dd1,iter, file = paste(f1, e1, nrec, iter, 'RData',sep='.'))
}

if( ! is.null(dtfile))
load(paste(f1, e1, nrec, iter, 'RData',sep='.'))

u1 = unique (dd1 [round (t, 1) > prob, f, e])
fe = matrix (c (u1$f, u1$e), ncol = 2)
colnames(fe) = c (lang1, lang2)
fe = fe [order (fe [,1]),]

date2 = as.POSIXlt (Sys.time(), "Iran")
##################################################################
mylist = list (time = date2 - date1, number_input = nrec, Value_prob = prob, iterIBM1 = iter, dictionary = fe)
##################################################################
return (mylist)
}
