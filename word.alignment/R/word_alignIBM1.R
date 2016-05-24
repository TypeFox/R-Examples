word_alignIBM1 <-
function(file_train1, file_train2, nrec = -1, iter = 4, minlen = 5, maxlen = 40, ul_s = FALSE, ul_t = TRUE, removePt = TRUE, all = FALSE, display = c("word1","number"), dtfile = NULL, f1 = 'fa',e1 = 'en', sym = FALSE, input = FALSE)
{
display = match.arg (display)

date1 = as.POSIXlt (Sys.time(), "Iran")
a = b = count0 = count = total = i = j = e = f = g = c ()

#-----------------------Translation:f to e ----------------------
aa = prepareData (file_train1, file_train2, nrec = nrec, minlen = minlen, maxlen = maxlen, ul_s = ul_s, ul_t = ul_t, all = all, removePt = removePt, word_align = TRUE)
n1 = aa[[1]]

aa = cbind(paste('null',aa[[2]][,1]),aa[[2]][,2])
len = nrow(aa)

if (is.null(dtfile))
{
date();b = apply (aa, 1, function (x) {Vt1 = strsplit (as.character (x [1]), ' ') [[1]]; Vt2 = strsplit (as.character (x[2]), ' ') [[1]];
Vt1 = Vt1 [Vt1 != '']; Vt2 = Vt2 [Vt2 != '']; cbind (Var1 = rep.int (Vt1, length (Vt2)), Var2 = rep (Vt2, each = length (Vt1)))})

cc = vapply (b,length,FUN.VALUE=0)/2

#-------------------------- main code ---------------------------       
dd1 = data.table (g = rep (1 : len, cc), f = unlist (sapply (b, function (x) x [,1])), e = unlist (sapply (b, function (x) x [,2])), t = as.numeric (rep (1 / cc, cc)));date()

rm (b, cc)
gc ()

iteration = 0
for (iiiii in 1 : iter)
{
iteration = iteration + 1
dd1 [, count0 := t / sum(t), by = paste (g, e)]
dd1 [, t := NULL]
dd1 [, count := sum (count0), by = paste (e, f)]
dd1 [, total := sum (count0), by = f]
dd1 [, t := count/total]
dd1 [, count0 := NULL]
dd1 [, count := NULL]
dd1 [, total := NULL]
}
save (dd1,iteration, file = paste(f1, e1, nrec, iter, 'RData', sep = '.'))
if (input) return (dd1)
}
# -------------------- Using saved file ------------------------
if (! is.null (dtfile))
{
load(paste(f1, e1, nrec, iter, 'RData',sep='.'))
if (input) return (dd1)
}
#--------------------- Best alignment --------------------------
word2 = strsplit(aa,' ')[1:len]
word2 = sapply (1 : len, function (x) word2 [[x]] [nzchar (word2[[x]])])

word3 = strsplit(aa,' ')[(len+1):(2*len)]
word3 = sapply (1 : len, function (x) word3 [[x]] [nzchar (word3[[x]])])

lf = vapply(word2 ,length,FUN.VALUE=0)
le = vapply(word3 ,length,FUN.VALUE=0) 

dd1 [, i := unlist (sapply (1 : len, function (x) rep (0 : (lf [x]-1), le [x])))]
dd1 [, j := unlist (sapply (1 : len, function (x) rep (1 : (le[x]), each = lf [x])))]

d1 = dd1 [, i [ which.max (t)], by = paste (g, j)] [[2]]

c1 = c (0, cumsum (le))

if (display == 'word1')  ef = sapply (1 : len, function (x) paste (word3 [[x]], word2[[x]] [d1 [ (c1 [x] + 1) : c1 [x + 1]] + 1], sep = ':', collapse='    '))

if (display == 'number') 
{
ef = sapply (1 : len, function (x) d1 [ (c1 [x] + 1) : c1 [x + 1]])
if (sym) return (ef)
}
#------------- Expected Length of both languages----------------

ex1 = mean (lf) - 1 
ex2 = mean (le)

#------------- Vocabulary size of both languages----------------

v.s1 = length (unique (unlist (word2)))
v.s2 = length (unique (unlist (word3)))

#----------------- Word Translation Probability ----------------
dd2 = unique (dd1 [, t, by = paste (e,f)])

#names(ef) = 1 : len 

date2=as.POSIXlt(Sys.time(), "Iran")
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
mylist = list (n1 = n1, n2 = len, time = date2 - date1,
iterIBM1 = iteration, expended_l_source = ex1, expended_l_target = ex2, VocabularySize_source = v.s1,
VocabularySize_target = v.s2, word_translation_prob = dd2, word_align = ef, aa = aa)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
attr(mylist, "class") <- "alignment"
return (mylist)
}


print.alignment <-
function(x, ...) 
{
print(x $ time)
cat("Number of input sentence pairs is", x[[1]], "\n")
cat("Number of used sentence pairs is", x[[2]], "\n")
cat("Number of iterations for IBM Model 1 is ", x[[4]], "\n")
cat("Word alignment for some sentence pairs are", "\n")
sapply(1:3,function(i){cat(paste(i,x $ aa[i,1],sep=': '),'\n');cat(x $ word_align[i],'\n')})
cat("            ", ".", "\n")
cat("            ", ".", "\n")
cat("            ", ".", "\n")
sapply((length(x $ word_align) - 2) : length(x $ word_align),function(i){cat(paste(i,x $ aa[i,1],sep=': '),'\n');cat(x $ word_align[i],'\n')})
}