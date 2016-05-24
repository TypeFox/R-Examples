Evaluation1 <-
function(file_train1, file_train2, nrec = -1, tst.set_sorc, tst.set_trgt, nlen = -1, minlen1 = 5, maxlen1 = 40, minlen2 = 5, maxlen2 = 40,ul_s = FALSE, ul_t = TRUE, removePt = TRUE, all = FALSE, iter = 3, method = c("fix","Excel"), agn = c("an.agn","my.agn"), null.tokens = TRUE, excel1 = 'gold.xlsx', excel2 = 'align.xlsx', fixfile_gld = NULL, fixfile_agn = NULL, dtfile = NULL, f1 = 'fa', e1 = 'en', alpha = 0.5)
{
date1 = as.POSIXlt (Sys.time(), "Iran")
e = f = g = c()

method = match.arg(method)
agn = match.arg(agn)

p1 = prepareData(tst.set_sorc, tst.set_trgt, nrec = nlen, minlen = minlen2, maxlen = maxlen2, ul_s = ul_s, ul_t = ul_t, removePt = removePt, all = all, word_align = FALSE)
len = p1 $ used
# ------------------------- using fix.gold function --------------------------------
if (method == 'fix')
{
#------------------ Constructing a gold standard for the first time ----------------
if (is.null(fixfile_gld))
{
sure = possible = c()
for (sn in 1 : len)
{
mm = fix.gold (tst.set_sorc, tst.set_trgt, nrec = nlen, , method = "gold", minlen = minlen2, maxlen = maxlen2, ul_s = ul_s, ul_t = ul_t, num = sn)
save(mm,file = paste(sn,'RData',sep='.'))
sure = c (sure, t(outer (mm [1, -1], mm [-1, 1], paste))[mm [-1, -1] == 1])
possible = c (possible, t(outer (mm [1, -1], mm [-1, 1], paste)) [mm [-1, -1] == 2])
}
}
#----------------- Constructing a gold standard for later time ---------------------
if( ! is.null(fixfile_gld))
{
sure = possible = c()
for (sn in 1 : len)
{
load( paste (sn, 'RData', sep='.'))
sure = c (sure, t(outer (mm [1, -1], mm [-1, 1], paste))[mm [-1, -1] == 1])
possible = c (possible, t(outer (mm [1, -1], mm [-1, 1], paste)) [mm [-1, -1] == 2])
}
}
#-------------- Constucting a word alignment using another software ---------------
if (agn == 'an.agn')
{
#-------------- for the first time ---------------
if (is.null (fixfile_agn))
{
align = c()
for(sn in 1 : len)
{
m = fix.gold(tst.set_sorc, tst.set_trgt, method = "aligns", num = sn)
save (m, file = paste('agn',sn,'RData',sep='.'))
align = c(align, t(outer(m[1,-1],m[-1, 1],paste))[m[-1,-1] == 3])
}
}
#---------------- for later times ----------------
if (! is.null (fixfile_agn))
{
align = c()
for(sn in 1 : len)
{
load ( paste('agn',sn,'RData',sep='.'))
align = c(align, t(outer(m[1, -1],m[-1,1],paste))[m[-1,-1] == 3])
}
}
}
}
# ------------------------- using consExcel function --------------------------------
if (method == 'Excel')
{
#readline('Please ensure that you create the excel file(s) by consExcel function in this package and you fill the "excel1" file by "1|2" for sure|possible and then press "Enter" to continue.')
# ------------------------Constructing a gold standard ------------------------------
sure = possible = c()
for(sn in 1 : len)
{
df1 <- read.xlsx (xlsxFile = excel1, sheet = sn)
mm = as.matrix (df1)
sure = c (sure, t(outer (mm[1, -1], mm [-1, 1], paste)) [mm[-1, -1] == 1])
possible = c (possible, t(outer (mm [1, -1], mm [-1, 1], paste)) [mm [-1, -1] == 2])
}
#-------------- Constucting a word alignment using another softeware ---------------
if(agn == 'an.agn')
{
readline('Also please ensure that you fill the "excel2" file by "3" for alignments and then press "Enter" to continue.')
align = c()
for(sn in 1 : len)
{
df1 <- read.xlsx(xlsxFile = excel2, sheet = sn)
mm = as.matrix(df1)
align = c(align, t(outer(mm[1, -1],mm[-1, 1],paste))[mm[-1,-1] == 3])
}
}
}
#----- calculate word alignment based on my IBM model1 using word_alignIBM1 function -----
if(agn == 'my.agn')
{
#---------------------------- for the first time -----------------------------------
if (null.tokens)  p1[[3]] = sapply (1 : len, function(x) c('null', p1[[3]][[x]]))

if(is.null(dtfile))
{
dd1 = word_alignIBM1(file_train1,file_train2, iter = iter, nrec = nrec, minlen = minlen1, maxlen = maxlen1, ul_s = ul_s, ul_t = ul_t, input = TRUE, removePt = removePt, all = all)
save (dd1, iter,file = paste(f1, e1, nrec, iter,'RData', sep = '.'))
}
# ------------------------------ for later times ----------------------------------
if(! is.null(dtfile))
load(paste(f1, e1, nrec, iter,'RData', sep = '.'))
dd1[, g := NULL]
align=c()

for(i in 1:len) 
{
dd2 = unique (dd1[e %in% p1[[4]][[i]] & f %in% p1[[3]][[i]],])
dd2 = dd2[,f[which.max(t)],by = e]
align = c(align, paste(dd2[[1]],dd2[[2]]))
}
}

#---------- recall, precision and accuracy measures-----------
recall = sum(align %in% sure) / length(sure)
precision = sum(align %in% possible) / length(align)
AER = 1 - ((sum(align %in% sure) + sum(align %in% possible)) / (length(align) + length(sure)))
F_measure = 1 / (alpha/precision + (1 - alpha) / recall)
F_measure.PS = 1 / (alpha/precision + (1 - alpha) / recall)
F_measure.S = 1 / (alpha/ (sum(align %in% sure) / length(align))+ (1 - alpha) / recall)

date2 = as.POSIXlt (Sys.time(), "Iran")
#############################################################
list2 = list(time = date2 - date1, Recall = recall, Precision = precision, AER = AER, F_measure1 = F_measure.PS, F_measure2 = F_measure.S)
#############################################################
return(list2)
}
