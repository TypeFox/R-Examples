prepareData <-
function(file1, file2, nrec = -1, minlen = 5, maxlen = 40, ul_s = FALSE, ul_t = TRUE, all = FALSE, removePt = TRUE, word_align = TRUE)
{
s_sen = t_sen = aa = t = c()

s_sen = readLines (con <- file(file1), encoding = 'UTF-8', n = nrec, warn = FALSE)
close(con)

t_sen = readLines (con <- file(file2), encoding = 'UTF-8', n = nrec, warn = FALSE)
close(con)

if (length(s_sen) == length(t_sen))
{
for (k1 in 1 : length (s_sen)) if (s_sen[k1] == '') {t_sen [k1+1] = paste (t_sen [k1], t_sen [k1+1]); t_sen [k1] = ''}

for (k2 in 1 : length (t_sen)) if (t_sen[k2] == '') {s_sen [k2+1] = paste (s_sen [k2], s_sen [k2+1]); s_sen [k2] = ''}
}

s_sen = s_sen [nzchar (s_sen)]
t_sen = t_sen [nzchar (t_sen)]

aa = cbind(s_sen,t_sen)
len1 = nrow(aa)

#------------------------- Tokenization --------------------------
    
if (ul_s) aa[,1] = culf (aa [,1], lower = all)
if (ul_t) aa[,2] = culf (aa [,2], lower = all)

rm (s_sen, t_sen)
gc ()

aa = tokenize(aa, removePunct = removePt)
len2 = length(aa) / 2

aa4 = aa [1 : len2]
aa5 = aa [ (len2+1) : (2 * len2)]

aa= cbind (aa4, aa5)
aa = aa [apply(aa, 1, function(x) length (x[1] $ aa4) * length (x [2] $ aa5) != 0)]
len2 = length(aa) / 2

word4 = cbind (word2 = aa[1:len2], word3 = aa[ (len2 + 1) : (2 * len2)])
word4 = word4 [apply (word4, 1, function(x) length(x[1] $ word2) * length (x[2] $ word3) != 0)]
 
word2 = word4[1 : (length(word4) / 2)]
word3 = word4[(length(word4) / 2 + 1) : length(word4)]  
 
aa = cbind( sapply(word2 , paste, collapse = ' '), sapply(word3, paste, collapse = ' '))

aa = aa [apply (aa, 1, function(x) prod (vapply (strsplit (x, ' '), length, FUN.VALUE=0) >= minlen)& prod (vapply (strsplit (x, ' '), length, FUN.VALUE=0) <= maxlen) == 1) ,]

if(word_align) 
{
aa = list (len1, aa)
return(aa)
}

len2 = length(aa) / 2
aa = strsplit(aa,' ')

list1 = list (initial = len1, used = len2, sorc.tok = aa [1 : len2], trgt.tok = aa[ (len2 + 1) : (2 * len2)] )

return (list1)
}
