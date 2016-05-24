cluster_words = function(wordscluster,n){words = readLines("word_table.txt");
checkword = strsplit(words, "\"", fixed = T);
checkword2 = unlist(lapply(checkword, function(x) {return(x[2])}));
final = lapply(n,  function(z)         
{y = wordscluster[[z]];
res11 =  unlist(lapply(y,function(x){temp222 = 0;  if (nchar(x) > nchar(y[1]))  { temp111 = which(checkword2 == x);  if (length(temp111) != 0)  temp222 =  as.numeric(unlist(strsplit(checkword[[temp111]][3], "\t", fixed=T))[2])};return(temp222)}));
res22 = sort(res11, decreasing = TRUE, index.return=TRUE);
res33 = res22$ix;
res44 = list(Words = y[res33], Frequencies = res22$x)
return(res44)}); attr(final,"cluster_numbers_user_given_input") = n ;  return(final)}
