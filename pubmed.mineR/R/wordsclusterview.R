wordsclusterview=function (words_cluster, all) 
{if (missing(all)) all = FALSE;
    write("wordclusterview", file = "words_cluster_view.txt")
    for (i in 1:length(words_cluster)) {tempccc = unlist(words_cluster[[i]]);tempddd = 1:length(tempccc);
        if (length(tempddd) >= 15) {write(c(i, length(tempccc), tempccc[1:5], tempccc[trunc(median(tempddd) - 2):trunc(median(tempddd) + 2)], tempccc[(length(tempccc) - 4):length(tempccc)]), file = "words_cluster_view.txt", ncolumns = (length(tempccc) + 2), append = T, sep = "\t"); write("   ", file = "words_cluster_view.txt", append = T)}
        else if (all) {write(c(i, length(tempccc), tempccc), file = "words_cluster_view.txt", sep = "\t", append = T, ncolumns = (length(tempccc) + 2)); write("   ", file = "words_cluster_view.txt", append = T)}}}
