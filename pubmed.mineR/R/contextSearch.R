setGeneric("contextSearch", function(object,y) standardGeneric("contextSearch"));
setMethod("contextSearch", "Abstracts", function(object, y)
{
if (length(y) > 1){cat("\\documentclass[11pt]{article}", paste("\\title {", "combined", "search}", sep = " "), "\\begin{document}",  "\\maketitle", "%\\section{}","%\\subsection{}","Results of context search",   file = paste("combined", ".tex", sep=""), sep = "\n", append = FALSE)} else {cat("\\documentclass[11pt]{article}", paste("\\title {", y[1], "search}", sep = " "), "\\begin{document}",  "\\maketitle", "%\\section{}","%\\subsection{}","Results of context search",   file = paste(y[1], ".tex", sep=""), sep = "\n", append = FALSE)    };
temp1 = unlist(lapply(object@Abstract, "SentenceToken"));
temp7 = NULL; for (j in 1:length(y))
{temp2 = regexpr(y[j], temp1, fixed=TRUE);
temp3 = which(temp2 != -1);
if(length(temp3)== 0)
{print(paste("No result found for the term",y[j],sep =" "));
write(paste("No result found for the term",y[j],sep =" "), file = "dataout.txt", append = T)}
else {temp4 = temp1[temp3];
           temp5 = temp2[temp3];temp6 = NULL;
           for (i in 1:length(temp4)){temp6 = c(temp6,   paste("\\newline", "[", i, "]",substr(temp4[i], start = 1, stop = temp5[i] - 1), "\\textbf {", y[j], "}", substr(temp4[i], start = temp5[i] + nchar(y[j]), stop = nchar(temp4[i])), sep = "  "))   };temp7 = c(temp7,temp6);              
                if (length(y) > 1){cat(temp6, file = paste("combined", ".tex", sep=""), sep = "\n", append = TRUE)} else {cat(temp6, file = paste(y[1], ".tex", sep=""), sep = "\n", append = TRUE)};
                write(temp4, file = "companion.txt", sep = "\n", append = TRUE);}}
  if (length(y) > 1) {cat("\\end{document}", file = paste("combined", ".tex", sep=""), sep = "\n", append = TRUE)} else {cat("\\end{document}", file = paste(y[1], ".tex", sep=""), sep = "\n", append = TRUE)};})
