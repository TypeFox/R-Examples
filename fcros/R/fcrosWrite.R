fcrosWrite <-
function(af, file = "fcrosResults.txt", values = TRUE, thr = 1) {
   n <- length(af$ri);
   if (values) {
      id <- 1:n
      idx <- id[af$p.value <= thr];
      idnames <- af$idnames[idx];
      ri <- af$ri[idx];
      FC2 <- af$FC2[idx];
      fVal <- af$f.value[idx];
      pVal <- af$p.value[idx];
      if (nrow(summary(af)) > 8) {
         FC <- af$FC[idx];
         af.values <- matrix(c(as.character(idnames),ri, FC, FC2, fVal, pVal), ncol = 6);
         colnames(af.values) <- c("idnames", "ri", "FC", "FC2", "f.value", "p.value");
      }
      else {
           af.values <- matrix(c(as.character(idnames), ri, FC2, fVal, pVal), ncol = 5);
           colnames(af.values) <- c("idnames", "ri", "FC2", "f.value", "p.value");
      }
      write.table(af.values, file, quote = FALSE, sep = "\t", eol = "\n", col.names = TRUE, row.names = FALSE);

   } else {
      af.params <- matrix(c((af$bounds)[1], (af$bounds)[2], (af$params)[1], (af$params)[2], (af$params)[3]), ncol = 5);
      colnames(af.params) <- c("lb","ub","delta","mean","sd");

      write.table(af.params, file, quote = FALSE, sep="\t", eol = "\n", col.names = TRUE, row.names = FALSE);

   }
}
