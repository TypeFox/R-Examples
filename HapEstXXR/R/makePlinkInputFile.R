makePlinkInputFile <-
  function (famid, patid, fid, mid, sex, trait, CHR, SNP, POS, 
            geno.matrix, linkage.file, map.file, cov.file) {
    N <- length(famid)
    ns <- length(CHR)
    if (!is.matrix(geno.matrix)) 
      geno.matrix <- as.matrix(geno.matrix)
    ## Achtung: Missing als -1 kodiert (!)
    geno.matrix <- alleleRto1 (geno.matrix, miss.val = -1) 
    # Missing unten als -1 erwartet wird.
    if ((N != length(patid)) || (N != length(fid)) || (N != length(mid)) ||
          (N != length(sex)) ||  (N != length(trait))) {
      stop ("Error: length differ.")
    }
    
    if ((ns != length(SNP)) || (ns != length(POS))) {
      stop ("Dimension of CHR, SNP and POS don't match.")
    }
    if ((dim(geno.matrix)[1] != N) || (dim(geno.matrix)[2] != ns)) {
      stop ("Dimension of geno.matrix don't match to famid or CHR.")
    }
    
    geno.matrix <- replace(geno.matrix, is.na(geno.matrix), -1)
    ## make PED file
    file.create(linkage.file) 
    for (i in 1:N)      {
      geno <- (as.numeric(geno.matrix[i, ]))
      ## Achtung: Missing al -1 kodiert (!)
      geno <- factor(geno, levels = c(-1, 0, 1, 2),
                     labels = c("0\t0", "1\t1", "2\t1", "2\t2"))
      txt1 <- paste(c(famid[i], 
                      patid[i],
                      fid[i], 
                      mid[i],
                      sex[i],
                      trait[i]),
                    collapse = "\t")
      txtgeno <- paste (geno, collapse = "\t")
      txt <- paste (txt1, txtgeno, sep = "\t")
      cat (txt, "\n", 
           sep = "", 
           file = linkage.file,
           append = TRUE)
      if ((i %% 1000) == 0) cat (i, "\n")
    }
    ## make MAP file
    file.create(map.file)
    for (i in 1:ns) {
      txtmap <- paste (as.integer(CHR[i]), 
                       as.character(SNP[i]), 
                       as.integer(POS[i]))
      cat (txtmap, "\n", 
           sep = "", 
           file = map.file, 
           append = TRUE)
    }
    ## make Covariate files
    if (!all(is.na(covariates))) {
      if (!is.matrix(covariates)) covariates <- as.matrix(covariates)
      file.create(cov.file)
      for (i in 1:N) {
        txt <- paste (c(famid[i], 
                        patid[i], 
                        covariates[i, ]), 
                      collapse = "\t")
        cat (txt,
             "\n", 
             sep = "",
             file = cov.file, 
             append = TRUE)
        if ((i %% 1000) == 0) cat (i, "\n", sep="")
      }
    }
  }
