#----------------------------------
# multipoint functions
#----------------------------------

solar_multipoint <- function(dir, out, out.dir, out.chr)
{
  ### check arguments
  stopifnot(file.exists(dir))

  ### var  
  tmp.dir <- !missing(out.dir)

  trait.dir <- paste(out$traits, collapse = ".")
  model.dir <- file.path(trait.dir, out$solar$model.filename)

  ### `dir.poly`  
  dir.poly <- out$multipoint$dir.poly

  out.path <- file.path(dir, out.dir)
  stopifnot(dir.create(out.path))
  out.path <- normalizePath(out.path)

  ### create a temp. dir.
  dir.multipoint <- file.path(dir, paste0("solar_multipoint_", out.dir))
  stopifnot(dir.create(dir.multipoint))
  
  files.dir <- list.files(dir.poly, include.dirs = TRUE, full.names = TRUE)
  stopifnot(file.copy(from = files.dir, to = dir.multipoint, recursive = TRUE))
    
  ### make `cmd`
  cmd <- c(paste("load model", model.dir),
    paste("mibddir", out$multipoint$mibddir), 
    paste("chromosome", out.chr),
    paste("interval", out$multipoint$interval),
    out$multipoint$multipoint.settings,
    paste("multipoint -overwrite", out$multipoint$multipoint.options))

  ret <- solar(cmd, dir.multipoint, result = FALSE) 

  ### copy result files
  files.dir <- list.files(file.path(dir.multipoint, trait.dir), include.dirs = TRUE, full.names = TRUE)
  stopifnot(file.copy(from = files.dir, to = out.path, recursive = TRUE))
  
  ### `tab.dir`
  tab.dir <- out.path
  tab.file <- file.path(tab.dir, "multipoint.out")
  solar.ok <- file.exists(tab.file)

  ### clean
  unlink(dir.multipoint, recursive = TRUE)
    
  ### return  
  out <-  list(solar = list(cmd = cmd, solar.ok = solar.ok), tab.dir = tab.dir)

  return(out)
}

#----------------------------------
# mibd functions
#----------------------------------

# @export
solarMIBD <- function(mibddir, verbose = 0, chr, nmibd, cores = 1)
{
  ### inc
  stopifnot(requireNamespace("Matrix", quietly = TRUE))
  
  ### arg
  stopifnot(file.exists(mibddir))
  
  ### parallel
  parallel <- (cores > 1)

  if(parallel) {
    # load required R package doParallel
    stopifnot(requireNamespace("doParallel", quietly = TRUE))
    
    doParallel::registerDoParallel(cores = cores)
  }

  
  ### var
  files <- list.files(mibddir, full.names = TRUE)
  
  ### pass 1: chr, cM
  out <- llply(files, function(f) {
    fname <- basename(f)
    
    parts <- strsplit(fname, "\\.")[[1]]
    stopifnot(length(parts) >= 4)
    
    nparts <- length(parts)
    chr <- as.integer(parts[nparts - 2])
    cM <- as.integer(parts[nparts - 1])

    if(verbose > 1) cat(" * solarMIBD: file name ", fname, ", chr ", chr, ", cM ", cM, "\n", sep = "")
    
    stopifnot(!is.na(chr))
    stopifnot(!is.na(cM))
    stopifnot(chr > 0)
    stopifnot(chr < 23)    
    
    list(file = f, chr = chr, cM = cM)
  })
  
  ### pass 2: fitler by chr
  if(!missing(chr)) {
    chr.par <- chr
    
    ind <- laply(out, function(x) x$chr %in% chr.par)
    
    out <- out[ind]
    stopifnot(length(out) > 0)
  }
  
  ### pass 3.1: filter by `nmibd`
  if(!missing(nmibd)) {
    stopifnot(length(out) >= nmibd)
    
    out <- out[1:nmibd]
  }
  
  ### pass 3.2: read MIBD matrices
  out <- llply(out, function(x) {
   if(verbose > 0) cat(" * solarMIBD: reading", x$file, "\n")
   
    mf <- read_mibd_csv_gz(x$file)
    mat <- mf2mat(mf)
    mat <- Matrix::Matrix(mat)
    
    c(x, list(mibd = mat))
  }, .parallel = parallel)
  
  ### return
  return(out)  
}

get_info_mibd <- function(mibddir)
{
  stopifnot(file.exists(mibddir))

  files <- list.files(mibddir, pattern = ".gz", full.names = TRUE)
  stopifnot(length(files) > 0)
  file1 <- files[1]
  
  ### `mibddir.format`
  mibddir.format <- ifelse(grepl(",", readLines(file1, n = 1)), "csv", "pedindex")
  
  ### `mibddir.ids`  
  mibddir.ids <- NULL
  if(mibddir.format == "csv") {
    tab <- read.table(file1, sep = ",", header = TRUE, colClasses = "character")
    stopifnot(all(c("id1", "id2") %in% names(tab)))
    
    mibddir.ids <- unique(with(tab, c(id1, id2)))
  }
  
  out <- list(mibddir.format = mibddir.format, mibddir.ids = mibddir.ids)
  
  ### return
  return(out)
}

convert_mibd <- function(indir, outdir, pedindex.out, verbose = 1)
{
  stopifnot(file.exists(indir))
  stopifnot(dir.create(outdir))
  
  ### infiles
  infiles <- list.files(indir, pattern = ".gz", full.names = TRUE)
  
  ### read pedindex
  pf <- read_pedindex(pedindex.out)
  N <- nrow(pf)
  
  for(i in 1:length(infiles)) {
    f <- infiles[i]
    
    if(verbose) {
      cat(" *", i, "/", length(infiles), "file", f, "\n")
    }
    
    mf <- read_mibd_gz(f)
    N.diag <- with(mf, sum(IBDID1 == IBDID2))
    stopifnot(N == N.diag)
    
    mf <- kf_match_pedindex(mf, pf)

    mf2 <- subset(mf, select = c("ID1", "ID2", "matrix1", "matrix2"))
    mf2 <- rename(mf2, c(ID1 = "id1", ID2 = "id2"))
    
    # order
    ord <- with(mf2, order(as.integer(id1), as.integer(id2)))
    mf2 <- mf2[ord, ]
    
    ### write file
    of <- file.path(outdir, basename(f))
    ret <- write.table(mf2, gzfile(of), quote = FALSE,
      row.names = FALSE, col.names = TRUE, sep = ",")
  }
  
  return(invisible())
}

#
read_mibd_gz <- function(mibd.gz)
{
  mf <- read.table(gzfile(mibd.gz), colClasses = c("character", "character", "numeric", "numeric"))
  names(mf) <- c("IBDID1", "IBDID2", "matrix1", "matrix2")
###
#         1          1   1.000000   1.000000
#         2          2   1.000000   1.000000
#         3          1   0.500000   0.000000
#         3          2   0.500000   0.000000

  stopifnot(ncol(mf) == 4)
   
  return(mf)
}

read_mibd_csv_gz <- function(mibd.gz)
{
  mf <- read.table(gzfile(mibd.gz), sep = ",", header = TRUE,
    colClasses = c("character", "character", "numeric", "numeric"))
  names(mf) <- c("ID1", "ID2", "matrix1", "matrix2")
#
#    ID1   ID2 matrix1 matrix2
#1 01101 01101     1.0       1
#2 01102 01102     1.0       1
#3 01202 01101     0.5       0
#
  stopifnot(ncol(mf) == 4)
   
  return(mf)
}

mf2mat <- function(mf)
{
  stopifnot(class(mf) == "data.frame")
  stopifnot(all(c("ID1", "ID2", "matrix1") %in% names(mf)))
  
  ids <- unique(c(mf$ID1, mf$ID2))
  N <- length(ids)
  
  mat <- matrix(0, nrow = N, ncol = N)
  rownames(mat) <- ids
  colnames(mat) <- ids
  
  for(i in 1:nrow(mf)) {
    mat[mf$ID1[i], mf$ID2[i]] <- mf$matrix1[i]
    mat[mf$ID2[i], mf$ID1[i]] <- mf$matrix1[i]
  }
  
  return(mat)
}
  
#----------------------------------
# Read LOD functions
#----------------------------------

read_multipoint_lod <- function(dir, num.traits, mode = "none")  
{
  stopifnot(!missing(num.traits)) 
 
  switch(mode,
    "none" = {
      switch(as.character(num.traits),
        "1" = read_multipoint_lod_univar(dir),
        "2" = read_multipoint_lod_bivar(dir),
        stop("error in switch by `num.traits`"))
    },
    "gxed" = read_multipoint_lod_gxed(dir),
    stop("error in switch by `mode`"))
}

read_multipoint_lod_gxed <- function (dir)
{
  out <- list()
  
  multipoint.files <- list.files(dir, "multipoint[1-9].out", full.names = TRUE)
  num.passes <- length(multipoint.files)
  stopifnot(num.passes <= 1)
  
  names.tab <- list("1" = c("Model", "LOD", "Loglike", "xgsd", "ygsd", "rhog", "xqsd1", "yqsd1", "rhoq1"))
  ncol.tab <- list("1" = 12)
  col.tab <- list("1" = c(2, 4, 5:12))
  colnames.tab <- list("1" = c("chr", "pos", "LOD", "Loglike", "xgsd", "ygsd", "rhog", "xqsd1", "yqsd1", "rhoq1"))

### pass 1
#       Model            LOD        Loglike      xgsd      ygsd      rhog     xqsd1     yqsd1     rhoq1  
#------------------- -----------  -----------  --------  --------  -------- --------- --------- ---------
#chrom 02  loc     7      3.3134      673.687  0.031174  0.074905  1.000000  0.083059  0.040770  1.000000 
#chrom 02  loc     9      3.2409      673.520  0.032504  0.074223  1.000000  0.080781  0.042219  1.000000 
#...

### pass 2
# ...

  for(i in 1:num.passes) {
    f <- multipoint.files[i]
    
    tnames <- names.tab[[i]]
    tncol <- ncol.tab[[i]]
    tcol <- col.tab[[i]]
    tcolnames <- colnames.tab[[i]]
    
    # names
    names <- unlist(strsplit(readLines(f, n = 1), "\\s+"))
    names <- names[names != ""]
    stopifnot(all(names == tnames))
    
    # table 
    tab <- read.table(f, skip = 2)
    stopifnot(ncol(tab) == tncol)
    tab <- tab[, tcol]
    colnames(tab) <- tcolnames
    
    out <- c(out, list(tab))
  }

  if(num.passes > 0) {
    names(out) <- paste("df", 1:num.passes, sep = "")
    names(out)[1] <- "df"
  }
  
  out$num.passes <- num.passes
  
  return(out)  
}

read_multipoint_lod_bivar <- function(dir)  
{
  out <- list()
  
  multipoint.files <- list.files(dir, "multipoint[1-9].out", full.names = TRUE)
  num.passes <- length(multipoint.files)
  stopifnot(num.passes <= 1)
  
  names.tab <- list("1" = c("Model", "LOD", "Loglike", "H2r_trait1", "H2r_trait2", 
    "RhoG", "H2q1_trait1", "H2q1_trait2", "RhoQ1"))
  ncol.tab <- list("1" = 12)
  col.tab <- list("1" = c(2, 4, 5:12))
  colnames.tab <- list("1" = c("chr", "pos", "LOD", "Loglike", "H2r_trait1", 
    "H2r_trait2", "RhoG", "H2q1_trait1", "H2q1_trait2", "RhoQ1"))

### pass 1
#       Model            LOD        Loglike   H2r(trait1) H2r(trait2)    RhoG   H2q1(trait1) H2q1(trait2)   RhoQ1   
#------------------- -----------  ----------- ----------- -----------  -------- ------------ ------------ --------- 
#chrom 05  loc     0      9.7291    -2612.669    0.485964    0.513254  0.930689     0.386615     0.236398  1.000000  
#chrom 05  loc     1     12.0168    -2607.213    0.364941    0.421572  0.915142     0.503510     0.324686  1.000000  
#chrom 05  loc     2     13.7087    -2603.198    0.285072    0.356037  0.898079     0.579351     0.385922  1.000000 
#  ...

  for(i in 1:num.passes) {
    f <- multipoint.files[i]
    
    tnames <- names.tab[[i]]
    tncol <- ncol.tab[[i]]
    tcol <- col.tab[[i]]
    tcolnames <- colnames.tab[[i]]
    
    # names
    names <- unlist(strsplit(readLines(f, n = 1), "\\s+"))
    names <- names[names != ""]
    stopifnot(all(names[1:3] == tnames[1:3])) # just first three names, as further naming depend on particular names of traits 
    
    # table 
    tab <- read.table(f, skip = 2)
    stopifnot(ncol(tab) == tncol)
    tab <- tab[, tcol]
    colnames(tab) <- tcolnames
    
    out <- c(out, list(tab))
  }

  if(num.passes > 0) {
    names(out) <- paste("df", 1:num.passes)
    names(out)[1] <- "df"
  }
  
  out$num.passes <- num.passes
  
  return(out)
}

read_multipoint_lod_univar <- function(dir)  
{
  out <- list()
  
  multipoint.files <- list.files(dir, "multipoint[1-9].out", full.names = TRUE)
  num.passes <- length(multipoint.files)
  stopifnot(num.passes <= 2)
  
  names.tab <- list("1" = c("Model", "LOD", "Loglike", "H2r", "H2q1"),
    "2" = c("Model", "LOD", "Loglike", "H2r", "H2q1", "H2q2"))
  ncol.tab <- list("1" = 8, "2" = 9)
  col.tab <- list("1" = c(2, 4, 5:8), "2" = c(2, 4, 5:9))
  colnames.tab <- list("1" = c("chr", "pos", "LOD", "Loglike", "H2r", "H2q1"),
    "2" = c("chr", "pos", "LOD", "Loglike", "H2r", "H2q1", "H2q2"))

### pass 1
#             Model         LOD       Loglike       H2r      H2q1  
#  -------------------   ---------  -----------  --------  --------
#  chrom 05  loc     0     11.0544    -1466.123  0.492049  0.386445
#  chrom 05  loc     1     13.2765    -1461.007  0.375424  0.499153
#  chrom 05  loc     2     14.9359    -1457.186  0.298530  0.572836
#  ...

### pass 2
#            Model         LOD       Loglike       H2r      H2q1      H2q2  
# -------------------   ---------  -----------  --------  --------  --------
# chrom 05  loc     0      0.7355    -1455.422  0.285594  0.439268  0.145092
# chrom 05  loc     1      0.5086    -1455.945  0.284433  0.404778  0.180871
# chrom 05  loc     2      0.3176    -1456.384  0.284086  0.301798  0.284157
# ...

  for(i in 1:num.passes) {
    f <- multipoint.files[i]
    
    tnames <- names.tab[[i]]
    tncol <- ncol.tab[[i]]
    tcol <- col.tab[[i]]
    tcolnames <- colnames.tab[[i]]
    
    # names
    names <- unlist(strsplit(readLines(f, n = 1), "\\s+"))
    names <- names[names != ""]
    stopifnot(all(names == tnames))
    
    # table 
    tab <- read.table(f, skip = 2)
    stopifnot(ncol(tab) == tncol)
    tab <- tab[, tcol]
    colnames(tab) <- tcolnames
    
    out <- c(out, list(tab))
  }

  if(num.passes > 0) {
    names(out) <- paste("df", 1:num.passes, sep = "")
    names(out)[1] <- "df"
  }
  
  out$num.passes <- num.passes
  
  return(out)
}

