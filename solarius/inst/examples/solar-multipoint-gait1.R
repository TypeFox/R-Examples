#library(gait)
load_all("~/git/ugcd/gait")

library(plyr)

dir <- tempfile(pattern = "solar-multipoint")
dir <- "solar"

### gait phenotypes
traits0 <- gait1.traits.cascade()
traits <- paste0("tr1_", traits0)

dat <- gait1.phen(transform = "tr1", traits = traits[1])

ids <- gait1.id.ibd()
dat <- mutate(dat, 
  MIBD = ifelse(ID %in% ids, 1, NA))

### linkage 
mibddir <- gait1.mibddir()

### expclude FAMIDs
#fams <- c("02", "17", "18")
#dat <- subset(dat, !(FAMID %in% fams))

### export data
df2solar(dat, dir) 

#f <- list.files("/home/andrey/Datasets/GAIT1/mibdPedigree/", full.names = TRUE)
#ret <- file.copy(f, dir)

### run solar 1
cmd <- c(
  "pedigree load dat.ped", 
  "phen load dat.phe",
  "model new", "trait tr1_FVII", 
  #"covariate MIBD()", 
  "polygenic")
ret <- solar(cmd, dir)

### run solar 2
cmd <- c("load model tr1_FVII/null0.mod",
  paste("mibddir", mibddir), 
  "chromosome 22", "interval 50", 
  "multipoint -overwrite")

ret <- solar(cmd, dir)

### clean
#unlink(dir, recursive = TRUE)

ibdids <- c("27", "28", "30", "31", "42", "43", "364", "365", "393")  
pf <- read_pedindex("solar/pedindex.out")
ids <- with(pf, ID[IBDID %in% ibdids])
stopifnot(length(ids) == length(ibdids))

mf.extra <- data.frame(ID1 = ids, ID2 = ids, matrix1 = 1, matrix2 = 1)

files <- head(list.files(mibddir, full.names = TRUE))
for(f in files) {
  mf <- read_mibd_csv_gz(f)
  mf <- rbind(mf, mf.extra)
  
  of <- file.path(dir, basename(f))
  
  z <- gzfile(of, "w")  # compressed file
  #cat("TITLE extra line", "2 3 5 7", "", "11 13 17", file = zz, sep = "\n")
  write.table(mf, z, sep = ",", quote = FALSE,
    row.names = FALSE, col.names = TRUE)
  close(z)
}

### run solar 3
cmd <- c("load model tr1_FVII/null0.mod",
  paste("mibddir", "."), 
  "chromosome 10", "interval 1", 
  "multipoint -overwrite")

ret <- solar(cmd, dir)



