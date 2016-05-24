
### run 1: data set from multic package (copy files)
dir <- tempfile(pattern = "solar-multipoint")

dir.create(dir)
files <- list.files(package.file("extdata", "solarOutput", package = "solarius"), full.names = TRUE)
file.copy(files, dir, recursive = TRUE)
mibddir <- "solarMibds"

### run solar 1
cmd <- c(
  "pedigree load simulated.ped", 
  "phen load simulated.phen",
  "model new", "trait trait1", "polygenic")
ret <- solar(cmd, dir)

### run solar 2
cmd <- c("load model trait1/null0.mod",
  paste("mibddir", mibddir), 
  "chromosome all", "interval 5", 
  "multipoint -overwrite")

ret <- solar(cmd, dir)

print(ret)

### clean
unlink(dir, recursive = TRUE)

### run 2: data set from multic package (data.frame + no order)
dir <- tempfile(pattern = "solar-multipoint")

dat <- loadMulticPhen()
mibddir <- package.file("extdata", "solarOutput", "solarMibds", package = "solarius")  

df2solar(dat, dir, sort.ped = FALSE) 

### run solar 1
cmd <- c(
  "pedigree load dat.ped", 
  "phen load dat.phe",
  "model new", "trait trait1", "polygenic")
ret <- solar(cmd, dir)

### run solar 2
cmd <- c("load model trait1/null0.mod",
  paste("mibddir", mibddir), 
  "chromosome all", "interval 5", 
  "multipoint -overwrite")

ret <- solar(cmd, dir)

print(ret)

### clean
unlink(dir, recursive = TRUE)

### run 3: data set from multic package (data.frame + order)
dir <- tempfile(pattern = "solar-multipoint")

dat <- loadMulticPhen()
mibddir <- package.file("extdata", "solarOutput", "solarMibds", package = "solarius")  

df2solar(dat, dir, sort.ped = TRUE) 

### run solar 1
cmd <- c(
  "pedigree load dat.ped", 
  "phen load dat.phe",
  "model new", "trait trait1", "polygenic")
ret <- solar(cmd, dir)

### run solar 2
cmd <- c("load model trait1/null0.mod",
  paste("mibddir", mibddir), 
  "chromosome all", "interval 5", 
  "multipoint -overwrite")

ret <- solar(cmd, dir)

print(ret)

### clean
unlink(dir, recursive = TRUE)
