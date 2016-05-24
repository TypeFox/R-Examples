library(multic)

### par
multic.dir <- "multic"

### dat
dat.dir <- system.file("inst", "extdata", "solarOutput", package = "solarius")
if(!file.exists(dat.dir)) {
  dat.dir <- system.file("extdata", "solarOutput", package = "solarius")
}
stopifnot(file.exists(dat.dir))

ped <- read.table(file.path(dat.dir, "simulated.ped"), header = TRUE, sep = ",")
phen <- read.table(file.path(dat.dir, "simulated.phen"), header = TRUE, sep = ",")
dat <- merge(ped, phen)

dat <- dat[order(dat$famid, dat$id), ] ## data needs to be sorted by famid

dat30 <- subset(dat, famid < 30)

### export into `multic` format
solar2multic(phi2 = file.path(dat.dir, "phi2.gz"), 
  pedigree.file = file.path(dat.dir, "simulated.ped"),
  pedindex.out = file.path(dat.dir, "pedindex.out"),
  pedindex.cde = file.path(dat.dir, "pedindex.cde"),
  ibd.directory = file.path(dat.dir, "solarMibds"), 
  output.directory = multic.dir)
  
### run `multic`
mod1 <- multic(trait1 ~ sex + age, data = dat[with(dat, order(famid)), ], subset = (famid < 30),
  famid, id, fa, mo, sex,
  #mloci.out = file.path(dirInput, "mloci.out"),
  share.out = file.path(multic.dir, "share.out"))
  
mod2 <- multic(trait1 ~ sex + age, data = dat[with(dat, order(famid)), ], subset = (famid < 30),
  famid, id, fa, mo, sex,
  mloci.out = file.path(multic.dir, "mloci.out"),
  share.out = file.path(multic.dir, "share.out"))

### plots
p2 <- ggplot(mod2$log.liks, aes(distance, lod.score)) + geom_point() + geom_line()
