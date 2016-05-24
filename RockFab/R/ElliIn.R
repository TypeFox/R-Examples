ElliIn <-
function(elli.files){
  #For each file extract required ellipsoid data
  for(j in 1:length(elli.files)){
    elli.lines <- readLines(con = elli.files[j], warn = FALSE)
x.mag <- as.numeric(unlist(strsplit(elli.lines[18], split = '\t'))[2])
y.mag <- as.numeric(unlist(strsplit(elli.lines[18], split = '\t'))[3])
z.mag <- as.numeric(unlist(strsplit(elli.lines[18], split = '\t'))[4])
    x.az <- as.numeric(unlist(strsplit(elli.lines[21], split = '\t'))[2])
y.az <- as.numeric(unlist(strsplit(elli.lines[21], split = '\t'))[3])
z.az <- as.numeric(unlist(strsplit(elli.lines[21], split = '\t'))[4])
    x.in <- as.numeric(unlist(strsplit(elli.lines[22], split = '\t'))[2])
y.in <- as.numeric(unlist(strsplit(elli.lines[22], split = '\t'))[3])
z.in <- as.numeric(unlist(strsplit(elli.lines[22], split = '\t'))[4])
xy.strike <- as.numeric(unlist(strsplit(elli.lines[24], split = '\t'))[2])
xy.dip <- as.numeric(unlist(strsplit(elli.lines[24], split = '\t'))[3])
x.rake <- as.numeric(unlist(strsplit(elli.lines[24], split = '\t'))[5])
n.samp <- as.numeric(unlist(strsplit(elli.lines[6], split = '\t'))[2])
f.index <- as.numeric(unlist(strsplit(elli.lines[7], split = '\t'))[2])
oss <- (sqrt((log(x.mag) - log(y.mag))^2 + (log(y.mag) - log(z.mag))^2 + (log(z.mag) - log(x.mag))^2) / sqrt(3))
lp <- (2 * log(y.mag) - log(x.mag) - log(z.mag)) / (log(x.mag) - log(z.mag))
#Create data frame object on first iteration
if(j == 1){
  my.data <- data.frame(
    id = j,
    elli.file = elli.files[j],
    x.mag = x.mag,
y.mag = y.mag,
z.mag = z.mag,
x.az = x.az,
y.az = y.az,
z.az = z.az,
x.in = x.in,
y.in = y.in,
z.in = z.in,
xy.strike = xy.strike,
xy.dip = xy.dip,
x.rake = x.rake,
n.samp = n.samp,
f.index = f.index,
oss = oss,
lp = lp, stringsAsFactors = FALSE)
#Populate data frame rows for additional files
}else(my.data[j,] <- list(j, elli.files[j], x.mag, y.mag, z.mag, x.az, y.az, z.az, x.in, y.in, z.in, xy.strike, xy.dip, x.rake, n.samp, f.index, oss, lp))
    remove(elli.lines)
  }
  return(my.data)
}
