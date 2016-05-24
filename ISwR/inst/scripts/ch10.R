age <- subset(juul, age >= 10 & age <= 16)$age
range(age)
agegr <- cut(age, seq(10,16,2), right=F, include.lowest=T)
length(age)
table(agegr)
agegr2 <- cut(age, seq(10,16,2), right=F)
table(agegr2)
q <- quantile(age, c(0, .25, .50, .75, 1))
q
ageQ <- cut(age, q, include.lowest=T)
table(ageQ)
levels(ageQ) <- c("1st", "2nd", "3rd", "4th")
levels(agegr) <- c("10-11", "12-13", "14-15") 
pain <- c(0,3,2,2,1)
fpain <- factor(pain,levels=0:3, 
       labels=c("none","mild","medium","severe"))
text.pain <-  c("none","severe", "medium", "medium", "mild") 
factor(text.pain)
ftpain <- factor(text.pain)
ftpain2 <- factor(ftpain,
                  levels=c("none", "mild", "medium", "severe"))
ftpain3 <- ftpain2
levels(ftpain3) <- list(
        none="none",
        intermediate=c("mild","medium"),
        severe="severe")
ftpain3
ftpain4 <- ftpain2
levels(ftpain4) <- c("none","intermediate","intermediate","severe")
ftpain4
stroke <- read.csv2(
  system.file("rawdata","stroke.csv", package="ISwR"),
  na.strings=".") 
names(stroke) <- tolower(names(stroke))
head(stroke)
stroke <- transform(stroke,
   died = as.Date(died, format="%d.%m.%Y"),
   dstr = as.Date(dstr, format="%d.%m.%Y"))
summary(stroke$died) 
summary(stroke$dstr)
summary(stroke$died - stroke$dstr)
head(stroke$died - stroke$dstr)
o <- options(width=60) # minor cheat for visual purposes
stroke <- transform(stroke,
  end = pmin(died,  as.Date("1996-1-1"), na.rm = T),
  dead = !is.na(died) & died < as.Date("1996-1-1"))
head(stroke)
options(o); rm(o)
stroke <- transform(stroke,
  obstime = as.numeric(end - dstr, units="days")/365.25)  
rawstroke <- read.csv2(
  system.file("rawdata","stroke.csv", package="ISwR"),
  na.strings=".") 
ix <- c("DSTR", "DIED")
rawstroke[ix] <- lapply(rawstroke[ix], 
                        as.Date, format="%d.%m.%Y") 
head(rawstroke) 
ix <- 6:9
rawstroke[ix] <- lapply(rawstroke[ix], 
                        factor, levels=0:1, labels=c("No","Yes"))   
strokesub <- ISwR::stroke[1:10,2:3]
strokesub
strokesub <- transform(strokesub, 
  event = !is.na(died))
strokesub <- transform(strokesub, 
 obstime = ifelse(event, died-dstr, as.Date("1996-1-1") - dstr))
strokesub
juulgrl <- subset(juul, sex==2, select=-c(testvol,sex))
juulboy <- subset(juul, sex==1, select=-c(menarche,sex))
juulgrl$sex <- factor("F")
juulgrl$testvol <- NA
juulboy$sex <- factor("M")
juulboy$menarche <- NA
juulall <- rbind(juulboy, juulgrl)
names(juulall)
levels(juulall$sex)
head(nickel)
head(ewrates)
nickel <- transform(nickel, 
  agr = trunc(agein/5)*5, 
  ygr = trunc((dob+agein-1)/5)*5+1)  
mrg <- merge(nickel, ewrates, 
  by.x=c("agr","ygr"), by.y=c("age","year"))
head(mrg,10)
head(alkfos)
a2 <- alkfos
names(a2) <- sub("c", "c.", names(a2))
names(a2)
a.long <- reshape(a2, varying=2:8, direction="long")
head(a.long)
tail(a.long)
o <- with(a.long, order(id, time))
head(a.long[o,], 10)
a.long2 <- na.omit(a.long)
attr(a.long2, "reshapeLong") <- NULL
a.wide2 <- reshape(a.long2, direction="wide", v.names="c",
                 idvar="id", timevar="time")
head(a.wide2)
l <- split(a.long$c, a.long$id)
l[1:3] 
l2 <- lapply(l, function(x) x / x[1])
a.long$c.adj <- unsplit(l2, a.long$id)
subset(a.long, id==1)
a.long$c.adj <- ave(a.long$c, a.long$id, 
    FUN = function(x) x / x[1])
all.equal(unsplit(l2, a.long$id),  a.long$c.adj)
l <- split(a.long, a.long$id)
l2 <- lapply(l, transform, c.adj = c / c[1])
a.long2 <- unsplit(l2, a.long$id)
all.equal(a.long2$c.adj,  a.long$c.adj)
head(nickel)
entry <- pmax(nickel$agein, 60)
exit <- pmin(nickel$ageout, 65)
valid <- (entry < exit) 
entry <- entry[valid]
exit  <- exit[valid]
cens <- (nickel$ageout[valid] > 65)  
nickel60 <- nickel[valid,]
nickel60$icd[cens] <- 0
nickel60$agein <- entry
nickel60$ageout <- exit
nickel60$agr <- 60
nickel60$ygr <- with(nickel60, trunc((dob+agein-1)/5)*5+1)
head(nickel60)
trim <- function(start)
{
  end   <- start + 5
  entry <- pmax(nickel$agein, start)
  exit  <- pmin(nickel$ageout, end)
  valid <- (entry < exit) 
  cens  <- (nickel$ageout[valid] > end)  
  result <- nickel[valid,]
  result$icd[cens] <- 0
  result$agein <- entry[valid]
  result$ageout <- exit[valid]
  result$agr <- start
  result$ygr <- with(result, trunc((dob+agein-1)/5)*5+1)
  result
}
head(trim(60))
nickel.expand <- do.call("rbind", lapply(seq(20,95,5), trim))
head(nickel.expand)
subset(nickel.expand, id==4)  
nickel.expand <- merge(nickel.expand, ewrates, 
  by.x=c("agr","ygr"), by.y=c("age","year"))
head(nickel.expand)
all.equal(nickel.expand, ISwR::nickel.expand)
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
