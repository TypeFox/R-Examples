"EnterMeta" <-
function () 
 { 
d <- matrix(,ncol=7)
d <- data.frame(d)
names(d) <-c("study", "Rxy", "n", "Rxx", "Ryy", "u", "moderator")
d$study <- as.factor(d$study)
d$Rxy <- as.numeric(d$Rxy)
d$n <- as.numeric(d$n)
d$Rxx <- as.numeric(d$Rxx)
d$Ryy <- as.numeric(d$Ryy)
d$u <- as.numeric(d$u)
d$moderator <- as.factor(d$moderator)
meta <- edit(d)
}

