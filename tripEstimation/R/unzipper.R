"unzipper" <-
function(px) {


## px is a list of pimg.lists
epoch <- ISOdatetime(1970, 1, 1, 0, 0, 0, tz = "GMT")
tms <- unlist(lapply(px, function(x) attr(x, "times"))) 
Zs <- unlist(lapply(px, function(x) attr(x, "Z")))
if (sum(Zs) > 0 && !all(Zs)) stop("mixture of X and Z pimg.lists not supported")
  tms <- 12*3600*(round(unclass(tms)/(12*3600))) + epoch

xb <- px[[1]]$xbound
yb <- px[[1]]$ybound
#unique.times <- unique(sort(unique(tms)) + epoch)


p <- pimg.list(tms, xlim = c(xb[1], xb[2]),  ylim = c(yb[1], yb[2]), img.dim = c(xb[3], yb[3]), Z = all(Zs))


l <- list()


for (ix in px) l <- c(l, ix)

for (ip in 1:length(l)) p[[ip]] <- l[[ip]]

names(p) <- names(tms)
p
#attr(l, "times") <- tms
## px is all the pimgs matching tms

}

