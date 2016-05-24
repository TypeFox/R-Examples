simPedHS <- function(s, d, n, uniqueDname = TRUE, prefix = NULL){
  if(n<2){ stop("must have more than 1 offspring per family (n > or = 2)")
   }
  sires <- paste("s", seq(1, s, 1), sep="")
  if(uniqueDname) {
    dams <- paste("d", seq(1, (s*d), 1), sep="")
    } else{
      dams <- rep(paste("d", seq(1,d, 1), sep=""), s)
      }
  offspring <- paste("o", seq(1, (s*d*n), 1), sep="")

  if(is.null(prefix)){
    ped <- data.frame(id = c(sires, dams, offspring),
	dam = c(rep(NA, length(sires)), rep(NA, length(dams)), rep(dams, each=n)),
	sire = c(rep(NA, length(sires)), rep(NA, length(dams)), rep(sires, each=d*n)),
	sex = c(rep("M", length(sires)), rep("F", length(dams)), rep(c("M","F"), length.out = length(offspring))))
  } else{
       p <- as.character(prefix)
       ped <- data.frame(id = paste0(p, c(sires, dams, offspring)),
	dam = c(rep(NA, length(sires)), rep(NA, length(dams)), rep(paste0(p, dams), each=n)),
	sire = c(rep(NA, length(sires)), rep(NA, length(dams)), rep(paste0(p, sires), each=d*n)),
	sex = c(rep("M", length(sires)), rep("F", length(dams)), rep(c("M","F"), length.out = length(offspring))))     
    }
ped
}

