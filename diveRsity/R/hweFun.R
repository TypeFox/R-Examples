#' fisher's exact tests for HWE
#' Kevin Keenan, 2014
hweFun <- function(infile, mcRep = 2000){
  #infile <- "./test_files/mono.gen"
  #rgp <- diveRsity:::rgp
  #dat <- rgp(infile)
  #mcRep = 10000
  dat <- rgp(infile)
  #locus chisq calc
  chiLoc <- lapply(dat$genos, function(x){
    apply(x, 2, function(y){
      al1 <- na.omit(y[,1])
      al2 <- na.omit(y[,2])
      if(identical(al1, al2)){
        list(p = NA, chisq = NA, method = NA)
      } else {
        lev <- unique(c(al1, al2))
        chi_tab <- table(factor(al1, levels = lev), factor(al2, levels = lev))
        chi_tab <- 0.5*(chi_tab + t(chi_tab))
        res <- chisq.test(chi_tab, simulate.p.value = T, B = mcRep)
        list(p = res$p.value, chisq = res$statistic, method = res$method)
      }
    })
  })
  pval <- lapply(chiLoc, function(x){
    sapply(x, "[[", "p")
  })
  chiAll <- lapply(pval, function(x){
    df <- 2*length(x)
    chi <- -2*sum(log(x), na.rm = TRUE)
    p <- pchisq(chi, df = df, lower.tail = FALSE)
    list(p = p, chisq = chi)
  })
  list(locus = chiLoc, multilocus = chiAll)
}