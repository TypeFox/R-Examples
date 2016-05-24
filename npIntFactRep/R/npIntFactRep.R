#' @title Nonparametric Interaction Tests for Factorial Designs with Repeated Measures
#' @description Nonparametric aligned rank tests for interactions in two-way factorial designs with repeated measures (e.g., Beasley & Zumbo, 2006; Higgings & Tashtoush, 1994). Five ANOVAs are reported: (1) PARAMETRIC on the original data, (2) on the interaction alignments as a CHECK, (3) on the aligned REGULAR, (4) FRIEDMAN, and (5) KOCH ranks. In the rank tests, only the resulting values for the interaction are relevant.
#' @usage npIntFactRep(dat)
#' @references
#' Higgins, J.J., & Tashtoush, S. (1994). "An aligned rank transform test for interaction". Nonlinear World, 1, 201-211.
#' Beasley, T.M., & Zumbo, B.D. (2009). "Aligned rank tests for interactions in split-plot designs: Distributional assumptions and stochastic homogeneity". Journal of Modern Applied Statistical Methods, 8, 16-50.
#' @details Returns five ANOVA tables using the 'ezANOVA' function (from the 'ez' package by Michael A. Lawrence), which is well-suited for factorial designs with repeated measures. The Higgins & Tashtoush formula for repeated measures (1994, p. 208) is used for the interaction alignments. The ANOVAs include F- and p-values, along with generalized eta squared (ges) effect size statistics and the sphericity test. Type II sum of squares are calculated, which are appropriate for the alignments check with unbalanced designs. The first ANOVA is the PARAMETRIC one on the original data. The second one is on the interaction alignments as a CHECK (the p-values for 'group' and for 'rep' should both be = 1). The three last ANOVAs are aligned rank tests: one on the aligned REGULAR, one on the FRIEDMAN, and one on the KOCH ranks. In these ranks tests, ONLY the resulting values for the INTERACTION group:rep are relevant
#' @author Jos Feys
#' @param dat The name of the R data set in wide (one row per subject) format. This set should NOT have any missing values (NA). Missing data should be replaced (e.g., with 'mi', 'mice' or 'Amelia') before running. The data set MUST have subj as name (header) for the (numeric or character values of the) units of observations (subjects/blocks) factor and MUST have group as name for the (numeric or character values of the) between factor. Other variables must be excluded, because the number of levels of the repeated measures (within) factor is calculated based upon the remaining columns after eliminating the columns group and subj. The repeated measures factor is named rep.
#' @examples
#' \dontrun{
#' dat1 <- read.csv (file="c:/R/wide.csv", head=T)
#' npIntFactRep(dat1)
#' }
#' # example with character values for both group and subj variables (columns)
#' dat2 <- read.table(header = TRUE, text = "
#' subj group x1 x2 x3 x4
#'  p1 a 1 2 3 4
#'  p2 a 2 2 3 3
#'  p3 a 1 3 3 4
#'  p4 b 8 6 4 2
#'  p5 b 6 6 4 4
#'  p6 b 8 6 6 2
#'  p7 c 4 4 4 3
#'  p8 c 3 3 3 4
#'  p9 c 3 4 4 2
#'  ")
#' npIntFactRep(dat2)
#' # example with numeric values for the subj variable and charater values for the group variable
#' dat3 <- read.table(header = TRUE, text = "
#' subj group 1 2 3 4 5
#' 1 a1 1 2 3 4 5
#' 2 a1 2 2 3 3 3
#' 3 a1 1 3 3 4 2
#' 4 a2 8 6 4 2 1
#' 5 a2 6 6 4 4 4
#' 6 a2 8 6 6 2 3
#' ")
#' npIntFactRep(dat3)
#' @import ez
#' @import plyr
#' @import stats
#' @return 5 ANOVA summary tables
#' @export

npIntFactRep <- function (dat){
  dep <- NULL # Setting the variables to NULL first
  colnames(dat) <- tolower(colnames(dat))

  l <- ncol(dat)
  dat_long <- reshape(dat, direction="long",
                      varying=list(names(dat)[3:l]),
                      v.names="dep",timevar = "rep",
                      idvar=c("subj","group"))

  #sort by subj, group

  dat_long <- dat_long[order(dat_long$subj,dat_long$group,dat_long$rep),]
  # dat_long

  ## Convert variables to factor
  dat_long <- within(dat_long, {
    group <- factor(group)
    rep <- factor(rep)
    subj <- factor(subj)
  })


  Parametric <- ezANOVA(data=dat_long, dv = .(dep), wid=.(subj), between = .(group),
                        within=.(rep), type=2, detailed=FALSE)
  # Parametric

  group <- dat["group"]
  subj  <- dat["subj"]
  # remove columns "subj" & "group"
  # only  rep. meas. columns are left
  # one does not need the know the number
  colsdel <- c("group","subj")
  y <- dat[,!(names (dat) %in% colsdel)]

  # number of columns (rep.meas.) & observations (subj's)
  k <- ncol(y)
  n <- nrow(y)
  # initializations (probably not necessary, but ...)
  ad <- matrix(0,n,k)
  ar <- matrix(0,n,k)
  fr <- matrix(0,n,k)
  q  <- matrix(0,n,k)

  # Higgins & Tashtoush formula (1994, p. 108)
  rmean <- colMeans(y)
  rmean <- matrix(rep(rmean, each=n), ncol=k)
  pmean <- rowMeans(y)
  gmean <- mean (pmean)
  #aligned data
  ad <- (y - rmean - pmean) + gmean
  # aligned data
  ad_full <- cbind(subj, group, ad)
  # ad_full
  l <- ncol(ad_full)
  ad_long <- reshape(ad_full, direction="long",
                     varying=list(names(ad_full)[3:l]),
                     v.names="dep", timevar="rep",
                     idvar=c("subj","group"))
  ## Convert variables to factor
  ad_long <- within(ad_long, {
    group <- factor(group)
    rep <- factor(rep)
    subj <- factor(subj)
  })

  AlignCheck <- ezANOVA(data=ad_long, dv = .(dep), wid=.(subj), between = .(group),
                        within=.(rep), type=2, detailed=FALSE)


  ar <- matrix(rank(ad), ncol=k) ################ gaf fout ties='average'
  # voor adjustment
  # ar
  # div <- ((n*k)+1)
  # ar <- ar/ div
  # na adjustment
  # ar
  ar_full <- cbind(subj, group, ar)

  l <- ncol(ar_full)
  ar_long <- reshape(ar_full, direction="long",
                     varying=list(names(ar_full)[3:l]),
                     v.names="dep", timevar="rep",
                     idvar=c("subj","group"))

  ## Convert variables to factor
  ar_long <- within(ar_long, {
    group <- factor(group)
    rep <- factor(rep)
    subj <- factor(subj)
  })


  RegularRanks <- ezANOVA(data=ar_long, dv = .(dep), wid=.(subj), between = .(group),
                          within=.(rep), type=2, detailed=FALSE)


  #FRIEDMAN
  for(i in 1:n){
    fr[i,] <- matrix(rank(ad[i,]), ncol=k)
  }
  #before adjustment
  #fr
  # part1 = (k+1)/2
  # part2 = ((k**2)-1)/12
  # fr <- (fr-part1)/part2
  #after adjustment
  #fr
  fr_full <- cbind(subj, group, fr)

  l <- ncol(fr_full)
  fr_long <- reshape(fr_full, direction="long",
                     varying=list(names(fr_full)[3:l]),
                     v.names="dep", timevar="rep",
                     idvar=c("subj","group"))

  ## Convert variables to factor
  fr_long <- within(fr_long, {
    group <- factor(group)
    rep <- factor(rep)
    subj <- factor(subj)
  })

  FriedmanRanks <- ezANOVA(data=fr_long, dv = .(dep), wid=.(subj), between = .(group),
                           within=.(rep), type=2, detailed=FALSE)


  #KOCH
  dif <- array(0, dim=c(n,k,k))
  rdif <- array(0, dim=c(n,k,k))
  for(j in 1:k){
    for(i in 1:k){
      dif[,i,j] <- y[,j] - y[,i]
      rdif[,i,j] <- matrix(rank(dif[,i,j]), nrow=n)
    }
  }
  som <- matrix(0,n,k)
  for(i in 1:k){
    som[,i] <- matrix(rowSums (rdif[,,i]))
  }
  #voor adjustment
  #som
  # deel1 = (n+1)/2
  # deel2 = (k-1)*(n+1)
  #na adjustment
  # q <- (som-deel1)/deel2
  q <- som
  q_full <- cbind(subj, group, q)

  l <- ncol(q_full)
  q_long <- reshape(q_full, direction="long",
                    varying=list(names(q_full)[3:l]),
                    v.names="dep", timevar="rep",
                    idvar=c("subj","group"))
  ## Convert variables to factor
  q_long <- within(q_long, {
    group <- factor(group)
    rep <- factor(rep)
    subj <- factor(subj)
  })

  KochRanks <- ezANOVA(data=q_long, dv = .(dep), wid=.(subj), between = .(group),
                       within=.(rep), type=2, detailed=FALSE)
  ANOVAs <- cat(format("\nIn the AlignCheck ANOVA, the p-values for 'group' and for 'rep' should be = 1.\n(Ignore ANOVA 'NULL')\n", sep="", width=80, justify="centre"))
  my_list <- list(ANOVAs, Parametric, AlignCheck, RegularRanks, FriedmanRanks, KochRanks)
  names(my_list) <- c("ANOVAs", "Parametric", "AlignCheck", "RegularRanks", "FriedmanRanks", "KochRanks")


  # my_list
  return(my_list)

}



