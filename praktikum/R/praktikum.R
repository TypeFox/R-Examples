#' \tabular{ll}{
#' Package: \tab praktikum\cr
#' Type: \tab Package\cr
#' Version: \tab 0.1\cr
#' Date: \tab 2014-02-14\cr
#' License: \tab GPL (>= 2)\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' Kasulikud funktsioonid kvantitatiivsete mudelite kursuse (SHPH.00.004) jaoks / Functions used in the course "Quantitative methods in behavioural sciences" (SHPH.00.004), University of Tartu
#'
#'
#' @name praktikum-package
#' @aliases praktikum
#' @docType package
#' @title Lisamoodul 'praktikum'
#' @author Kenn Konstabel \email{kenn.konstabel AT tai.ee}
#' @keywords package
#' @examples NULL
NA


#' Leia kirjeldavad statistikud terve andmetabeli kohta / Table of descriptive statistics for numeric variables in a data frame
#'
#' Funktsioonile tuleb ette anda andmetabel või selle osa. Tulemuseks on
#' tabel (data.frame) arvuliste tunnuste keskmise, standardhälbe, miinimumi ja
#' maksimumiga, samuti N-iga (st mitte-puuduvate väärtuste arv)
#' Näiteks: kirjelda(X, "Pikkus Kaal Vanus")
#' 
#'
#' @param D andmetabel
#' @param C tunnuste nimed ühe pika tekstijoruna nt "Pikkus Kaal Vanus"
#' @return Tabel (data.frame) leitud statistikute väärtustega
#' @encoding utf-8
#' @export
kirjelda <- function(D, C=NULL) {
  if(!is.null(C)) cols <- c_(C) else cols <- TRUE
  D <- D[, cols]
  safe <- function(F, x) if(is.numeric(x)) F(x, na.rm=TRUE) else NA
  combo <- function(..) data.frame(
      mean = safe(mean,..),
      sd = safe(sd,..),
      min = safe(min, ..),
      max = safe(max, ..),
      N = sum(!is.na(..))
    )
  do.call(rbind, lapply(D,combo))  
}


#' T-testide tabel koos keskmiste ja standardhälvetega / Table of t-tests with means and standard deviations
#'
#' Argumentideks on andmetabel, sõltuvate tunnuste nimekiri ja grupitunnuse nimi.
#' Grupitunnusel peab olema kindlasti 2 väärtust!
#' nt TT.test(D, "Pikkus Kaal Vanus", "Sugu")
#'
#' @param D andmetabel
#' @param V tunnuste nimed kas vektorina või ühe pika tekstijoruna nt "Pikkus Kaal Vanus"
#' @param G grupitunnus
#' @param ... ekstra parameetrid edastamiseks funktsioonile t.test nt var.equal
#' @return Tabel (data.frame) keskmiste, standardhälvete ja t-testi tulemustega
#' @encoding utf-8
#' @export
TT.test <- function(D, V, G, ...){
  GR <- D[[G]]
  if(length(unique(GR))!=2) stop("Grouping factor must have exactly 2 levels")
  cols <- c_(V)
  D <- D[, cols]
  ttt <- lapply(D, function(..) t.test(..~GR, ...))
  res <- lapply(ttt, "[", c_("statistic parameter p.value estimate"))
  res <- do.call(rbind, lapply(lapply(lapply(res, unlist), as.list), as.data.frame))
  N <- names(table(GR))
  N2 <- c("t", "df", "p", paste0("Mean[",paste(G, N, sep=":"), "]"))
  sds <- do.call(rbind, lapply(lapply(lapply(D, tapply, GR, sd), as.list), as.data.frame))
  N2 <- c(N2, paste0("SD[", paste(G, N, sep=":"), "]"))
  res <- setNames(cbind(res, sds), N2)
  res[, c(4, 6, 5, 7, 1:3)]
}

#' Korrelatsioonimaatriks koos p-väärtuste jm jamaga / Correlation matrix with p-values
#'
#' Argumendiks on andmetabel ja tunnuste nimekiri.
#' Tulemuseks on korrelatsioonimaatriks koos p-väärtustega
#' nt kor.tabel(D, "Pikkus Kaal Vanus", N=FALSE)
#'
#' @param D andmetabel
#' @param V tunnuste nimed ühe pika tekstijoruna nt "Pikkus Kaal Vanus"
#' @param N kas lisada tulp N-idega st valimi suurus iga tunnusteaari kohta, TRUE või FALSE
#' @param ... ekstra parameetrid edastamiseks funktsioonile cor.test 
#' @return Tabel  korrelatsioonide, p-väärtuste ja N-idega
#' @encoding utf-8
#' @export
kor.tabel <- function(D, V, N=TRUE, ...){
  cols <- c_(V)
  D <- D[, cols]
  asd <- as.data.frame
  rrapply2 <- function(A, FUN, ...) mapply(function(a, B) lapply(B, function(x) FUN(a, x, ...)), a = A, MoreArgs = list(B=A))
  res <- rrapply2(D, cor.test, ...)
  ps <- asd(apply(res, 1:2, function(x) x[[1]]$p.value))
  rs <- asd(apply(res, 1:2, function(x) x[[1]]$estimate))
  Ns <- asd(apply(res, 1:2, function(x) x[[1]]$parameter)+2)
  if(N) res2 <- mapply(function(r,p,N) data.frame(r=r, p=p, N=N), rs, ps, Ns, SIMPLIFY=FALSE) else res2 <- mapply(function(r,p) data.frame(r=r, p=p), rs, ps, SIMPLIFY=FALSE)
  
  res2 <- data.frame(res2)
  rownames(res2) <- rownames(res)
  res2
  }


#' Split a string into a character vector
#'
#' Splits a string: c_("a b c") is equivalent to c("a", "b", "c")
#' A wrapper around strsplit
#'
#' @param str a string
#' @param sep separator
#' @return a vector
#' @export
c_ <- function(str, sep=" "){
  strsplit(str, sep)[[1]]
}
