#' Deviation from the mode (DM)
#'
#' Computes the deviation from the mode for a vector of frequencies of categories.
#'
#' According to Wilcox (1973, p. 327), 'the measure can be thought of as an index of deviation from the modal frequency, analogous to the variance as a measure of deviation from the mean'. The formula for the DM is:
#' \deqn{1 - \frac{\sum_{i = 1}^k (f_m - f_i)}{N(K-1)}}
#'
#'
#' @param x a vector of frequencies
#' @param na.rm if TRUE, missing values are removed. If FALSE, NA is returned if there is any NA value.
#' @return The value of the DM statistics, which is normalised (varies between 0 and 1).
#' @references Wilcox, Allen R. 'Indices of Qualitative Variation and Political Measurement.' \emph{The Western Political Quarterly} 26, no. 2 (1 June 1973): 325-43. doi:10.2307/446831.
#' @examples
#' x <- rmultinom(1, 100, rep_len(0.25, 4))
#' x <- as.vector(t(x))
#' names(x) <- c("a", "b", "c", "d")
#' DM(x)
#'
#' df <- rmultinom(10, 100, rep_len(0.25, 4))
#' df <- as.data.frame(t(df))
#' names(df) <- c("a", "b", "c", "d")
#' apply(df, 1, DM)
#' @importFrom stats na.omit
#' @export
DM <- function(x, na.rm = TRUE) {
  if (na.rm)
    x <- na.omit(x)
  mode <- max(x)
  1 - (sum(mode - x) / (sum(x) * (length(x) - 1)))
}

#' Average Deviation Analog (ADA)
#'
#' Computes the average deviation analog (ADA) for a vector of frequencies of
#' categories.
#'
#' According to Wilcox (1973, p. 328), the ADA is 'an analog of the average or
#' mean deviation'. The formula for the ADA is:
#' \deqn{1 - \frac{\sum_{i=1}^k \left| f_i - \frac{N}{K}\right|}{2 \frac{N}{K}(K-1)}}
#'
#' @param x a vector of frequencies
#' @param na.rm if TRUE, missing values are removed. If FALSE, NA is returned if there is any NA value.
#' @return The value of the ADA statistics, which is normalised (varies between 0
#'   and 1).
#' @references Wilcox, Allen R. 'Indices of Qualitative Variation and Political
#'   Measurement.' \emph{The Western Political Quarterly} 26, no. 2 (1 June
#'   1973): 325-43. doi:10.2307/446831.
#' @examples
#' x <- rmultinom(1, 100, rep_len(0.25, 4))
#' x <- as.vector(t(x))
#' ADA(x)
#'
#' df <- rmultinom(10, 100, rep_len(0.25, 4))
#' df <- as.data.frame(t(df))
#' apply(df, 1, ADA)
#' @importFrom stats na.omit
#' @export
ADA <- function(x, na.rm = TRUE) {
  if (na.rm)
    x <- na.omit(x)
  nk <- mean(x)
  1 - (sum(abs(x - nk)) / (2 * nk * (length(x) - 1)))
}

#' Mean Difference Analog (MDA)
#'
#' Computes the mean difference analog (MDA) for a vector of frequencies of
#' categories.
#'
#' According to Wilcox (1973, p. 328), the MDA is 'an analog of the mean
#' difference, a measure of variation that is discussed and used much less
#' frequently than the average deviation or the standard deviation. It is
#' defined as "the average of the differences of all the possible pairs of
#' variate-values, taken regardless of sign"'. The formula for the MDA is:
#' \deqn{1 - \frac{\sum_{i=1}^{k-1} \sum_{j=i+1}^k |f_i - f_j|}{N(K-1)}}
#'
#' @param x a vector of frequencies
#' @param na.rm if TRUE, missing values are removed. If FALSE, NA is returned if there is any NA value.
#' @return The value of the MDA statistics, which is normalised (varies between
#'   0 and 1).
#' @references Wilcox, Allen R. 'Indices of Qualitative Variation and Political
#'   Measurement.' \emph{The Western Political Quarterly} 26, no. 2 (1 June
#'   1973): 325-43. doi:10.2307/446831.
#' @examples
#' x <- rmultinom(1, 100, rep_len(0.25, 4))
#' x <- as.vector(t(x))
#' MDA(x)
#'
#' df <- rmultinom(10, 100, rep_len(0.25, 4))
#' df <- as.data.frame(t(df))
#' apply(df, 1, MDA)
#' @importFrom stats na.omit
#' @export
MDA <- function(x, na.rm = TRUE) {
  if (na.rm)
    x <- na.omit(x)
  K <- length(x)
  N <- sum(x)
  s <- 0
  for (i in 1:(K-1)) {
    for (j in (i+1):K) {
      s <- s + abs(x[i]-x[j])
    }
  }
  1 - s/(N*(K-1))
}

#' Variance Analog (VA)
#'
#' Computes the variance analog (VA) for a vector of frequencies of categories.
#'
#' According to Wilcox (1973, p. 329), the VA is 'based on the variance, which
#' is defined as the arithmetic mean of the squared differences of each value
#' from the mean'. The formula for the VA is:
#' \deqn{1 - \frac{\sum_{i=1}^k \left(f_i - \frac{N}{K}\right)^2}{\frac{N^2(K-1)}{K}}}
#'
#' @param x a vector of frequencies
#' @param na.rm if TRUE, missing values are removed. If FALSE, NA is returned if there is any NA value.
#' @return The value of the VA statistics, which is normalised (varies between
#'   0 and 1).
#' @references Wilcox, Allen R. 'Indices of Qualitative Variation and Political
#'   Measurement.' \emph{The Western Political Quarterly} 26, no. 2 (1 June
#'   1973): 325-43. doi:10.2307/446831.
#' @examples
#' x <- rmultinom(1, 100, rep_len(0.25, 4))
#' x <- as.vector(t(x))
#' VA(x)
#'
#' df <- rmultinom(10, 100, rep_len(0.25, 4))
#' df <- as.data.frame(t(df))
#' apply(df, 1, VA)
#' @importFrom stats na.omit
#' @export
VA <- function(x, na.rm = TRUE) {
  if (na.rm)
    x <- na.omit(x)
  N <- sum(x)
  K <- length(x)
  1 - (sum((x - N/K)^2) / (N^2 * (K - 1) / K))
}

#' HREL
#'
#' Computes the HREL index for a vector of frequencies of categories.
#'
#' According to Wilcox (1973, p. 329), and following Senders (1958), the HREL is 'a measure originally developed by engineers for use in specifying the properties of communications channels. The rationale for HREL is presented in terms of guessing by Virginia Senders (supplementing the mode as best guess): "What we need is a measure of uncertainty, or of 'poorness of a guess,' which will be high when the number of alternative possibilities is high, and low when some ofthe possibilities are much more likely than others. One possible measure is the average number of questions we have to ask to specify the correct alternative'. The formula for the HREL is:
#' \deqn{- \frac{\sum_{i=1}^k \frac{f_i}{N} \log_2 \frac{f_i}{N}}{\log_2 K}}
#'
#' @param x a vector of frequencies
#' @param na.rm if TRUE, missing values are removed. If FALSE, NA is returned if there is any NA value.
#' @return The value of the HREL statistics, which is normalised (varies between
#'   0 and 1).
#' @references Wilcox, Allen R. 'Indices of Qualitative Variation and Political
#'   Measurement.' \emph{The Western Political Quarterly} 26, no. 2 (1 June
#'   1973): 325-43. doi:10.2307/446831.
#'   Senders, Virginia L. Measurement and Statistics. New York: Oxford University Press, 1958.
#' @examples
#' x <- rmultinom(1, 100, rep_len(0.25, 4))
#' x <- as.vector(t(x))
#' HREL(x)
#'
#' df <- rmultinom(10, 100, rep_len(0.25, 4))
#' df <- as.data.frame(t(df))
#' apply(df, 1, HREL)
#' @importFrom stats na.omit
#' @export
HREL <- function(x, na.rm = TRUE) {
  if (na.rm)
    x <- na.omit(x)
  N <- sum(x)
  K <- length(x)
  - sum(x/N * log(x/N, base = 2)) / (log(K, base = 2))
}

#' B (modified geometric mean)
#'
#' Computes the B index for a vector of frequencies of categories.
#'
#' According to Wilcox (1973, p. 330), and following Kaiser (1968), the B index
#' relies on the geometric mean, but corrects it for undesirable properties.
#' The formula for the B index is:
#' \deqn{1 - \sqrt{1 - \left(\sqrt[K]{\prod_{i=1}^k \frac{f_i K}{N}}\right)^2}}
#'
#' @param x a vector of frequencies
#' @param na.rm if TRUE, missing values are removed. If FALSE, NA is returned if there is any NA value.
#' @return The value of the B statistics, which is normalised (varies between 0
#'   and 1).
#' @references Wilcox, Allen R. 'Indices of Qualitative Variation and Political
#'   Measurement.' \emph{The Western Political Quarterly} 26, no. 2 (1 June
#'   1973): 325-43. doi:10.2307/446831.
#'   Kaiser, Henry F. 'A Measure of the Population Quality of Legislative Apportionment.' \emph{The American Political Science Review} 62, no. 1 (March 1968): 208. doi:10.2307/1953335.
#' @examples
#' x <- rmultinom(1, 100, rep_len(0.25, 4))
#' x <- as.vector(t(x))
#' B(x)
#'
#' df <- rmultinom(10, 100, rep_len(0.25, 4))
#' df <- as.data.frame(t(df))
#' apply(df, 1, B)
#' @importFrom stats na.omit
#' @export
B <- function(x, na.rm = TRUE) {
  if (na.rm)
    x <- na.omit(x)
  N <- sum(x)
  K <- length(x)
  1 - sqrt(1 - ((prod(x * K / N))^(1/K))^2)
}
