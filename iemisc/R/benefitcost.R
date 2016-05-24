#' Benefit-Cost Ratio (Engineering Economics)
#'
#' Compute the benefit-cost ratio between two alternatives
#'
#' Benefit is expressed as
#'
#' 	\deqn{Benefit = AB\left[\frac{\left(1 + i\right)^n - 1}{i\left(1 + i\right)^n}\right]}
#'
#' \describe{
#'	\item{\emph{Benefit}}{the present equivalent benefit}
#'	\item{\emph{AB}}{the annual benefit}
#'	\item{\emph{i}}{the "effective interest rate" per year}
#'	\item{\emph{n}}{the number of years}
#' }
#'
#'
#'
#' Cost is expressed as
#'
#' 	\deqn{Cost = PC + OM\left[\frac{\left(1 + i\right)^n - 1}{i\left(1 + i\right)^n}\right] - S\left[\frac{1}{\left(1 + i\right)^n}\right]}
#'
#' \describe{
#'	\item{\emph{Cost}}{the present equivalent cost}
#'	\item{\emph{PC}}{the present or initial cost}
#'	\item{\emph{OM}}{the annual operations & maintenance cost}
#'	\item{\emph{S}}{the salvage value}
#'	\item{\emph{i}}{the "effective interest rate" per year}
#'	\item{\emph{n}}{the number of years}
#' }
#'
#'
#'
#' Benefit-Cost ratio is expressed as
#'
#' 	\deqn{BC = \frac{B_2 - B_1}{C_2 - C_1} \geq 1}
#'
#' \describe{
#'	\item{\emph{BC}}{the present equivalent cost}
#'	\item{\emph{\eqn{B_1}}}{the benefit for alternative 1}
#'	\item{\emph{\eqn{B_2}}}{the benefit for alternative 2}
#'	\item{\emph{\eqn{C_1}}}{the cost for alternative 1}
#'	\item{\emph{\eqn{C_2}}}{the cost for alternative 2}
#' }
#'
#'
#'
#' @param ic1 numeric vector that contains the initial cost for option 1
#' @param n1 numeric vector that contains the useful life (years) for option 1
#' @param ac1 numeric vector that contains the annual cost [operations &
#'    maintenance (O&M)] for option 1
#' @param i1 numeric vector that contains the effective interest rate per
#'    period as a percent for option 1
#' @param ab1 numeric vector that contains the annual benefits for option 1
#' @param salvage1 numeric vector that contains the salvage value for option 1
#' @param option1 character vector that contains the name of option for
#'    option 1
#' @param ic2 numeric vector that contains the initial cost for option 2
#' @param n2 numeric vector that contains the useful life (years) for option 2
#' @param ac2 numeric vector that contains the annual cost [operations &
#'    maintenance (O&M)] for option 2
#' @param i2 numeric vector that contains the effective interest rate per
#'    period as a percent for option 2
#' @param ab2 numeric vector that contains the annual benefits for option 2
#' @param salvage2 numeric vector that contains the salvage value for option 2
#' @param option2 character vector that contains the name of option for
#'    option 2
#' @param table character vector that contains the table output table format
#'    (ptable, rtable, or both)
#'
#'
#' @return data.frame with character vectors with the monetary values having
#' thousands separator in a pretty table (ptable) & message with the best
#' option, data.frame with numeric vectors without the thousands separator in
#' regular table (rtable) & a message with the best option, or both options
#' combined in a list
#'
#'
#' @references
#' \enumerate{
#' \item Michael R. Lindeburg, PE, \emph{EIT Review Manual}, Belmont, California: Professional Publications, Inc., 1996, page 14-2, 14-4.
#' \item William G. Sullivan, Elin M. Wicks, and C. Patrick Koelling, \emph{Engineering Economy}, Fourteenth Edition, Upper Saddle River, New Jersey: Pearson/Prentice Hall, 2009, page 133, 142, 442-443, 452-453.
#' }
#'
#' @encoding UTF-8
#'
#'
#' @examples
#' library(iemisc)
#' # Example from Lindeburg Reference text (page 14-4)
#' benefitcost(ic1 = 300000, n1 = 10, ac1 = 45000, ab1 = 150000, i1 = 10,
#' salvage1 = 0, ic2 = 400000, n2 = 10, ac2 = 35000, ab2 = 200000, i2 = 10,
#' salvage2 = 10000, option1 = "A", option2 = "B", table = "rtable")
#'
#'
#' # This is useful for saving the results as the named data.frame rtable
#' rtable <- benefitcost(ic1 = 300000, n1 = 10, ac1 = 45000, ab1 = 150000,
#' i1 = 10, salvage1 = 0, ic2 = 400000, n2 = 10, ac2 = 35000, ab2 = 200000,
#' i2 = 10, salvage2 = 10000, option1 = "A", option2 = "B", table = "rtable")
#'
#'
#' # This is useful for saving the results as the named data.frame ptable
#' ptable <- benefitcost(ic1 = 300000, n1 = 10, ac1 = 45000, ab1 = 150000,
#' i1 = 10, salvage1 = 0, ic2 = 400000, n2 = 10, ac2 = 35000, ab2 = 200000,
#' i2 = 10, salvage2 = 10000, option1 = "A", option2 = "B", table = "ptable")
#'
#'
#' # This is useful for saving the results as the named list of 2 data.frames
#' # called both
#' both <- benefitcost(ic1 = 300000, n1 = 10, ac1 = 45000, ab1 = 150000,
#' i1 = 10, salvage1 = 0, ic2 = 400000, n2 = 10, ac2 = 35000, ab2 = 200000,
#' i2 = 10, salvage2 = 10000, option1 = "A", option2 = "B", table = "both")
#'
#'
#'
#'
#' # Example 10-8 from the Sullivan Reference text (page 452-453)
#' project <- benefitcost(ic1 = 750000, n1 = 35, ac1 = 120000, ab1 = 245000,
#' i1 = 9, salvage1 = 0, ic2 = 625000, n2 = 25, ac2 = 110000, ab2 = 230000,
#' i2 = 9, salvage2 = 0, option1 = "Project I", option2 = "Project II",
#' table = "rtable")
#'
#'
#'
#'
#' @export
benefitcost <- function (ic1, n1, ac1, ab1, i1, salvage1, ic2, n2, ac2, ab2, i2, salvage2, option1, option2, table = c("ptable", "rtable", "both")) {

# option 1
# cost 1
i1 <- i1 / 100

cost1 <- ic1 + (ac1 * (((1 + i1) ^ n1 - 1) / (i1 * ((1 + i1) ^ n1)))) - (salvage1 * ((1 / (1 + i1)) ^ n1))

# benefit 1
benefit1 <- ab1 * (((1 + i1) ^ n1 - 1) / (i1 * ((1 + i1) ^ n1)))

b_c1 <- benefit1 / cost1


# option 2
# cost 2
i2 <- i2 / 100

cost2 <- ic2 + (ac2 * (((1 + i2) ^ n2 - 1) / (i2 * ((1 + i2) ^ n2)))) - (salvage2 * ((1 / (1 + i2)) ^ n2))

# benefit 2
benefit2 <- ab2 * (((1 + i2) ^ n2 - 1) / (i2 * ((1 + i2) ^ n2)))

b_c2 <- benefit2 / cost2


# Benefit-Cost ratio of Option 2 to Option 1
b_c_rank <- (benefit2 - benefit1) / (cost2 - cost1)
b_c_rank_choose <- ifelse(b_c_rank >= 1, paste0("choose", " ", option2, "."), paste0("choose", " ", option1, "."))


if (table == "ptable") {

ptable <- data.frame(c(NA_character_, "Benefit", "Cost", "Benefit-Cost Ratio"), c(option1, formatC(benefit1, big.mark = ",", format = "f", digits = 2), formatC(cost1, big.mark = ",", format = "f", digits = 2), round(b_c1, digits = 2)), c(option2, formatC(benefit2, big.mark = ",", format = "f", digits = 2), formatC(cost2, big.mark = ",", format = "f", digits = 2), round(b_c2, digits = 2)), stringsAsFactors = FALSE)
col.names <- as.character(ptable[1, ])
colnames(ptable) <- col.names
ptable <- ptable[-1, ]

cat("\n", paste("The Benefit-Cost ratio of", option2, "to", option1, "is", round(b_c_rank, digits = 2), "thus", b_c_rank_choose), "\n\n")

# code block below modified from data.frame function
attr(ptable, "col.names") <- col.names
attr(ptable, "class") <- "data.frame"
ptable


} else if (table == "rtable") {

rtable <- data.frame(c("Benefit", "Cost", "Benefit-Cost Ratio"), c(round(benefit1, digits = 2), round(cost1, digits = 2), round(b_c1, digits = 2)), c(round(benefit2, digits = 2), round(cost2, digits = 2), round(b_c2, digits = 2)), stringsAsFactors = FALSE)
col.names <- c(NA_character_, option1, option2)
colnames(rtable) <- col.names

cat("\n", paste("The Benefit-Cost ratio of", option2, "to", option1, "is", round(b_c_rank, digits = 2), "thus", b_c_rank_choose), "\n\n")

# code block below modified from data.frame function
attr(rtable, "col.names") <- col.names
attr(rtable, "class") <- "data.frame"
rtable


} else if (table == "both") {

# ptable
ptable <- data.frame(c(NA_character_, "Benefit", "Cost", "Benefit-Cost Ratio"), c(option1, formatC(benefit1, big.mark = ",", format = "f", digits = 2), formatC(cost1, big.mark = ",", format = "f", digits = 2), round(b_c1, digits = 2)), c(option2, formatC(benefit2, big.mark = ",", format = "f", digits = 2), formatC(cost2, big.mark = ",", format = "f", digits = 2), round(b_c2, digits = 2)), stringsAsFactors = FALSE)
col.names <- as.character(ptable[1, ])
colnames(ptable) <- col.names
ptable <- ptable[-1, ]

cat("\n", paste("The Benefit-Cost ratio of", option2, "to", option1, "is", round(b_c_rank, digits = 2), "thus", b_c_rank_choose), "\n\n")

# code block below modified from data.frame function
attr(ptable, "col.names") <- col.names
attr(ptable, "class") <- "data.frame"


# rtable
rtable <- data.frame(c("Benefit", "Cost", "Benefit-Cost Ratio"), c(round(benefit1, digits = 2), round(cost1, digits = 2), round(b_c1, digits = 2)), c(round(benefit2, digits = 2), round(cost2, digits = 2), round(b_c2, digits = 2)), stringsAsFactors = FALSE)
col.names <- c(NA_character_, option1, option2)
colnames(rtable) <- col.names

# code block below modified from data.frame function
attr(rtable, "col.names") <- col.names
attr(rtable, "class") <- "data.frame"

return(list(rtable, ptable))
}
}
