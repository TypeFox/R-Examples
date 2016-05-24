#' Compares FSMs
#'
#' \code{compare_fsm} uses a specified distance measure to compare FSMs.
#'
#' Compares a user-defined FSM to a decoded estimated FSM. If you have have FSMs
#' that may have values in the matrices that arent all simple integers, you can
#' use the distance metric that is most appropriate. Euclidean does
#' sqrt(sum((x_i - y_i)^2)) - the L2 norm. Manhattan takes abs diff between
#' them - the L1 norm. Binary treats non-zero elements as "on" and zero
#' elements as "off" and distance is the proportion of bits in which only one is
#' on amongst those in which at least one is on.
#'
#' @param users Numeric vector or numeric matrix with a predefined FSM
#' @param gas Numeric vector or numeric matrix with an evolved FSM
#' @param comparison Character string of length one with either "manhattan",
#'   "euclidean", or "binary".
#'
#' @return Numeric vector of length one for the distance between the two
#'   supplied FSMs, calculated according to the comparison argument.
#'
#' @export

compare_fsm <- function(users, gas, comparison = "manhattan"){
  gas <- as.vector(gas)
  users <- as.vector(users)
  both <- rbind(gas[-c(2, 3, 6, 7)], users[-c(2, 3, 6, 7)]) # dont compare the non-identifiable parameters
  # these are non-identifiable when there are 4 inputs and the second row, 1st and 3rd cols of
  # the state matrix cannot be reached by a deterministic strategy.

  as.numeric(stats::dist(x=both, method = comparison))
}
################################################################################


# ################################################################################
# ## compare to tft  #############################################################
# ################################################################################
# # states=2; inputs=4; actions=2
# # use this matrix to visualize the state matrix
#  tft.state <- matrix(c(1,1,2,2, 1, 1, 2, 2),
#                   nrow = 2, ncol = 4, byrow = FALSE)
# tft.action <- c(1,2)
# # now include the action.vec to the string for comparison to empirical results
# tft <- c(1,2, 1,1,2,2, 1, 1, 2, 2)
# #
# # complexity(evolved.models$decoded[[2]]) #28.34225
# # complexity(tft) # 28.16787
# #
# # compare.fsm(tft, evolved.models$decoded[[2]]) # 5 (makes sense, there are 5 diff entries)
# ################################################################################


################################################################################
## compare to many  ############################################################
################################################################################

# create a list of user-defined GAs and then loop through that list, comparing
# the GA-derived FSM to each user-defined GA and seeing which it is closest to.

################################################################################
## to calculate the Kolmogorov complexity of the state matrix ##################
################################################################################
# see \citep{soler-toscano_calculating_2014} for paper on the mathematics of this process:
# Calculating Kolmogorov Complexity from the Output Frequency Distributions of Small Turing Machines

# require(acss)
#
# complexity <- function(fsm, possible.symbols = 2){
#         # takes a decoded state matrix of a fsm
#         # acss cannot handle large strings > 12 (which is why we
#         # need to use a decoded string, which is smaller).
#
#         fsm <- unlist(fsm)
#         # now convert to a char vec of length 1 by collapsing the vector
#         fsm <- paste(fsm, collapse="")
#
#         # returns matrix in which the row name is the string entered and the
#         # columns are (1) algorithmic complexity and (2) algorithmic posterior
#         # probability that a random process created the string.
#         acss(string = fsm, alphabet= possible.symbols)
# }
################################################################################






