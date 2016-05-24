#' @name choice_matrix
#' @title Choice Matrix Data
#' @description This data set contains the subject's choices and point values for the difficult questions.
#' @docType data
#' @usage data(choice_matrix)
#' @format A data frame with 3780 observations on the following 5 variables.
#' \describe{
#'   \item{\code{subject_id}}{Research Participant Subject ID. There are 102 IDs and each ID has 15 observations.}
#'   \item{\code{hard_q_id}}{The item ID of the hard question assigned to the student [16-30]}
#'   \item{\code{easy_q_id}}{The item ID of the easy question assigned to the student [1-15]}
#'   \item{\code{choose_hard_q}}{Selected either: Difficult Question (1) or Easy Question (0)}
#'   \item{\code{high_value}}{Range of values associated with Difficult Question that span from 12 to 16, repeated three times per subject}
#'   \item{\code{low_value}}{Range of values associated with Easy Question that span from 4 to 6, repeated five times per subject}
#'   \item{\code{is_correct_choice}}{Did the user select an item that was answered correctly?}
#' }
#' @source Choice38 Experiment at UIUC during Spring 2014 - Fall 2014
#' @author Steven Culpepper and James Balamuta
NULL