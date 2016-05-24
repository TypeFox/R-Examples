#' Number of words in quoted retweets
#'
#' A anonymized dataset containing the retweeting activities of 894 microblog users on 15 tweets
#'
#' @format A data frame 13410 rows and 6 columns
#' \describe{
#'   \item{tweet.id}{The id of a status posted on microblog}
#'   \item{user.id}{The id of a user on microblog}
#'   \item{fans}{The number of fans of the user, on the log scale}
#'   \item{tweets}{The number of tweets of the user, on the log scale}
#'   \item{isRetweet}{Whether the user retweets the given tweet, boolean}
#'   \item{num.words}{Number of words attached while retweeting. NA if doesn't retweet}
#' }
#' @source collected by the author of the package on microblog
"rt"