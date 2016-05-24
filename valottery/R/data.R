#' Mega Millions (Before 10/22/13)
#'
#' Historical data for the Mega Millions game. On October 22, 2013, the format changed
#' from 5/56 + 1/46 to the current 5/75 + 1/15 format. Game play: Pick five different
#' numbers from 1 through 75; then select one Mega Ball number from 1 through 15.
#' @format A data frame with 1,713 rows and 7 variables:
#' \describe{
#'   \item{date}{date of draw}
#'   \item{megaball}{megaball result}
#'   \item{N1}{1st number in order}
#'   \item{N2}{2nd number in order}
#'   \item{N3}{3rd number in order}
#'   \item{N4}{4th number in order}
#'   \item{N5}{5th number in order}
#' }
#' @source \url{https://www.valottery.com}
#'
"mega.mill.1"

#' Mega Millions (10/22/13 and beyond)
#'
#' Historical data for the Mega Millions game. On October 22, 2013, the format changed
#' from 5/56 + 1/46 to the current 5/75 + 1/15 format. Game play: Pick five different
#' numbers from 1 through 56; then select one Mega Ball number from 1 through 46.
#' @format A data frame with 194 rows and 7 variables:
#' \describe{
#'   \item{date}{date of draw}
#'   \item{megaball}{megaball result}
#'   \item{N1}{1st number in order}
#'   \item{N2}{2nd number in order}
#'   \item{N3}{3rd number in order}
#'   \item{N4}{4th number in order}
#'   \item{N5}{5th number in order}
#' }
#' @source \url{https://www.valottery.com}
#'
"mega.mill.2"

#' Power Ball
#'
#' Historical data for the Power Ball game. Game play: Pick
#' five different numbers from 1 through 59; then select one
#' Powerball number from 1 through 35.
#' @format A data frame with 582 rows and 7 variables:
#' \describe{
#'   \item{date}{date of draw}
#'   \item{powerball}{powerball result}
#'   \item{N1}{1st number in order}
#'   \item{N2}{2nd number in order}
#'   \item{N3}{3rd number in order}
#'   \item{N4}{4th number in order}
#'   \item{N5}{5th number in order}
#' }
#' @source \url{https://www.valottery.com}
#' @examples
#' ## According to game rules, the powerball is numbered 1 - 35,
#' ## but apparently there were times when it went up to 39
#' i <- power.ball$powerball > 35
#' any(i)
#' sum(i)
#' power.ball$powerball[i]
"power.ball"

#' CASH4LIFE
#'
#' Historical data for the CASH4LIFE game. Game play: Pick five different numbers from 1 through 60;
#' then select one Cash Ball number from 1 through 4.
#' @format A data frame with 34 rows and 7 variables:
#' \describe{
#'   \item{date}{date of draw}
#'   \item{cashball}{cash ball result}
#'   \item{N1}{1st number in order}
#'   \item{N2}{2nd number in order}
#'   \item{N3}{3rd number in order}
#'   \item{N4}{4th number in order}
#'   \item{N5}{5th number in order}
#' }
#' @source \url{https://www.valottery.com}
#' @examples
#' ## Check numbers drawn are uniformly distributed
#' x <- qunif(ppoints(nrow(cash.4.life)*5),1,60)
#' y <- sort(unlist(cash.4.life[,3:7]))
#' qqplot(x,y)
"cash.4.life"

#' $1,000,000 Money Ball
#'
#' Historical data for the $1,000,000 Money Ball game. Game Play: pick five numbers
#' 1 - 35, the Lottery then selects five numbered balls. If the Gold Million Dollar
#' Money Ball is drawn before all five numbers have been selected, the top prize
#' jumps to $1,000,000. Note: This game was discontinued 8/29/15.
#' @format A data frame with 100 rows and 7 variables:
#' \describe{
#'   \item{date}{date of draw}
#'   \item{moneyball}{money ball result: yes or no}
#'   \item{N1}{1st number in order}
#'   \item{N2}{2nd number in order}
#'   \item{N3}{3rd number in order}
#'   \item{N4}{4th number in order}
#'   \item{N5}{5th number in order}
#' }
#' @source \url{https://www.valottery.com}
#' @examples
#' ## probability of drawing money ball before first 5 balls
#' (1/36) + (1/35) + (1/34) + (1/33) + (1/32)
#' ## observed money ball results
#' prop.table(table(money.ball$moneyball))
#' ## simulate money ball draws before first 5 draws
#' set.seed(123)
#' mean(replicate(1000, any(sample(c(1:35,"mb"),5)=="mb")))
"money.ball"

#' Decades of Dollars
#'
#' Historical data for the Decades of Dollars game. You pick six (6) numbers,
#' the Lottery will then select six (6) numbered balls.
#' Note: This game was discontinued in favor of CASH4LIFE.
#' @format A data frame with 443 rows and 7 variables:
#' \describe{
#'   \item{date}{date of draw}
#'   \item{N1}{1st number in order}
#'   \item{N2}{2nd number in order}
#'   \item{N3}{3rd number in order}
#'   \item{N4}{4th number in order}
#'   \item{N5}{5th number in order}
#'   \item{N6}{6th number in order}
#' }
#' @source \url{https://www.valottery.com}
#'
"decades.of.dollars"

#' Cash 5 (once daily)
#'
#' Historical data for the Cash 5 once daily game. Game Play: Pick five
#' numbers from 1 through 34. Note: On April 11, 1999, Cash 5 switched
#' to twice daily drawings.
#' @format A data frame with 1187 rows and 6 variables:
#' \describe{
#'   \item{date}{date of draw}
#'   \item{N1}{1st number in order}
#'   \item{N2}{2nd number in order}
#'   \item{N3}{3rd number in order}
#'   \item{N4}{4th number in order}
#'   \item{N5}{5th number in order}
#' }
#' @source \url{https://www.valottery.com}
#'
"cash.5.1xday"

#' Cash 5 (twice daily)
#'
#' Historical data for the Cash 5 twice daily game. Game Play: Pick five
#' numbers from 1 through 34. Note: On April 11, 1999, Cash 5 switched
#' to twice daily drawings.
#' @format A data frame with 11,164 rows and 7 variables:
#' \describe{
#'   \item{date}{date of draw}
#'   \item{time}{time of drawing: day or night}
#'   \item{N1}{1st number in order}
#'   \item{N2}{2nd number in order}
#'   \item{N3}{3rd number in order}
#'   \item{N4}{4th number in order}
#'   \item{N5}{5th number in order}
#' }
#' @source \url{https://www.valottery.com}
#' @examples
#' max.days <- apply(subset(cash.5.2xday,time=="day",-(1:2)),1,max)
#' max.nights <- apply(subset(cash.5.2xday,time=="night",-(1:2)),1,max)
#' op <- par(mfrow=c(1,2))
#' hist(max.days)
#' hist(max.nights)
#' par(op)
"cash.5.2xday"

#' Pick 4 (once daily)
#'
#' Historical data for the Pick 4 once daily game. Game play: Pick
#' a four digit number from 0000 through 9999. Note: On January 30, 1995,
#' Pick 4 switched to twice daily drawings.
#' @format A data frame with 1,041 rows and 5 variables:
#' \describe{
#'   \item{date}{date of draw}
#'   \item{N1}{1st digit}
#'   \item{N2}{2nd digit}
#'   \item{N3}{3rd digit}
#'   \item{N4}{4th digit}
#' }
#' @source \url{https://www.valottery.com}
#' @examples
#' ## Any Pick 4 happen more than once?
#' results <- apply(pick.4.1xday[,-1],1,function(x)paste(x,collapse = ""))
#' any(table(results) > 1)
#' ## Which numbers?
#' i <- which(table(results) > 1,useNames = FALSE)
#' sort(table(results)[i],decreasing = TRUE)
"pick.4.1xday"

#' Pick 4 (twice daily)
#'
#' Historical data for the Pick 4 twice daily game. Game play: Pick
#' a four digit number from 0000 through 9999. Note: On January 30, 1995,
#' Pick 4 switched to twice daily drawings.
#' @format A data frame with 13,788 rows and 6 variables:
#' \describe{
#'   \item{date}{date of draw}
#'   \item{time}{time of drawing: day or night}
#'   \item{N1}{1st digit}
#'   \item{N2}{2nd digit}
#'   \item{N3}{3rd digit}
#'   \item{N4}{4th digit}
#' }
#' @source \url{https://www.valottery.com}
#'
"pick.4.2xday"

#' Pick 3 (once daily)
#'
#' Historical data for the Pick 3 once daily game. Game Play:
#' Pick a three digit number from 000 through 999. Note: On January 30, 1995,
#' Pick 3 switched to twice daily drawings.
#' @format A data frame with 1,777 rows and 4 variables:
#' \describe{
#'   \item{date}{date of draw}
#'   \item{N1}{1st digit}
#'   \item{N2}{2nd digit}
#'   \item{N3}{3rd digit}
#' }
#' @source \url{https://www.valottery.com}
#' @examples
#' lapply(pick.3.1xday[,-1],function(x)round(prop.table(table(x)),2))
"pick.3.1xday"

#' Pick 3 (twice daily)
#'
#' Historical data for the Pick 3 twice daily game. Game Play:
#' Pick a three digit number from 000 through 999. Note: On January 30, 1995,
#' Pick 3 switched to twice daily drawings.
#' @format A data frame with 13,790 rows and 5 variables:
#' \describe{
#'   \item{date}{date of draw}
#'   \item{time}{time of drawing: day, night1 or night2 (Note: two nightly drawings were held
#'   on 10/30/08 and 11/09/08)}
#'   \item{N1}{1st digit}
#'   \item{N2}{2nd digit}
#'   \item{N3}{3rd digit}
#' }
#' @source \url{https://www.valottery.com}
#'
"pick.3.2xday"
