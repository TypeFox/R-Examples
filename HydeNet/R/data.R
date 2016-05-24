#' Blackjack Dealer Outcome Probabilities
#'
#' A dataset containing the conditional probability of various dealer
#' outcomes given the "upcard". (The dealer and player each get two
#' cards; only one of the dealer's cards is shown, and this is called
#' the "upcard")
#'
#' @format A data frame with 70 rows and 3 variables:
#' \describe{
#'   \item{dealerUpcard}{dealer upcard}
#'   \item{dealerOutcome}{outcome of dealer's hand, under the rule that
#'   cards are drawn until the dealer's hand total is at least 17}
#'   \item{probability}{conditional probability of dealerOutcome given
#'   dealerUpcard}
#'   ...
#' }
#' @source \url{https://www.blackjackinfo.com/dealer-outcome-probabilities}
"BJDealer"

#' Black Jack Hybrid Decision Network
#'
#' An object of class \code{HydeNetwork} establishing a graphical model for a game
#' of Black Jack.
#'
#' @format A \code{HydeNetwork} object constructed using the code shown in 
#' the example.  The network has seven random nodes, three ten deterministic
#' nodes, three decision nodes, and one utility node.  This is (almost) the same 
#' network used in `vignette("DecisionNetworks", package="HydeNet")`.
#' 
#' @examples
#' \dontrun{
#' BlackJack <- 
#'        HydeNetwork(~ initialAces | card1*card2
#'                    + initialPoints | card1*card2
#'                    + highUpcard | dealerUpcard
#'                    + hit1 | initialPoints*highUpcard
#'                    + acesAfterCard3 | initialAces*card3
#'                    + pointsAfterCard3 | card1*card2*card3*acesAfterCard3
#'                    + hit2 | pointsAfterCard3*highUpcard
#'                    + acesAfterCard4 | acesAfterCard3*card4
#'                    + pointsAfterCard4 | card1*card2*card3*card4*acesAfterCard4
#'                    + hit3 | pointsAfterCard4*highUpcard
#'                    + acesAfterCard5 | acesAfterCard4*card5
#'                    + pointsAfterCard5 | card1*card2*card3*card4*card5*acesAfterCard5
#'                    + playerFinalPoints | initialPoints*hit1*pointsAfterCard3
#'                    *hit2*pointsAfterCard4*hit3*pointsAfterCard5
#'                    + dealerFinalPoints | dealerUpcard
#'                    + payoff | playerFinalPoints*dealerFinalPoints)
#' cardProbs  <- c(rep(1/13,8), 4/13, 1/13)  # probs. for 2, 3, ..., 9, (10-K), A
#' 
#' BlackJack <- setNode(BlackJack, card1, nodeType="dcat",  
#'                      pi=vectorProbs(p=cardProbs, card1))
#' BlackJack <- setNode(BlackJack, card2, nodeType="dcat",  
#'                      pi=vectorProbs(p=cardProbs, card2))
#' BlackJack <- setNode(BlackJack, card3, nodeType="dcat",  
#'                      pi=vectorProbs(p=cardProbs, card3))
#' BlackJack <- setNode(BlackJack, card4, nodeType="dcat",  
#'                      pi=vectorProbs(p=cardProbs, card4))
#' BlackJack <- setNode(BlackJack, card5, nodeType="dcat",  
#'                      pi=vectorProbs(p=cardProbs, card5))
#' 
#' BlackJack <- setNode(BlackJack, dealerUpcard, nodeType="dcat",
#'                pi=vectorProbs(p=cardProbs, dealerUpcard))
#' 
#' #Note: node dealerFinalPoints will be defined below, following some discussion 
#' #      about its conditional probability distribution.
#' 
#' #####################################
#' # Deterministic Nodes
#' #####################################
#' BlackJack <- setNode(BlackJack, highUpcard,     
#'                "determ", define=fromFormula(),
#'                nodeFormula = highUpcard ~ ifelse(dealerUpcard > 8, 1, 0))
#' 
#' BlackJack <- setNode(BlackJack, initialAces,    
#'                "determ", define=fromFormula(),
#'                nodeFormula = initialAces ~ ifelse(card1==10,1,0) + 
#'                               ifelse(card2==10,1,0))
#' 
#' BlackJack <- setNode(BlackJack, acesAfterCard3, 
#'                "determ", define=fromFormula(),
#'                nodeFormula = acesAfterCard3 ~ initialAces + ifelse(card3==10,1,0))
#' 
#' BlackJack <- setNode(BlackJack, acesAfterCard4, 
#'                "determ", define=fromFormula(),
#'                nodeFormula = acesAfterCard4 ~ acesAfterCard3 + ifelse(card4==10,1,0))
#' 
#' BlackJack <- setNode(BlackJack, acesAfterCard5, 
#'                "determ", define=fromFormula(),
#'                nodeFormula = acesAfterCard5 ~ acesAfterCard4 + ifelse(card5==10,1,0))
#' 
#' BlackJack <- setNode(BlackJack, initialPoints, 
#'                "determ", define=fromFormula(),
#'                nodeFormula = initialPoints ~ card1+card2+2)
#' 
#' BlackJack <- setNode(BlackJack, pointsAfterCard3, "determ", define=fromFormula(),
#'                nodeFormula = pointsAfterCard3 ~
#'                  ifelse(acesAfterCard3 == 3,
#'                         13,
#'                         ifelse(acesAfterCard3 == 2,
#'                                card1 + card2 + card3 + 3 - 10,
#'                                ifelse(acesAfterCard3 == 1,
#'                                       ifelse(card1 + card2 + card3 + 3 > 22,
#'                                              card1 + card2 + card3 + 3 - 10,
#'                                              card1 + card2 + card3 + 3),
#'                                       card1 + card2 + card3 + 3
#'                                )
#'                         )
#'                  )
#' )
#' 
#' BlackJack <- setNode(BlackJack, pointsAfterCard4, "determ", define=fromFormula(),
#'                nodeFormula = pointsAfterCard4 ~
#'                  ifelse(acesAfterCard4 == 4,
#'                         14,
#'                         ifelse(acesAfterCard4 == 3,
#'                                ifelse(card1 + card2 + card3 + card4 + 4 > 38,
#'                                       card1 + card2 + card3 + card4 + 4 - 30,
#'                                       card1 + card2 + card3 + card4 + 4 - 20
#'                                ),
#'                                ifelse(acesAfterCard4 > 0,
#'                                       ifelse(card1 + card2 + card3 + card4 + 4 > 22,
#'                                              card1 + card2 + card3 + card4 + 4 - 10,
#'                                              card1 + card2 + card3 + card4 + 4
#'                                       ),
#'                                       card1 + card2 + card3 + card4 + 4
#'                                )
#'                         )
#'                  )
#' )
#' 
#' BlackJack <- 
#'   setNode(BlackJack, pointsAfterCard5, "determ", define=fromFormula(),
#'           nodeFormula = pointsAfterCard5 ~ 
#'             ifelse(acesAfterCard5 == 5,
#'               15,
#'               ifelse(acesAfterCard5 == 4,
#'                 ifelse(card1 + card2 + card3 + card4 + card5 + 5 > 51,
#'                   card1 + card2 + card3 + card4 + card5 + 5 - 40,
#'                   card1 + card2 + card3 + card4 + card5 + 5 - 30
#'                 ),
#'                 ifelse(acesAfterCard5 == 3,
#'                   ifelse(card1 + card2 + card3 + card4 + card5 + 5 > 51,
#'                     card1 + card2 + card3 + card4 + card5 + 5 - 30,
#'                     card1 + card2 + card3 + card4 + card5 + 5 - 20
#'                   ),
#'                   ifelse(acesAfterCard5 == 2,
#'                     ifelse(card1 + card2 + card3 + card4 + card5 + 5 > 31,
#'                       card1 + card2 + card3 + card4 + card5 + 5 - 20,
#'                       card1 + card2 + card3 + card4 + card5 + 5 - 10
#'                     ),
#'                     ifelse(acesAfterCard5 > 0,
#'                       ifelse(card1 + card2 + card3 + card4 + card5 + 5 > 22,
#'                         card1 + card2 + card3 + card4 + card5 + 5 - 10,
#'                         card1 + card2 + card3 + card4 + card5 + 5
#'                       ),
#'                       card1 + card2 + card3 + card4 + card5 + 5
#'                     )
#'                   )
#'                 )
#'               )
#'             )
#' )
#' 
#' BlackJack <- setNode(BlackJack, playerFinalPoints, "determ", define=fromFormula(),
#'                nodeFormula = playerFinalPoints ~ 
#'                  ifelse(hit1 == 0,
#'                         initialPoints,
#'                         ifelse(hit2 == 0,
#'                                pointsAfterCard3,
#'                                ifelse(hit3 == 0, pointsAfterCard4, pointsAfterCard5)
#'                         )
#'                  )
#' )
#' 
#' BlackJack <- setDecisionNodes(BlackJack, hit1, hit2, hit3)
#' BlackJack <- setUtilityNodes(BlackJack, payoff)
#' }
#' 
"BlackJack"

#' Black Jack Network Training Dataset
#'
#' These are simulated data on 1,000 Black Jack hands.
#'
#' @format A data frame with 10000 rows and 7 variables:
#' \describe{
#'   \item{dealerUpcard}{The card in the dealer's hand visible to all players}
#'   \item{card1}{Value of the first card}
#'   \item{card2}{Value of the second card}
#'   \item{initialPoints}{Total points with the two cards}
#'   \item{hit1}{Binary variable indicating if a hit was taken}
#'   \item{card3}{Value of the third card}
#'   \item{pointsAfterCard3}{Total points with three cards}
#'   \item{hit2}{Binary variable indicating if a hit was taken}
#'   \item{card4}{Value of the fourth card}
#'   \item{pointsAfterCard4}{Total points with four cards}
#'   \item{hit3}{Binary variable indicating if a hit was taken}
#'   \item{card5}{Value of the fifth card}
#'   \item{pointsAfterCard5}{Total points with five cards}
#' }
#' @source 
#' Bicycle Cards, "Blackjack," 
#' Retrieved from http://www.bicyclecards.com/card-games/rule/blackjack

"BlackJackTrain"

#' JAGS Probability Distributions.
#'
#' A dataset listing the JAGS probability distributions and their parameters
#'
#' @format A data frame with 30 rows and 7 variables:
#' \describe{
#'   \item{DistName}{Distribution Name}
#'   \item{FnName}{Function Name}
#'   \item{xLow}{Minimum value for x, the random variable}
#'   \item{xHigh}{Maximum value for x, the random variable}
#'   \item{Parameters}{Names of the parameters}
#'   \item{paramLimit}{Limits on the parameter}
#'   \item{paramLogic}{The text of a logical check used in \code{setNode} to 
#'     ensure stated parameters are valid.}
#' }
#' @source \url{http://people.math.aau.dk/~kkb/Undervisning/Bayes14/sorenh/docs/jags_user_manual.pdf}
"jagsDists"

#' JAGS Functions Compatible with R.
#'
#' A dataset listing the JAGS functions and their R equivalents.
#'
#' @format A data frame with 30 rows and 3 variables:
#' \describe{
#'   \item{jags_function}{JAGS function name}
#'   \item{r_function}{R function Name}
#'   \item{r_package}{R package where the function is found.}
#' }
#' @source \url{http://people.math.aau.dk/~kkb/Undervisning/Bayes14/sorenh/docs/jags_user_manual.pdf}
"jagsFunctions"

#' Pulmonary Embolism Dataset
#'
#' These are simulated data on 10,000 cases with suspected pulmonary embolism at a hospital.
#'
#' @format A data frame with 10000 rows and 7 variables:
#' \describe{
#'   \item{wells}{Wells score (integer ranging from 1 to 10 indicating the degree to which PE is suspected based on clinical review of symptoms)}
#'   \item{pregnant}{Factor indicating pregnancy (No, Yes)}
#'   \item{pe}{Factor indicating pulmonary embolism has occurred (No,Yes)}
#'   \item{angio}{Result of pulmonary angiography test (Negative, Positive)}
#'   \item{d.dimer}{Numeric result of diagnostic blood test called D-Dimer.}
#'   \item{treat}{Factor indicating whether or not treatment for PE was administered (No,Yes)}
#'   \item{death}{Factor indicating patient mortality (No,Yes)}
#' }
#' @source Simulated data - not from real patients.
"PE"

#' Example Conditional Probability Table Resulting from the \code{inputCPT} function.
#'
#' This is an example of the output generated by the \code{inputCPT} function as 
#' illustrated in the article being submitted to JSS.  It is saved as an object 
#' named \code{h} in the article.
#' 
#' @source No Source.  It's really just made up.
"inputCPTExample"

#' Conditional Probability Table for side effects as a function of emesis and drug.
#'
#' This is a conditional probability table used in the emesis example of the JSS article.
#' 
"SE.cpt"

#' Conditional Probability Table for resolution of side effects as a function drugs and emesis.
#'
#' This is a conditional probability table used in the emesis example of the JSS article.
#' 
"Resolution.cpt"