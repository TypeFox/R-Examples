## ---- echo=FALSE, message=FALSE------------------------------------------
library(HydeNet)

## ---- fig.width=7, eval = 1----------------------------------------------
net <- HydeNetwork(~ initialAces | card1*card2
                   + initialPoints | card1*card2
                   + highUpcard | dealerUpcard
                   + hit1 | initialPoints*highUpcard
                   + acesAfterCard3 | initialAces*card3
                   + pointsAfterCard3 | card1*card2*card3*acesAfterCard3
                   + hit2 | pointsAfterCard3*highUpcard
                   + acesAfterCard4 | acesAfterCard3*card4
                   + pointsAfterCard4 | card1*card2*card3*card4*acesAfterCard4
                   + hit3 | pointsAfterCard4*highUpcard
                   + acesAfterCard5 | acesAfterCard4*card5
                   + pointsAfterCard5 | card1*card2*card3*card4*card5*acesAfterCard5
                   + playerFinalPoints | initialPoints*hit1*pointsAfterCard3
                                         *hit2*pointsAfterCard4*hit3*pointsAfterCard5
                   + dealerOutcome | dealerUpcard
                   + payoff | playerFinalPoints*dealerOutcome)

plot(net)

## ------------------------------------------------------------------------
#####################################
# Random Variable Nodes
#####################################

# Note: The calls to 'setNode' below are going to generate messages
#       that validation has been ignored.  'setNode' assumes that when
#       you pass a character string as a parameter definition, you 
#       are passing JAGS code.  See the Validation section of the 
#       'setNode' documentation for more information.
#       After the first call, we'll turn off validation

cardProbs  <- c(rep(1/13,8), 4/13, 1/13)  # probs. for 2, 3, ..., 9, (10-K), A

net <- setNode(net, card1, nodeType="dcat",  pi=vectorProbs(p=cardProbs, card1))
net <- setNode(net, card2, nodeType="dcat",  pi=vectorProbs(p=cardProbs, card2),
               validate=FALSE)
net <- setNode(net, card3, nodeType="dcat",  pi=vectorProbs(p=cardProbs, card3),
               validate=FALSE)
net <- setNode(net, card4, nodeType="dcat",  pi=vectorProbs(p=cardProbs, card4),
               validate=FALSE)
net <- setNode(net, card5, nodeType="dcat",  pi=vectorProbs(p=cardProbs, card5),
               validate=FALSE)

net <- setNode(net, dealerUpcard, nodeType="dcat",
               pi=vectorProbs(p=cardProbs, dealerUpcard),
               validate=FALSE)

#Note: node dealerOutcome will be defined below, following some discussion 
#      about its conditional probability distribution.

#####################################
# Deterministic Nodes
#####################################
net <- setNode(net, highUpcard,     "determ", define=fromFormula(),
               nodeFormula = highUpcard ~ ifelse(dealerUpcard > 8, 1, 0))

net <- setNode(net, initialAces,    "determ", define=fromFormula(),
               nodeFormula = initialAces ~ ifelse(card1==10,1,0) + ifelse(card2==10,1,0))

net <- setNode(net, acesAfterCard3, "determ", define=fromFormula(),
               nodeFormula = acesAfterCard3 ~ initialAces + ifelse(card3==10,1,0))

net <- setNode(net, acesAfterCard4, "determ", define=fromFormula(),
               nodeFormula = acesAfterCard4 ~ acesAfterCard3 + ifelse(card4==10,1,0))

net <- setNode(net, acesAfterCard5, "determ", define=fromFormula(),
               nodeFormula = acesAfterCard5 ~ acesAfterCard4 + ifelse(card5==10,1,0))
               
net <- setNode(net, initialPoints, "determ", define=fromFormula(),
               nodeFormula = initialPoints ~ card1+card2+2)

net <- setNode(net, pointsAfterCard3, "determ", define=fromFormula(),
               nodeFormula = pointsAfterCard3 ~
                 ifelse(acesAfterCard3 == 3,
                    13,
                    ifelse(acesAfterCard3 == 2,
                       card1 + card2 + card3 + 3 - 10,
                       ifelse(acesAfterCard3 == 1,
                          ifelse(card1 + card2 + card3 + 3 > 22,
                             card1 + card2 + card3 + 3 - 10,
                             card1 + card2 + card3 + 3),
                          card1 + card2 + card3 + 3
                       )
                    )
                 )
              )

net <- setNode(net, pointsAfterCard4, "determ", define=fromFormula(),
               nodeFormula = pointsAfterCard4 ~
                 ifelse(acesAfterCard4 == 4,
                    14,
                    ifelse(acesAfterCard4 == 3,
                       ifelse(card1 + card2 + card3 + card4 + 4 > 38,
                          card1 + card2 + card3 + card4 + 4 - 30,
                          card1 + card2 + card3 + card4 + 4 - 20
                       ),
                       ifelse(acesAfterCard4 > 0,
                          ifelse(card1 + card2 + card3 + card4 + 4 > 22,
                                 card1 + card2 + card3 + card4 + 4 - 10,
                                 card1 + card2 + card3 + card4 + 4
                          ),
                          card1 + card2 + card3 + card4 + 4
                       )
                    )
                 )
              )

net <- setNode(net, pointsAfterCard5, "determ", define=fromFormula(),
               nodeFormula = pointsAfterCard5 ~ 
                 ifelse(acesAfterCard5 == 5,
                    15,
                    ifelse(acesAfterCard5 == 4,
                       ifelse(card1 + card2 + card3 + card4 + card5 + 5 > 51,
                          card1 + card2 + card3 + card4 + card5 + 5 - 40,
                          card1 + card2 + card3 + card4 + card5 + 5 - 30
                        ),
                       ifelse(acesAfterCard5 == 3,
                          ifelse(card1 + card2 + card3 + card4 + card5 + 5 > 51,
                             card1 + card2 + card3 + card4 + card5 + 5 - 30,
                             card1 + card2 + card3 + card4 + card5 + 5 - 20
                          ),
                          ifelse(acesAfterCard5 == 2,
                             ifelse(card1 + card2 + card3 + card4 + card5 + 5 > 31,
                                card1 + card2 + card3 + card4 + card5 + 5 - 20,
                                card1 + card2 + card3 + card4 + card5 + 5 - 10
                             ),
                             ifelse(acesAfterCard5 > 0,
                                ifelse(card1 + card2 + card3 + card4 + card5 + 5 > 22,
                                   card1 + card2 + card3 + card4 + card5 + 5 - 10,
                                   card1 + card2 + card3 + card4 + card5 + 5
                                ),
                                card1 + card2 + card3 + card4 + card5 + 5
                             )
                          )
                       )
                    )
                 )
              )

net <- setNode(net, playerFinalPoints, "determ", define=fromFormula(),
               nodeFormula = playerFinalPoints ~ 
                 ifelse(hit1 == 0,
                    initialPoints,
                    ifelse(hit2 == 0,
                           pointsAfterCard3,
                           ifelse(hit3 == 0, pointsAfterCard4, pointsAfterCard5)
                    )
                 )
              )


## ------------------------------------------------------------------------
data(BJDealer)
dealerOutcome.cpt <- cpt(dealerOutcome ~ dealerUpcard,
                         data = BJDealer,
                         wt   = BJDealer$probability)
round(dealerOutcome.cpt,3)
net <- setNodeModels(net, dealerOutcome.cpt)

## ------------------------------------------------------------------------
net <- setDecisionNodes(net, hit1, hit2, hit3)
net <- setUtilityNodes(net, payoff)

c(net$nodeDecision$hit2, net$nodeUtility$payoff)

## ---- fig.width=7, eval=FALSE--------------------------------------------
#  plot(net)

## ------------------------------------------------------------------------
data(BlackJackTrain)
BlackJackTrain$highUpcard <- as.character(BlackJackTrain$dealerUpcard)
BlackJackTrain$highUpcard <- factor(BlackJackTrain$highUpcard %in% c("10-K","A"),
                                    c(FALSE, TRUE), c("9 or lower", "10 or higher"))

# glm.hit1 <- glm(hit1 ~ initialPoints + I(dealerUpcard %in% c("9","10-K","A")),
#                 data = BlackJackTrain, family="binomial")

glm.hit1 <- glm(hit1 ~ initialPoints+highUpcard,
                data = BlackJackTrain, family="binomial")
glm.hit2 <- glm(hit2 ~ pointsAfterCard3+highUpcard,
                data = BlackJackTrain, family="binomial")
glm.hit3 <- glm(hit3 ~ pointsAfterCard4+highUpcard,
                data = BlackJackTrain, family="binomial")

net <- setNodeModels(net, glm.hit1, glm.hit2, glm.hit3)

## ------------------------------------------------------------------------
net <- setNode(net, payoff, "determ", define=fromFormula(),
         nodeFormula = payoff ~
                         ifelse(playerFinalPoints > 21, -1,
                           ifelse(playerFinalPoints == 21,
                             ifelse(dealerOutcome == "Blackjack", 0,
                               ifelse(dealerOutcome == 7, 0, 1)),
                             ifelse(dealerOutcome == "Bust",
                               ifelse(playerFinalPoints < 22, 1, -1),
                               ifelse(dealerOutcome == "17",
                                 ifelse(playerFinalPoints == 17, 0,
                                 ifelse(playerFinalPoints > 17, 1, -1)),
                                 ifelse(dealerOutcome == "18",
                                   ifelse(playerFinalPoints == 18, 0,
                                     ifelse(playerFinalPoints > 18, 1, -1)),
                                   ifelse(dealerOutcome == "19",
                                     ifelse(playerFinalPoints == 19, 0,
                                       ifelse(playerFinalPoints > 19, 1, -1)),
                                     ifelse(dealerOutcome == "20",
                                       ifelse(playerFinalPoints == 20, 0,
                                         ifelse(playerFinalPoints > 20, 1, -1)),
                                       ifelse(playerFinalPoints == 21, 0, -1)))))))))


## ------------------------------------------------------------------------
trackedVars <- c("dealerUpcard","playerFinalPoints","dealerOutcome","payoff")
evidence <- list(card1 = 3)
compiledNet <- compileJagsModel(net, data = evidence,
                                n.chains = 3,
                                n.adapt = 5000)

post <- HydePosterior(compiledNet,
                      variable.names = trackedVars,
                      n.iter=10000)

dplyr::sample_n(post, 20)

## ------------------------------------------------------------------------
table(post$payoff)
mean(post$payoff)

## ------------------------------------------------------------------------
policies <- data.frame(hit1 = c(0,1,1,1),
                       hit2 = c(0,0,1,1),
                       hit3 = c(0,0,0,1))

## ---- echo=FALSE---------------------------------------------------------
set.seed(39482820)

## ------------------------------------------------------------------------
compiledNets <- compileDecisionModel(net, policyMatrix = policies)

samples <- lapply(compiledNets,
                  HydePosterior,
                  variable.names = trackedVars,
                  n.iter=10000)

lapply(samples, head)

## ------------------------------------------------------------------------
#summary of expected utility under each policy
lapply(samples, function(l) mean(l$payoff))

