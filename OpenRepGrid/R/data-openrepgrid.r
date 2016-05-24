###############################################################################
###		       data that comes along with the OpenRepGrid package			  	    ###
###############################################################################


###############################################################################

#' Grid data from Bell (2010).
#'
#' Grid data orginated (but is not shown in the paper) from a study by Haritos, 
#' Gindinis, Doan and Bell (2004) on element role titles. It was used to 
#' demonstrate the effects of construct alignment in Bell (2010, p. 46).
#'
#' @name data-bell2010
#' @aliases bell2010
#' @docType data
#' @references  Bell, R. C. (2010).	A note on aligning constructs.
#'              \emph{Personal Construct Theory and Practice}, 7, 43-48.
#'
#'              Haritos, A., Gindidis, A., Doan, C., & Bell, R. C. (2004). 
#'                The effect of element role titles on construct structure 
#'                and content. \emph{Journal of constructivist 
#'                psychology, 17}(3), 221-236.
#'
#' @keywords data
#'
NULL

# args <- list(
#   name=c( "self", "closest friend of the same sex", "the unhappiest person you know",
#           "A person of the opposite sex that you don't get along with",
#           "A teacher you respected", 
#           "Mother (or the person who filled that kind of role)",
#           "A person of the opposite sex that you like more than you dislike",
#           "The most confident person you know", 
#           "A person you work well with (at job or on sports team etc)",
#           "A teacher you did not respect"),
#   l.name=c( "relaxed", "not so smart (academically)", "dislikes sports",
#             "not interactive", "not transparent", "insensitive", "fearful&timid",
#             "rough", "accept as it is"),
#   r.name=c( "worried & tense", "smart (academically)", "loves sports",
#             "loves people", "transparent", "sensitive", "fearless", "gentle", 
#             "loves to argue"),
#   scores= c(4,4,6,5,3,6,5,2,2,6,
#             6,7,6,5,7,4,6,7,4,7,
#             6,3,7,6,4,4,2,3,6,3,
#             6,7,5,6,6,5,6,7,7,4,
#             6,4,5,7,3,7,6,5,6,3,
#             4,6,5,4,4,6,5,3,4,2,
#             5,4,4,6,5,3,5,6,5,5,
#             5,6,6,4,5,7,7,3,5,6,
#             5,5,6,7,4,4,6,7,5,5)
# )
# bell2010 <- makeRepgrid(args)
# bell2010 <- setScale(bell2010, 1, 7)
# save("bell2010", file="../data/bell2010.RData")


###############################################################################

#' Grid data from Bell and McGorry (1992).
#'
#' The grid data set is used in Bell's technical report "Using SPSS to Analyse
#' Repertory Grid Data" (1997, p. 6). Originally, the data comes from a study
#' by Bell and McGorry (1992).
#'
#' @name data-bellmcgorry1992
#' @aliases bellmcgorry1992
#' @docType data
#' @references  Bell, R. C. (1977). \emph{Using SPSS to 
#'                Analyse Repertory Grid Data}. Technical Report, 
#'                University of Melbourne.
#' 
#'              Bell, R. C., & McGorry, P. (1992). The analysis of repertory 
#'                grids used to monitor the perceptions of recovering psychotic 
#'                patients. In A. Thomson & P. Cummins (Eds.), \emph{European Perspectives 
#'                in Personal Construct Psychology} (p. 137-150). 
#'                Lincoln, UK: European Personal Construct Association.
#'
#' @keywords data
#'
NULL


# abbreviations not yet included
# c("bipolar", "schiz", "psychiat", "criminal", "average",
#   "aids", "diabetes", "cancer", "stress", "usualme", "menow",
#   "meth", "staffme", "idealme", "conlab")
# args <- list(
#   name= c('person with manic depressive illness', 'person with schizophrenia',
#           'psychiatric patient','convicted criminal','average person',
#           'AIDS patient','person with diabetes','person with cancer',
#           'person under stress','myself as I usually am','myself as I am now',
#           'myself as I will be in six months','myself as the staff see me',
#           'my ideal self'),
#   l.name=c("good", "dependable", "safe", "clearheaded", "stable", 
#            "predictable", "intelligent", "free", "healthy", "honest", 
#            "rational", "independent", "calm", "understood"),
#   r.name=rep("", 14),
#   scores =c(3,5,3,7,4,2,1,1,7,1,1,1,2,1,
#             4,5,4,7,3,3,1,2,7,2,1,1,3,1,
#             6,7,4,7,3,7,1,2,6,1,1,1,2,1,
#             7,7,4,7,2,7,2,4,7,1,1,1,2,1,
#             7,7,5,7,2,7,1,4,7,2,2,1,3,1,
#             7,7,5,7,6,5,4,4,7,7,7,3,7,7,
#             4,5,4,7,2,1,2,2,7,1,1,1,2,1,
#             7,7,5,7,1,7,1,1,1,1,7,1,7,1,
#             6,7,4,7,1,7,7,7,7,1,4,1,4,1,
#             5,5,4,7,1,1,3,4,4,1,1,1,3,1,
#             5,7,5,6,1,2,2,4,5,1,1,1,3,1,
#             5,6,5,7,1,5,3,4,5,1,1,1,2,1,
#             6,5,5,7,2,7,3,4,7,1,1,1,1,1,
#             6,7,6,7,1,7,3,4,7,3,7,1,7,1)
# )
# 
# bellmcgorry1992 <- makeRepgrid(args)
# bellmcgorry1992 <- setScale(bellmcgorry1992, 1, 7)
# save("bellmcgorry1992", file="../data/bellmcgorry1992.RData")


###############################################################################

#' Grid data from Boeker (1996).
#'
#' Grid data from a schizophrenic patient undergoing psychoanalytically 
#' oriented psychotherapy. The data was taken during the last stage of 
#' therapy (Boeker, 1996, p. 163).
#'
#' @name data-boeker
#' @aliases boeker
#' @docType data
#' @references    Boeker, H. (1996). The reconstruction of the self in 
#'                the psychotherapy of chronic schizophrenia: a case study 
#'                with the Repertory Grid Technique. In: Scheer, J. W., 
#'                Catina, A. (Eds.): \emph{Empirical Constructivism in Europe - 
#'                The Personal Construct Approach} (p. 160-167).
#'                Giessen: Psychosozial-Verlag. 
#'                
#' @keywords data
NULL

# Heinz Boeker in Scheer & Catina (1996, p.163)
# scale used: 1 to 6

# args <- list(
#  name=c("self", "ideal self", "mother", "father", "kurt", "karl",
#         "george", "martin", "elizabeth", "therapist", "irene",
#         "childhood self", "self before illness", "self with delusion" ,
#         "self as dreamer"),
#  l.name=c("balanced", "isolated", "closely integrated", "discursive",
#           "open minded", "dreamy", "practically oriented", "playful",
#           "socially minded", "quarrelsome", "artistic", "scientific",
#           "introvert", "wanderlust"),
#  r.name=c("get along with conflicts", "sociable", "excluded", "passive",
#           "indifferent", "dispassionate", "depressed", "serious", 
#           "selfish", "peaceful", "technical", "emotional",
#           "extrovert", "home oriented"),
#  scores= c(1, 4, 2, 2, 3, 5, 2, 5, 4, 2, 6, 2, 2, 3, 3,
#            3, 6, 3, 5, 5, 4, 5, 4, 5, 4, 4, 4, 2, 2, 3,
#            2, 2, 2, 3, 5, 3, 2, 3, 2, 3, 3, 4, 4, 5, 3,
#            4, 1, 3, 1, 2, 4, 2, 3, 3, 2, 3, 3, 3, 5, 4,
#            2, 1, 2, 1, 2, 4, 4, 2, 4, 2, 6, 3, 2, 2, 3,
#            4, 5, 3, 5, 4, 5, 4, 5, 4, 4, 6, 3, 3, 3, 2,
#            2, 1, 3, 2, 3, 3, 3, 2, 2, 3, 2, 3, 3, 3, 3,
#            4, 5, 4, 3, 4, 3, 2, 3, 4, 4, 5, 3, 2, 4, 3,
#            2, 1, 3, 2, 4, 5, 4, 1, 3, 2, 6, 3, 3, 3, 3,
#            5, 5, 5, 5, 5, 2, 5, 2, 4, 4, 1, 6, 5, 5, 5,
#            5, 1, 2, 4, 3, 5, 3, 2, 4, 3, 3, 4, 4, 4, 4,
#            2, 1, 5, 3, 4, 4, 5, 3, 4, 1, 6, 4, 2, 3, 3,
#            4, 5, 4, 6, 5, 3, 5, 3, 5, 2, 5, 2, 2, 2, 3,
#            1, 1, 4, 2, 4, 5, 2, 5, 5, 3, 6, 1, 1, 2, 1)
# )
# boeker <- makeRepgrid(args)
# boeker <- setScale(boeker, 1, 6)
# save("boeker", file="../data/boeker.RData")


###############################################################################

#' Grid data from Fransella, Bell and Bannister (2003).
#'
#' A dataset used throughout the book "A Manual for Repertory Grid Technique" 
#' (Fransella, Bell and Bannister, 2003, p. 60). 
#'
#' @name data-fbb2003
#' @aliases fbb2003
#' @docType data
#' @references  Fransella, F., Bell, R. & Bannister, D. (2003). A Manual for Repertory 
#'              Grid Technique (2. Ed.). Chichester: John Wiley & Sons.
#' @keywords data
#'
NULL


# args <- list(
#   name= c("self", "my father", "an old flame", "an ethical person", "my mother", 
#           "a rejected teacher", "as I would love to be", "a pitied person"),
#   l.name=c("clever", "disorganized", "listens", "no clear view", "understands me",
#            "ambitious", "respected", "distant", "rather aggressive"),
#   r.name=c("not bright", "organized", "doesn't hear", "clear view of life", "no understanding",
#            "no ambition", "not respected", "warm", "not aggressive"),
#   scores =c(2,1,6,3,5,7,1,5,
#             6,6,4,5,2,2,5,2,
#             3,1,6,3,3,7,1,4,
#             5,6,3,3,3,5,7,3,
#             3,2,6,2,2,6,2,5,
#             6,3,5,4,7,3,3,5,
#             2,2,4,2,5,6,1,4,
#             3,3,7,3,5,1,6,5,
#             1,3,3,3,5,2,5,7)
# )
# 
# fbb2003 <- makeRepgrid(args)
# fbb2003 <- setScale(fbb2003, 1, 7)
# save("fbb2003", file="../data/fbb2003.RData")


###############################################################################

#' Grid data from Feixas and Saul (2004). 
#'
#'A desription by the authors:
#' "When Teresa, 22 years old, was seen by the second author (LAS) at the 
#' psychological services of the University of Salamanca, she was in the 
#' final year of her studies in chemical sciences. Although Teresa proves 
#' to be an excellent student, she reveals serious doubts about her self worth. 
#' She cries frequently, and has great difficulty in meeting others, 
#' even though she has a boyfriend who is extremely supportive. Teresa 
#' is anxiously hesitant about accepting a new job which would involve moving 
#' to another city 600 Km away from home." (Feixas & Saul, 2004, p. 77).
#'
#' @name data-feixas2004
#' @aliases feixas2004
#' @docType data
#' @references  Feixas, G., & Saul, L. A. (2004). The Multi-Center Dilemma Project: 
#'              an investigation on the role of cognitive conflicts in health. 
#'              \emph{The Spanish Journal of Psychology, 7}(1), 69-78.
#'
#' @keywords data
#'
NULL

# args <- list(
#   name= c("Self now", "Mother", "Father", "Brother", "Boyfriend", "Friend 1",
#           "Friend 2", "Non-grata", "Friend 3", "Cousin", "Godmother", "Friend 4", 
#           "Ideal Self"),
#   l.name= c("Pessimistic", "Self-demanding", "Fearful", "Lives to work", 
#             "Imposes his/her wishes", "Teasing", "Appreciates others", "Aggressive",
#             "Concerned about others", "Avaricious", "Sensitive", "Cheeky", "Hypocritical",
#             "Blackmailer", "Appears stronger than is", "Does not look after the friendship",
#             "Non Accessible", "Introverted", "Gets depressed easily", 
#             "Tries to find the good in things"),
#   r.name= c("Optimistic", "Takes it easy", "Enterprising", "Works to live", 
#             "Tolerant with others", "Touchy", "Does not appreciate others",
#             "Calm", "Selfish", "Generous", "Materialistic, superficial", "Respectful",
#             "Sincere", "Non blackmailer", "Natural", "Looks after the friendship",
#             "Accessible", "Extroverted", "Does not get depressed easily", 
#             "Sees only the negative"),
#   scores= c(1, 1, 5, 2, 7, 3, 6, 2, 6, 4, 3, 2, 7,
#             1, 6, 6, 2, 2, 5, 6, 3, 5, 6, 4, 5, 4,
#             2, 2, 6, 2, 4, 5, 6, 5, 2, 3, 4, 5, 5,
#             5, 1, 2, 2, 6, 6, 6, 1, 6, 7, 6, 6, 7,
#             6, 2, 1, 1, 4, 3, 6, 1, 7, 3, 4, 2, 7,
#             2, 7, 1, 6, 4, 3, 4, 6, 3, 3, 5, 6, 3,
#             2, 6, 6, 6, 1, 5, 4, 7, 4, 2, 2, 5, 1,
#             6, 4, 2, 2, 7, 4, 6, 2, 6, 6, 6, 3, 7,
#             2, 2, 6, 7, 2, 3, 5, 7, 3, 3, 2, 2, 2,
#             6, 1, 1, 1, 7, 5, 5, 1, 6, 3, 3, 6, 7,
#             1, 5, 7, 7, 1, 4, 5, 7, 1, 4, 3, 4, 1,
#             6, 6, 5, 4, 6, 6, 6, 1, 6, 5, 6, 5, 7,
#             5, 4, 4, 2, 6, 5, 5, 1, 6, 6, 5, 4, 7,
#             3, 2, 2, 1, 5, 6, 6, 1, 6, 6, 6, 3, 7,
#             6, 3, 1, 2, 5, 2, 4, 2, 7, 6, 6, 5, 6,
#             6, 3, 3, 3, 6, 2, 1, 2, 4, 4, 6, 4, 7,
#             5, 2, 2, 1, 4, 2, 4, 1, 6, 3, 5, 2, 7,
#             1, 2, 6, 2, 4, 5, 7, 5, 2, 6, 6, 5, 5,
#             1, 2, 6, 3, 6, 3, 7, 6, 1, 3, 3, 3, 6,
#             6, 6, 4, 6, 1, 5, 2, 7, 6, 3, 3, 5, 1)
# )
# 
# feixas2004 <- makeRepgrid(args)
# feixas2004 <- setScale(feixas2004, 1, 7)
# save("feixas2004", file="../data/feixas2004.RData")









###############################################################################

#' Case as described by the authors:
#' "Sarah, aged 32, was referred with problems of depression and sexual 
#' difficulties relating to childhood sexual abuse. She had three 
#' children and was living with her male partner.
#' From the age of 9, her brother, an adult, had sexually abused Sarah. 
#' She attended a group for survivors of child sexual abuse and 
#' completed repertory grids prior to the group, immediately after the 
#' group and at 3- and 6-month follow-up." (Leach et al. 2001, p. 230).\cr \cr
#' 
#' \code{leach2001a} is the pre-therapy, \code{leach2001b}
#' is the post-therapy therapy datset. The construct and elements are 
#' identical.
#'
#' @title   Pre- and post therapy dataset from Leach et al. (2001).
#'
#' @name data-leach2001
#' @aliases leach2001a leach2001b
#' @docType data
#'
#' @references      Leach, C., Freshwater, K., Aldridge, J., & 
#'                  Sunderland, J. (2001). Analysis of repertory grids 
#'                  in clinical practice. \emph{The British Journal 
#'                  of Clinical Psychology, 40}, 225-248.
#'
#' @keywords data
#'
NULL

# name.abb <- c("CS", "SN", "WG", "MG", "Fa", "Pa", "IS", "Mo", "AC") # not included yet
# args <- list(
#   name= c("Child self", "Self now",  "Women in general", 
#           "Men in general", "Father", "Partner", "Ideal self",
#           "Mother", "Abuser in childhood"),
#   l.name= c("assertive", "confident", "does not feel guilty", "abusive",
#             "frightening", "untrustworthy", "powerful", "big headed",
#             "independent", "confusing", "guilty", "cold", "masculine",
#             "interested in sex"),
#   r.name= c("not assertive", "unconfident", "feels guilty", "not abusive", 
#             "not frightening", "trustworthy", "powerless", "not big headed", 
#             "dependent", "not confusing", "not guilty", "shows feelings", 
#             "feminine", "not interested in sex"),
#   scores= c(2,7,4,2,3,5,3,1,1,
#             1,7,3,2,2,4,2,1,1,
#             1,6,4,2,1,1,1,1,1,
#             7,7,4,6,7,6,7,3,1,
#             7,7,4,5,7,7,7,3,2,
#             7,7,6,5,7,7,7,3,1,
#             7,5,4,2,3,5,2,1,1,
#             7,5,4,2,6,6,4,2,1,
#             5,6,3,2,2,4,1,3,1,
#             7,2,4,4,7,6,7,1,2,
#             7,3,4,4,7,6,7,4,1,
#             7,3,5,4,7,7,6,2,6,
#             7,7,5,1,1,2,5,2,1,
#             7,5,3,1,1,1,2,7,1 ))
# leach2001a <- makeRepgrid(args)
# leach2001a <- setScale(leach2001a, 1, 7)
# save("leach2001a", file="../data/leach2001a.RData")


###############################################################################


# name.abb <- c("CS", "SN", "WG", "MG", "Fa", "Pa", "IS", "Mo", "AC") # not included yet
# args <- list(
#   name= c("Child self", "Self now",  "Women in general", 
#            "Men in general", "Father", "Partner", "Ideal self",
#            "Mother", "Abuser in childhood"),
#   l.name= c("assertive", "confident", "does not feel guilty", "abusive",
#             "untrustworthy", "guilty", "big headed", "frightening",
#             "cold", "powerful", "confusing", "not interested in sex",
#             "dependent", "masculine"),
#   r.name= c("not assertive", "unconfident", "feels guilty", "not abusive",
#             "trustworthy", "not guilty", "not big headed",
#             "not frightening", "shows feelings", "powerless", 
#             "not confusing", "interested in sex", "independent",
#             "feminine"),
#   scores= c( 4,5,5,3,6,6,2,1,1,
#              3,6,4,3,3,5,1,1,1,
#              2,4,4,2,1,2,1,1,1,
#              7,5,4,4,7,7,7,3,1,
#              7,7,6,5,7,7,7,3,1,
#              7,7,4,4,7,7,7,5,1,
#              6,6,4,4,7,5,7,4,1,
#              7,6,4,4,7,7,7,2,4,
#              5,6,6,4,7,7,7,2,5,
#              6,3,3,2,3,5,1,1,1,
#              7,3,6,6,6,6,7,1,3,
#              1,4,4,4,5,6,6,1,7,
#              3,2,6,5,5,4,6,3,6,
#              6,6,7,1,2,1,7,4,1 ))
# leach2001b <- makeRepgrid(args)
# leach2001b <- setScale(leach2001b, 1, 7)
# save("leach2001b", file="../data/leach2001b.RData")



###############################################################################

#' Grid data from Mackay (1992). Data set 'Grid C'-
#'
#' 
#' used in Mackay's paper on inter-element correlation 
#' (1992, p. 65).
#'
#' @name data-mackay1992
#' @aliases mackay1992
#' @docType data
#' @references  Mackay, N. (1992). Identification, reflection, 
#'              and correlation: Problems in the bases of repertory 
#'              grid measures. \emph{International Journal of Personal
#'              Construct Psychology, 5}(1), 57-75. 
#'
#' @keywords data
#'
NULL


# args <- list(
#   name= c("Self", "Ideal self", "Mother", "Father", "Spouse",
#           "Disliked person"),
#   l.name=c("Quick", "*Satisfied", "Talkative", "*Succesful", 
#            "Emotional", "*Caring"),
#   r.name=c("*Slow", "Bitter", "*Quiet", "Loser", "*Calm", 
#            "Selfish"),
#   scores =c(7,4,7,5,3,3,
#             7,7,3,5,2,3,
#             6,4,6,5,5,2,
#             6,7,2,2,3,2,
#             5,6,5,4,4,1,
#             4,7,6,4,5,1)
# )
# 
# mackay1992 <- makeRepgrid(args)
# mackay1992 <- setScale(mackay1992, 1, 7)
# save("mackay1992", file="../data/mackay1992.RData")


###############################################################################

#' Grid data from Raeithel (1998).
#' 
#' Grid data to demonstrate the use of Bertin diagrams (Raeithel, 1998, p. 223).
#' The context of its administration is unknown.
#' 
#' @name data-raeithel
#' @aliases raeithel
#' @docType data
#' @references   Raeithel, A. (1998). Kooperative Modellproduktion von Professionellen 
#'          und Klienten. Erlaeutert am Beispiel des Repertory Grid.
#'          In A. Raeithel (1998). Selbstorganisation, Kooperation, 
#'          Zeichenprozess. Arbeiten zu einer kulturwissenschaftlichen, 
#'          anwendungsbezogenen Psychologie (p. 209-254). Opladen: 
#'          Westdeutscher Verlag.
#'
#' @keywords data
NULL

# args <- list(
#  name=c("Freund", "Therapeut", "Ideal", "Arzt", "Ich vorher", "Virus", 
#         "Mutter", "Schwester", "Partner", "Vater", "Neg. Person", "Ich",
#         "Freundin"),
#  l.name=c("charakterlos", "uninteressiert", "sich gehen lassen", 
#           "unerfahren", "abschaetzbar", "mutig", "leichtfuessig", 
#           "kuehl", "egoistisch eigennuetzig", "angepasst", "offen",
#           "unvorsichtig", "bruederlich freundschaftlich"),      
#  r.name=c("Charakter haben", "interessiert", "zielstrebig", "erfahren", 
#           "unberechenbar", "feige", "schwerfaellig", "warmherzig", 
#           "lebensbejahend sozial", "unangepasst", "hinterhaeltig", 
#           "diszipliniert gesund", "vaeterlich autoritaer"),
#  scores= c(  1,  1,  1,  1,  1, -1, -1, -1,  1, -1, -1,  1,  1,
#              1,  1,  1,  1,  1,  1, -1, -1,  1,  1, -1,  1,  1, 
#              1,  1,  1,  1,  1,  1, -1, -1,  1, -1, -1,  1,  1, 
#              1,  1,  1,  1,  1,  1, -1, -1, -1,  0, -1,  1,  1,
#             -1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1,  1,  0, 
#             -1, -1, -1, -1,  1,  1,  1,  1, -1,  1,  1,  1, -1,
#             -1, -1, -1, -1, -1, -1,  1,  1, -1,  1,  0,  1, -1, 
#              1,  1,  1,  0, -1, -1, -1,  1,  1,  1,  1,  0, -1, 
#              1,  1,  1,  1, -1, -1, -1,  0,  1,  1, -1, -1,  0, 
#              1, -1,  1,  0,  1,  1, -1, -1, -1, -1, -1, -1, -1,
#             -1, -1, -1, -1,  1,  1,  1, -1, -1, -1,  1, -1, -1, 
#             -1, -1,  1,  1, -1,  1,  0,  1, -1,  1, -1,  1, -1,
#              0, -1, -1,  1, -1,  0,  0, -1,  0,  1,  0,  0, -1)
# )
# 
# raeithel <- makeRepgrid(args)
# raeithel <- setScale(raeithel, -1, 1)
# save("raeithel", file="../data/raeithel.RData")


###############################################################################

#' Drug addict's grid data set from Slater (1977, p. 32).
#'
#' @name data-slater1977a
#' @aliases slater1977a
#' @docType data
#' @references  Slater, P. (1977). \emph{The measurement of intrapersonal space 
#'              by grid technique}. London: Wiley.
#'
#' @keywords data
#'
NULL

# args <- list(
#   name=c("Alcohol", "Barbiturates", "Cannabis", "Cocain", "Drynomil", 
#          "Heroin", "L.S.D.", "Madryx", "Methedrine (injections)"),
#   r.name=c( "Makes me talk more",
#             "Makes me feel high",
#             "Makes me feel blocked",
#             "Makes me feel sleepy",
#             "Gives me a warminn feeling inside",
#             "Makes me feel drunk",
#             "Makes me imagine things",
#             "Makes me feel sick",
#             "Makes me do things without knowing what I'm doing",
#             "Hepls me enjoy things",
#             "Gives me a good buzz",
#             "Makes me tense",
#             "Makes me feel sexy",
#             "After taking it I may see or hear people whoe aren't really there"),
#   l.name=rep("", 14),
#   scores=c( 2,2,4,1,1,4,5,2,1,
#             3,1,2,1,1,2,1,1,2,
#             3,1,2,1,1,2,1,1,2,
#             3,1,2,5,5,2,5,1,5,
#             2,1,2,1,4,2,3,1,3,
#             2,1,2,5,5,2,5,1,4,
#             3,3,1,3,2,2,1,3,2,
#             3,3,2,3,3,3,3,3,3,
#             3,1,4,3,2,3,4,2,2,
#             3,2,2,2,1,2,1,2,4,
#             3,1,1,1,3,1,1,1,3,
#             3,5,2,1,1,5,3,5,1,
#             2,5,2,4,3,5,3,5,4,
#             3,3,3,2,1,3,2,3,2)
# )
# 
# slater1977a <- makeRepgrid(args)
# slater1977a <- setScale(slater1977a, 1, 5)
# save("slater1977a", file="../data/slater1977a.RData")


###############################################################################

#' Grid data from Slater (1977). 
#'
#' Grid data (ranked) from a seventeen year old female psychiatric patient 
#' (Slater, 1977, p. 110). She was depressed, anxious and took to cutting 
#' herself. The data was originally reported by Watson (1970).
#'
#' @name data-slater1977b
#' @aliases slater1977b
#' @docType data
#' @references  Slater, P. (1977). \emph{The measurement of intrapersonal space 
#'              by grid technique}. London: Wiley.
#'
#'              Watson, J. P. (1970). The relationship between a self-mutilating 
#'              patient and her doctor. \emph{Psychotherapy and Psychosomatics,
#'              18}(1), 67-73.
#'
#' @keywords data
#'
NULL

# args <- list(
#   name= c("Wanting to talk to someone and being unable to", 
#           "Having the same thoughts for a long time",
#           "Being in a crowd", "Seeing G.", "Being at home",
#           "Being in hospital", "Being with my mother",
#           "Being with Dr. W.", "Being with Mrs. M.", "Being with my father"),
#   r.name= c("Make me cut myself", "Make me think people are unfriendly",
#             "Make me feel depressed", "Make me feel angry", 
#             "Make me feel scared", "Make me feel more grown-up", 
#             "Make me feel more like a child", "Make me feel lonely", 
#             "Help me in the long run", "Make me feel cheerful"),
#   l.name= rep("", 10),
#   scores= matrix(c(2,1,3,6,4,5,7,8,10,9,
#                  1,3,6,2,4,7,5,8,10,9,
#                  2,5,3,1,4,6,7,10,9,8,
#                  1,2,4,3,7,6,10,5,8,9,
#                  2,3,5,1,9,7,10,4,6,8,
#                  5,4,7,1,2,6,3,9,10,8,
#                  2,7,5,9,1,8,3,6,10,4,
#                  9,8,7,1,5,6,4,2,10,3,
#                  5,9,10,1,8,2,6,3,4,7,
#                  9,8,10,1,5,6,3,4,7,2)))
# slater1977b <- makeRepgrid(args)
# slater1977b <- setScale(slater1977b, 1, 10)
# save("slater1977b", file="../data/slater1977b.RData")


###############################################################################











