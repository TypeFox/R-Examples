################################################################################
#
#  This file contains functions that build and manipulate deck data frames
#
#  Functions beginning with "iwp" are called by the stats functions
#
#  At this time no draw poker variants are supported, this data set and
#  the wpdeal program would need to be modified to capture discard options
#  per betting round, and establish a discard "player" option, leaving aside
#  all considerations of encoding a reasonable discard strategy
#
#  Note that games like Iron Cross or Omaha will limit which cards can
#  be used to make a hand, and games that specify wild cards, add wild cards
#  or add cards will have that enforced in the wpdeal program.
#


################################################################################## 
namevals <- function(rnames, flag = "Facevalue") {
  #  This function strips values from rownames assigns values for matching, in 
  #  situations where the actual value or suit is obscured by wildcard notation 
  #
  #  Args:    a vector of rownames   eg c("AH", "10S")
  #
  #  Returns: a vector of Values     eg c("14", "10"), 
  #        or a vector of Suits      eg c("H", "S") 
  #        or a vector of Facevalue  eg c("A", "10")
  #
  switch(flag,
    "Suits" = {
      rnames <- ifelse(substring(rnames,1,2) == "10",
                  substring(rnames,3,3),
                    substring(rnames,2,2))
    },
    "Values" =  { 
      rnames <- substring(rnames,1,1)
      rnames <- ifelse(rnames == "A", "14", 
                  ifelse(rnames == "K", "13",
                    ifelse(rnames == "Q", "12", 
                      ifelse(rnames == "J", "11", 
                        ifelse(rnames == "1", "10", rnames)))))
    },
    "Facevalue" = { # default is Facevalue
      rnames <- ifelse(substring(rnames,1,2) == "10",
                  substring(rnames,1,2),
                    substring(rnames,1,1))
    }
  ) # end suit parameter check
  rnames
} # end internal function namevals

################################################################################## 
valface <- function(rvals) {
  #  This function strips values from rownames assigns values for matching, in 
  #  situations where the actual value or suit is obscured by wildcard notation 
  #
  #  Args:    a vector of values         eg c("14", "11", "10", "8")
  #
  #  Returns: a vector of facevalues     eg c("A", "J", "10", "8")
  ifelse(rvals == "14", "A", 
    ifelse(rvals == "13", "K",
      ifelse(rvals == "12", "Q", 
        ifelse(rvals == "11", "J", rvals))))
} # end internal function vlface

################################################################################
iwpshuffle <- function(wcards = NULL, shuffle = TRUE, deck = iwpstddeck) {
  #  This function creates and shuffles a deck, flagging all wildcards
  #
  # Args:
  #   wcards: a vector with string values either matching the 52 cards
  #           in a normal deck, or several common designations.  Supported
  #           wildcard designations are:
  #           "Suicide King":    = "KH"
  #           "One Eyed Jacks":  = c("JH", "JS")
  #           "Deuces":          = c("2H", "2C", "2D", "2S")
  #           "Heinz 57":        = c("5H", "5C", "5D", "5S",
  #                                  "7H", "7C", "7D", "7S")
  #           "Pregnant Threes": = c("3H", "3C", "3D", "3S",
  #                                  "6H", "6C", "6D", "6S",
  #                                  "9H", "9C", "9D", "9S")
  #           "Dr Pepper":       = c("2H", "2C", "2D", "2S",
  #                                  "4H", "4C", "4D", "4S",
  #                                  "10H", "10C", "10D", "10S")
  #           other cards are designated by "number""suit"
  #             2:9, 10, J, Q, K, A are the "numbers"
  #             D = Diamond, H = Heart, S = Spade, C = Club    
  #  Jokers are not included as an option because they have no conversion
  #  to a strict number or suit, which some utility functions require  
  #                                
  # 
  #  Returns: a shuffled 52 card deck data frame with following format
  #   row name: card description (eg, "AS" = Ace of Spades)
  #   $value    2:10, J=11, Q=12, K=13, A=14, wildcards are all = "W"
  #   $suit     D, H, S, C as in Card description unless wild (= "W")
  #   $hand     Starts out as "Deck" for all values
  #
 
  
  # deck size might be less than 52 for "Stripped" deck game variants
  dsize <- dim(deck)[1]

  # Handle other wildcard variables
  if (length(wcards) > 0) {
    if (sum(wcards == "Suicide King") > 0) {
      wcards <- c(wcards[(wcards != "Suicide King")],"KH")
    } # end of Suicide King check 
    if (sum(wcards == "One Eyed Jacks") > 0) {
      wcards <- c(wcards[(wcards != "One Eyed Jacks")],"JH","JS")
    } # end of One Eyed Jacks check 
    if (sum(wcards == "Deuces") > 0) {
      wcards <- c(wcards[(wcards != "Deuces")],"2H","2S","2C","2D")
    } # end of Deuces check 
    if (sum(wcards == "Heinz 57") > 0) {
      wcards <- c(wcards[(wcards != "Heinz 57")],"5H", "5S", "5C", "5D",
                                                 "7H", "7C", "7D", "7S")
    } # end of Heinz 57 check 
    if (sum(wcards == "Pregnant Threes") > 0) {
      wcards <- c(wcards[(wcards != "Pregnant Threes")],"3H","3S","3C","3D",
                             "6H", "6C", "6D", "6S", "9H", "9C", "9D", "9S")
    } # end of Pregnant Threes check 
    if (sum(wcards == "Dr Pepper") > 0) {
      wcards <- c(wcards[(wcards != "Dr Pepper")],"2H","2S","2C","2D",
                "4H", "4C", "4D", "4S", "10H", "10C", "10D", "10S")
    } # end of Dr Pepper check 

    # remove duplicates and any unsupported values
    wcards <- suppressWarnings(wcards[is.na(as.numeric(wcards))])
    wcards <- unique(wcards[rownames(deck[wcards, ]) != "NA"])
    # set both value and suit to W for wildcards
    deck[wcards,1:2] <- "W"
  } # End of wilcard logic

  
  if (shuffle) { # return the shuffled deck
    deck[sample(1:dsize,dsize, replace=FALSE), ]
  } else {       # don't shuffle the deck  
    deck[1:dsize, ]   
  } # end deck shuffle check
} #  end of internal function iwpdeck

################################################################################

iwpdeal <- function(ngame, players, wcards = NULL, drounds = 5, deck = NULL) {
  #
  #  Args:
  #    ngame:    any rowname from ?wpsupportedgames 
  #    players: Number greater than zero, less than max players for a game
  #    wcard:   If additional wild cards are desired, see wpgame for
  #             how to populate this field
  #    drounds: Dealing round - will deal until this round is complete.
  #             Currently the maximum dealing rounds in any game is 5, and
  #             this limit is hardcoded into this function & data frame wpgames
  #    deck:    if null, calls wpshuffle using wcards.  if not null, just uses
  #             the provided deck.
  #
  #             Dealing rounds > 5 will be treated as 5, results less than 5
  #             result in all rounds "dealt" up to the chosen round.  The 
  #             program uses a simple "redeal everything" logic, so playing 
  #             games a round at a time will be slightly slower than playing
  #             the whole game at once.  drounds < 1 throw errors.
  #
  #             Player numbers are checked against game type to ensure there 
  #             are sufficient cards in the deck to deal all hands.
  #       
  #  Returns:   A "Deck" object with cards assigned to various players (P1-Pn), 
  #             Community Cards (PC) or to "Deck".  Some variants have 
  #             additional categories:
  #               Double Flop Hold-Em has extra community cards (PF)
  #               Good Bad Ugly has "Good" "Bad" "Ugly" & "Discard" categories
  #

  # Create the deck
  if (is.null(deck)) {
    deck <- iwpshuffle(wcards)
  } else { #validate the deck structure
    if (!((sum(names(deck) == names(iwpstddeck)) == dim(iwpstddeck)[2]))) {
      print("Structure of deck parameter")
      print(str(deck))
      print("Structure of the standard deck")
      print(str(iwpstddeck))
      stop("deck parameter does not match expected structure")
    } # end deck structure check 
  } # end of deck logic

  dsize <- dim(deck)[1]

  ############################################################# 
  # Validate parameters
  sgame <- iwpsupportedgames[ngame,]$Stats.Game
  sgame <- ifelse(is.na(sgame),ngame, sgame)
  if (sum(rownames(iwpgames) == sgame) != 1) {
    print(c("Supported Games:"))
    print(sort(rownames(iwpsupportedgames)))
    stop(c("ngame parameter: ", ngame, " is not supported"))
  } else {
    game <- iwpgames[sgame, ] 
  } # end game parameter validation check

  if (is.na(as.integer(players)) | as.integer(players) <= 1) {
    stop(c("players parameter: ", players, " must be number > 1")) 
  } # end players parameter validation check
  players <- as.integer(players)

  mxdeal <- sum(game[c(1, 3, 5, 7, 9)])*players +
            sum(game[c(2, 4, 6, 8, 10, 11)])
  if (dsize < mxdeal) {
    print(paste("Deck Size:  ", dsize, " cards", sep=""))
    print(paste("Max Cards for ", rownames(game), " with ",
            players, " players:  ", mxdeal, sep=""))
    stop("Too many players for chosen game.") 
  } # end players parameter validation check 
 
  if (is.na(as.integer(drounds)) | as.integer(drounds) < 1)  {
    stop(c("drounds (dealing rounds) parameter: ", 
                     drounds, " must be number >= 1"))
  } # end betting rounds parameter validation check

  # Enforce the Fixed Wildcards parameter game$wcfix
  if (game$wcfix != "None") {
    wcval <- unlist(strsplit(game$wcfix, split = " "))
    wcidx <- substr(wcval, nchar(wcval), nchar(wcval)) %in% 
             c("D", "C", "S", "H")
    wcstr <- NULL
    if (sum(wcidx) > 0) {
      wcstr <- c(wcstr, wcval[wcidx])
    } # end wcfix = individual wild card
    if (sum(!wcidx) > 0) {
      wcstr <- c(wcstr, paste(wcval[!wcidx], c("C"), sep=""),
                        paste(wcval[!wcidx], c("H"), sep=""),
                        paste(wcval[!wcidx], c("D"), sep=""),
                        paste(wcval[!wcidx], c("S"), sep=""))
    } # end wcfix = all cards of a single value
    deck[wcstr, 1:2] <- "W" 
  } # end of fixed wildcard check

  ################################################################
  # Initialize card count
  card <- 0 
  # hole cards - sh = starting hole cards, eh = ending hole cards
  sh <- game$shole
  eh <- game$ehole

  ###################################################################  
  # Extra Card rule setup
  excard <- NULL
  if (game$erule != "None") {
    # Elaborate structure here allows support of future similar rules
    if (strsplit(game$erule, " ")[[1]][1] == "Fixed" &
        !is.na(as.numeric(strsplit(game$erule, " ")[[1]][2]))) {
      # Good Bad Ugly logic preprocessor
      if (rownames(game) == "Good Bad Ugly") {
        deck[card + 1, ]$hand <- "1 Good"
        deck[card + 2, ]$hand <- "2 Bad"
        deck[card + 3, ]$hand <- "3 Ugly"
      } # End Good Bad Ugly preprocessor
      card <- card + as.numeric(strsplit(game$erule, " ")[[1]][2])
    } # End of Extra Fixed card assignment
    # Allows "extra card after card type" rule
    if (strsplit(game$erule, " ")[[1]][1] == "After" &
        !is.na(strsplit(game$erule, " ")[[1]][2])) {
      excard <- c(excard, strsplit(game$erule, " ")[[1]][2])
    } # End of Extra After Card assignment
  } # end of erule setup
  ##############################################################
  # Start dealing cards
  for (i in 1:ifelse( drounds > 5, 5, drounds)) {
    # pc = player cards, sc = shared cards
    pc <- as.numeric(game[2*i - 1])
    sc <- as.numeric(game[2*i])
    
    if ( pc != 0 | sc != 0) {
      #################################################################
      # deal player cards
      if (pc != 0) {
        for (k in 1:players) {
          pn <- paste("P", k, sep="")
          deck[(card+1):(card + pc), "hand"] <- pn
          card <- card + pc

          ############################################################
          # Hole Card Logic assumes only start and end hole cards exist
          if (sh > 0 & i == 1) {  # Starting hole cards begin at card#1
            deck[(card-pc+1):(card-pc+sh), ]$holec <- TRUE
          } # end of hole card check
          if (eh > 0 & i == 5) {  # Ending hole cards are the final cards
            deck[(card-eh+1):card, ]$holec <- TRUE
          } # end of hole card check

          ############################################################
          # Deal an extra following card(s) if current card matches value
          # and if the current card is not a hole card
          if (length(excard) > 0 & !deck[card, "holec"]) {
            for (j in 1:length(excard)) { 
              if (namevals(rownames(deck[card, ])) == excard[j]) {
                deck[card+1, "hand"] <- pn
                card <- card + 1
              } # end assign extra card if value matches
            } # end loop through extra cards
          } # end extra card rule logic
        } # end deal to players kiio
      } # end check if player cards are dealt
      
      ################################################################
      # deal shared (community) cards
      if (sc != 0) {
        if (game$mkhand == "2Flop") { # 2 community hands
          deck[(card + 1):(card + sc/2), "hand"] <- "PC"
          deck[(card + 1 + sc/2):(card + sc), "hand"] <- "PF"
        } else { # only one community hand
          deck[(card + 1):(card + sc), "hand"] <- "PC"
        } # end of 2Flop check
        card <- card + sc
      } # end check if shared cards are dealt
    } # end no more cards to deal check   
  } # end betting rounds loop 

  ################################################################
  # Enforce floating wildcard logic
  #   All wildcards that appear in play are handled after all dealing
  #   is done, although if drounds <5, they might not be in final state
    #########################
  if (game$wcscope =="Global") {
    switch(game$wctype,
      ###########
      "Community" = {
        wv <- deck[deck$hand == "PC", ]$value[as.numeric(game$wccard)]
        deck[deck$value == wv, 1:2] <- "W"
      },  # end Global Community logic
      "Player" = {
        wv <- players
        deck[deck$value == wv, 1:2] <- "W"
      },  # end Global Countdown logic
      ##########
      "Sequence" = {  # this will work for sequences wccard in 2:10
        upcards <- substr(deck$hand,1,1) == "P" & !deck$holec
        upnum   <- namevals(rownames(deck[upcards, ]), "Values")
        upcards <- namevals(rownames(deck[upcards, ]))
        first2 <- grep(game$wccard, upcards)[1]
        # No wildcards if minimum number is > 2
        if (!is.na(first2)) {
          seq <- as.numeric(game$wccard)
          if (first2 < length(upnum)) {
            for (i in (1 + first2):length(upnum)) {
              if (as.numeric(upnum[i]) == seq +1 ) {  
                seq <- seq + 1
              } # end sequence difference check
            } # end loop through upcards
          } # end of set sequence number
          deck[paste(valface(seq), c("D","S","H","C"), sep=""), 1:2] <- "W"                
        } # end of wildcard calulation
      },  # end Global Sequence logic
      ########
      "Follow" = {
        uphand <- substr(deck$hand,1,1) == "P" & !deck$holec
        upcards <- namevals(rownames(deck[uphand, ]))
        uploc <- length(upcards)
        mcard <- NULL
        ######
        if (game$wccard == "Pair") {
          for (i in players:1) {
            plup <- deck$hand == paste("P", i, sep="") & !deck$holec
            plup <- namevals(rownames(deck[plup, ]))
            ppair <- FALSE
            for (j in length(plup):2) {
              if (length(grep(plup[j], plup)) > 1) {
                if (players*j < uploc) {
                  mcard <- upcards[players*j + 1]
                } # end last card up check
                ppair <- TRUE
              } # end pair in hand
            } # end loop through player upcards
            if (ppair) {
              break
            } # end check for pair in hand
          } # end loop through players
        } else {
        ######
          if (game$wccard == "Match") {
            for (i in uploc:(players+1)) {
              if (length(grep(upcards[i], upcards)) > 1) {
                if (i + 1 <= uploc) {
                  mcard <- upcards[i + 1]
                } # end last card up check
                break
              } # end match upcard
            } # end loop through upcards
        ######
          } else { # Folow the Card logic
            cardlocs <- grep(game$wccard, upcards)
            if (sum(cardlocs) > 0 ) {
              lastcard <- cardlocs[length(cardlocs)]
              if (lastcard < uploc - 1 ) {
                mcard <- upcards[lastcard + 1]
              } # end card at end test
            } # end follow Card logic
          } # end follow Match test
        } # end follow Pair test
        if (length(mcard) > 0) {
          deck[paste(mcard, c("D","S","H","C"), sep=""), 1:2] <- "W" 
        } # end setting variable
      },  # end Global Follow logic
    ) # end global type switch 
  } else { 
  ################################
    if (game$wcscope =="Hand") { 
      switch(game$wctype,
       ##########
        "Card" = {
          wclist <- NULL
          for (i in 1:players) {
            phand <- rownames(deck[deck$hand == paste("P", i, sep=""), ])
            pwcval <- namevals(phand[as.numeric(game$wccard)])
            wclist <- c(wclist, phand[grep(pwcval, phand)])
          } # end loop through players
          deck[wclist, 1:2] <- "W" 
        }, # end Hand Card logic
        "Pick" = {
          wclist <- NULL
          for (i in 1:players) {
            phand <- rownames(deck[deck$hand == paste("P", i, sep=""), ])
            pwpick <- table(as.numeric(namevals(
                            phand[1:as.numeric(game$wccard)], "Values")))
            pwcval <- valface(min(as.numeric(names(
                            pwpick[pwpick == max(pwpick)]))))
            wclist <- c(wclist, phand[grep(pwcval, phand)])
          } # end loop through players
          deck[wclist, 1:2] <- "W" 
        }, # end Hand Card logic
      ###########
        "Hole" = {
          wclist <- NULL
          for (i in 1:players) {
            phand <- rownames(deck[deck$hand == paste("P", i, sep=""), ])
            phole <- deck[phand, ]
            phole <- rownames(phole[phole$holec, ])
            # this logic finds low card in the hole (wccard = LOW)
            pwild <- sort(as.numeric(namevals(phole, "Values")))[1]
            # Match means that the low card in the hole must
            # match at least one other card in hand or nothing is wild. 
            if (substr(game$wccard,1,3) == "Low" | (game$wccard == "Match" & 
                length(grep(valface(pwild), phand)) > 1)) {
              wclist <- c(wclist, phand[grep(valface(pwild), phand)])
            } # end Match check 
            if (game$wccard == "Low or K" & length(grep("K", phole)) > 0) {
              wclist <- c(wclist, phole[grep("K", phole)])
            } # end low or king in hole check    
          } # end loop through players
          deck[wclist, 1:2] <- "W" 
        }  # end Hand Hole Card logic
      ) # end hand type switch
    } # end Hand wildcard logic
  } # end Global wildcard logic

  ##############################################################################
  # Enforce logic from extra fixed cards
  if (strsplit(game$erule, " ")[[1]][1] == "Fixed") {
    # Good Bad Ugly logic postprocessor
    #   if drounds <4, this might give a skewed result
    if (rownames(game) == "Good Bad Ugly") {
      goodidx <- namevals(rownames(deck)) ==
                 namevals(rownames(deck[deck$hand == "1 Good", ])) 
      badidx  <- namevals(rownames(deck)) == 
                 namevals(rownames(deck[deck$hand == "2 Bad", ]))
      uglypl  <- namevals(rownames(deck)) == 
                 namevals(rownames(deck[deck$hand == "3 Ugly", ]))
      uglypl  <- deck[uglypl & !deck$holec & deck$hand != "3 Ugly"
                 & deck$hand != "Deck" & deck$hand != "Discard", ]$hand

      # Good cards are wild
      deck[goodidx, 1:2]   <- "W"    
      # Bad cards are discarded                 
      deck[badidx & deck$hand != "2 Bad", "hand"]  <- "Discard"
      # Players lose if they match ugly with a card up
      if (length(uglypl) > 0) {
        for (i in 1:length(uglypl)) {
          deck[deck$hand == uglypl[i], "hand"] <- paste(uglypl[i],"Ugly")  
        } # end loop through ugly
      } # end ugly player check
    } # end Good Bad Ugly postprocessor
  } # end fixed erule check
                
  deck
} # end internal function 1wpdeal

################################################################################
iwptable <- function(deck, cross = FALSE) {
  #  This function creates human readable strings showing the "hand" of each
  #  player, excluding the "Deck" hand, providing a visualization of what could
  #  be seen on the table if all hole cards were turned up.
  #
  #  This function would not normally be  called unless the state of the deck 
  #  at a particular point was of interest to humans analyzing the results.  
  #
  #  Args:
  #    deck:  A deck data frame of the sort generated by the "iwpdeck" function 
  #    cross: Forms a cross out of 5 community cards, if 5 exist
  #
  #  Returns:  List with the following elements
  #    $Extra:   Currently only used in "Good Bad Ugly"
  #              eg:  if Deuces are wild, "2C" would appear as "2Cw"
  #              eg:  ace of spades in the hole shows as "#AS"
  #    $Community:  Community cards with deck$hand == "PC" or "PF"
  #              if cross is TRUE, PC is a 3x3 matrix
  #    $pcards:  All deck$hand entries with form "P1", "P2"... "PN"
  #
  #    Displays one value per unique "hand" in the deck, with card names 
  #    separated by a space.  Cards that are wild have a small "w" 
  #    concatenated and hole cards are prefixed with a hash (#)
  #   
  if (!((sum(names(deck) == names(iwpstddeck))  == dim(iwpstddeck)[2]))) {
    print("Structure of deck parameter")
    print(str(deck))
    print("Structure of the standard deck")
    print(str(iwpstddeck))
    stop("deck parameter does not match expected structure")
  } # end deck structure check  
  #####
  # main processing
  gtable <- NULL
  for (i in sort(unique(deck$hand))) {
    if (i != "Deck") {
      gtable[i] <- paste(ifelse(deck[deck$hand == i,2] == "W", 
                     paste(rownames(deck[deck$hand == i, ]),"w",sep = ""),
                       rownames(deck[deck$hand == i, ])), collapse = " ")  
      gtable[i] <- paste(ifelse(deck[deck$hand == i, "holec"], 
                     paste("#", unlist(strsplit(gtable[i], " ")), sep = ""), 
                       unlist(strsplit(gtable[i], " "))), collapse=" ")                  
    }  # end check to see if it is the remaining deck
  } # end loop through unique hands
  ######
  # Generate the visual table coefficients 
  pcards <- gtable[substr(names(gtable), 1, 1) == "P" &
                  !(names(gtable) %in% c("PC" ,"PF")) ]
  ## Community Cards
  if (sum(names(gtable) %in% c("PC" ,"PF")) > 0) {
    ccards <- gtable[names(gtable) %in% c("PC" ,"PF")]
    ## Cross logic
    if (cross & (length(deck[deck$hand == "PC", "hand" ]) == 5)) {
      ccards <- unlist(strsplit(ccards, split = " "))
      ccards <- matrix(c("...", ccards[1], "...",
                     ccards[4], ccards[5], ccards[2],
                     "...", ccards[3], "..."), ncol=3, byrow = TRUE)
    } # end of cross logic
  } else {
    ccards <- "None"
  } # end of ccards logic
 
  ## Extra Cards
  if (sum(substr(names(gtable), 1, 1) != "P") > 0) {
    ecards <- gtable[substr(names(gtable), 1, 1) != "P"]
  } else {
    ecards <- "None"
  } # end of ecards logic

  list(Extra = ecards,  
       Community = ccards,
       Players = pcards)
} # End of internal function iwptable

################################################################################

iwpshowdown <- function(deck, ngame) {
  #  This function simulates the "showdown" portion of the game, where
  #  each player reveals their best hand.  Unlike iwptable, community cards
  #  are incorporated into each players hand.  
  #
  #
  #  Args:
  #    deck:  A deck data frame of the sort generated by the "iwpdeck" function
  #    ngame:  a legal game from "rownames(wpsupportedgames)"  
  #
  #  Returns: showdown data frame with one row per unique player in the deck
  #           rowname     = player name (P1, P2 etc)
  #          $mainraw     = text of hand in wptable format
  #          $mainhand    = text of hand used in wpmakehand format
  #          $maintype    = type of hand (eg 3-Kind or FullHouse)
  #          $mainscore   = raw score, used to break ties
  #          $splitraw    = text of hand in wptable format
  #          $splithand   = text of hand used in wpmakehand format
  #          $splittype   = type of hand (eg 3-Kind or FullHouse)
  #          $splitscore  = raw score, used to break ties
  #          $potpct      = percent of pot won by an individual player
  #
  #

  ######################
  # Validate parameters
  sgame <- iwpsupportedgames[ngame,]$Stats.Game
  sgame <- ifelse(is.na(sgame),ngame, sgame)
  if (sum(rownames(iwpgames) == sgame) != 1) {
    print(c("Supported Games:"))
    print(sort(rownames(iwpsupportedgames)))
    stop(c("ngame parameter: ", ngame, " is not supported"))
  } else {
    game <- iwpgames[sgame, ] 
  } # end game parameter validation check
  if (!((sum(names(deck) == names(iwpstddeck)) == dim(iwpstddeck)[2]))) {
    print("Structure of deck parameter")
    print(str(deck))
    print("Structure of the standard deck")
    print(str(iwpstddeck))
    stop("deck parameter does not match expected structure")
  } # end deck structure check  

  #######
  # Initialize the key variables
  pnam <- sort(unique(deck[substr(deck$hand,1,1) == "P", ]$hand))
  pnam <- pnam[pnam != "PC" & pnam != "PF"]
  mhl <- ifelse(game$hand == "LO", -1, 1) 
  shl <- ifelse(game$split == "LO", -1, 1)

  #####
  # Generate visual table coefficients for community hands and
  # the matrix of possible hands based on mkhand

  pcdlen <- sum(deck$hand == "P1")
  ccdlen <- sum(deck$hand == "PC")
  if (game$mkhand == "Cross") {
    phidx  <- matrix(c(1:pcdlen, pcdlen + 4, pcdlen + 5, pcdlen + 2,
                     1:pcdlen, pcdlen + 1, pcdlen + 5, pcdlen + 3),
                     nrow = 2, byrow = TRUE)
  } else {
    if (substr(game$mkhand,1,2) == "2H") {
      h2idx <- c(1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4)
      if (pcdlen == 5) {
        h2idx <- c(h2idx, 1, 5, 2, 5, 3, 5, 4, 5)
      } # end hole card=5 check
      h2idx <- matrix(h2idx, ncol = 2, byrow = TRUE)
      if (game$mkhand == "2H3C") {
        cnidx <- c(1, 2, 3, 1, 2, 4, 1, 2, 5,
                   2, 3, 4, 2, 3, 5, 2, 4, 5, 3, 4, 5)
        cnidx <- cnidx + pcdlen
        cnidx <- matrix(cnidx, ncol=3, byrow = TRUE)
      } else {
        cnidx <- matrix(c((pcdlen+1):(pcdlen+ccdlen)),
                        nrow = 1, byrow = TRUE)
      } # end hole card=5 check
      phidx <- matrix(replicate((dim(h2idx)[1] * dim(cnidx)[1]) *
                                (dim(h2idx)[2] + dim(cnidx)[2]), 0),
               ncol = (dim(h2idx)[2] + dim(cnidx)[2]), byrow = TRUE)
      i <- 0
      for (j in 1:dim(h2idx)[1]) {
        for (k in 1:dim(cnidx)[1]) {
          i <- i+ 1
          phidx[i, ] <- c(h2idx[j, ],cnidx[k, ]) 
        } # end cnloop
      } # end h2loop
    } else {
      phidx <- matrix(1:(pcdlen + ccdlen),
                      nrow = 1, byrow = TRUE)
    } # end of 2H3C/2H5C check
  } # end of cross format check

  ######
  # Generate showdown matrix which tracks player outcomes
  sdown <- data.frame( maintype    = replicate(length(pnam), ""),
                       mainhand    = replicate(length(pnam), ""),
                       mainscore   = replicate(length(pnam), "0"),
                       splittype   = replicate(length(pnam), ""),
                       splithand   = replicate(length(pnam), "0"),
                       splitscore  = replicate(length(pnam), 0),
                       winmain     = replicate(length(pnam), 0),
                       winsplit    = replicate(length(pnam), 0),
                       potpct      = replicate(length(pnam), 0),
                      stringsAsFactors = FALSE, row.names = pnam)

  for (i in 1:length(pnam)) {
    ##############################
    # Build the player hands   
    pmhand <- deck[deck$hand %in% c(pnam[i], "PC"), ]
    if (sgame == "Good Bad Ugly" ) {
      if (length(grep("Ugly", pnam[i])) > 0) {
        sdown[i, c("maintype", "mainhand")] <- c("NoHand", "Folded")
        sdown[i, c("splittype", "splithand")] <- c("None", "")
        next
      } # end ugly logic
    } # end Good Bad Ugly game logic
    if (sgame %in% c("Baseball", "Good Bad Ugly") ) {
      phidx <- matrix(1:length(pmhand$value), nrow = 1, byrow = TRUE)
    } # end variable hand size logic
    pshand <- pmhand
    if (game$mkhand == "2Flop") {
        pshand <- deck[deck$hand %in% c(pnam[i], "PF"), ]
    } # end of 2nd community hand check  
    ##############################
    # Build the players showdown hands
    sdown[i, ] <- c("NoHand", "", "0", "None", "", "0", 0, 0, 0)
    for (j in 1:dim(phidx)[1]) {
## Debug code ###
# assign("debugdeck", deck, envir = .GlobalEnv)
# assign("debughand", pmhand[phidx[j, ], ], envir = .GlobalEnv)
#################
      mainhand <- iwpmakehand(pmhand[phidx[j, ], ], 
                             sgame, split = FALSE)
      if (((as.numeric(mainhand[1]) * mhl) >= 
           (as.numeric(sdown[i,"mainscore"]) * mhl)) |
          (mainhand[1] != "0" & sdown[i,"mainscore"] == 0)) {
        sdown[i, c("mainscore", "mainhand", "maintype")] <- mainhand       
      } # end set main showdown hand
      if (game$split != "None") {
        splithand <- iwpmakehand(pshand[phidx[j, ], ], 
                                 sgame, split = TRUE)
        if (((as.numeric(splithand[1]) * shl) >= 
             (as.numeric(sdown[i,"splitscore"]) * shl)) |
            (splithand[1] != "0" & sdown[i,"splitscore"] == 0)) {
          sdown[i, c("splitscore", "splithand", "splittype")] <- splithand  
        } # end set split showdown hand
      } # end test for if split hand exists
    } # end loop through legal hand combinations 
  } # end loop through players
  ###################
  # Redeal logic
  if (game$redeal == "Both") {
    bidx <- sdown$mainscore == max(as.numeric(sdown$mainscore)) &
            sdown$splitscore == max(as.numeric(sdown$splitscore)) &
            sdown$mainscore > 0 & sdown$splitscore > 0
    if (sum(bidx) == 0) {
      sdown$maintype <- "Redeal"
      sdown$splittype <- "Redeal"
    } # end no winner check
    sdown$mainscore[!bidx] <- 0
    sdown$splitscore[!bidx] <- 0
  } # end "Both must win or redeal" logic
  if (game$redeal == "QS Up" & (sum(substr(deck$hand,1,1) == "P" 
       & deck$holec == FALSE & rownames(deck) == "QS") > 0)) {
    sdown$maintype <- "Redeal"
    sdown$mainscore <- 0
  } # end redeal if QS is up logic
  if (game$redeal == "No Follow Up" & 
      (sum(substr(deck$hand,1,1) == "P" 
           & deck$holec == FALSE 
           & substr(rownames(deck),1,1) == game$wccard
       ) == 0)) {
    sdown$maintype <- "Redeal"
    sdown$mainscore <- 0
  } # end redeal if no Queens are up logic
  ####################
  # Score the winner(s)
  sdown$mainscore <- as.numeric(sdown$mainscore)
  sdown$splitscore <- as.numeric(sdown$splitscore)
  sdown[, c("winmain", "winsplit", "potpct")] <- as.numeric(0)
  if (sum(!(sdown$maintype %in% c("None", "NoHand", "Redeal"))) > 0) { 
    mwidx <- sdown$mainscore
    mwidx <- mwidx[!(sdown$maintype %in% c("None", "NoHand"))]
    mwin <- ifelse(mhl == 1, max(mwidx), min(mwidx))
    sdown[sdown$mainscore == mwin, ]$winmain <- as.numeric(1)
  } # end of mainwin logic
  if (sum(!(sdown$splittype %in% c("None", "NoHand", "Redeal"))) > 0) { 
    swidx <- sdown$splitscore
    swidx <- swidx[!(sdown$splittype %in% c("None", "NoHand"))]
    swin <- ifelse(shl == 1, max(swidx), min(swidx))
    sdown[sdown$splitscore == swin, ]$winsplit <- as.numeric(1)
  } # end of splitwin logic
  mwinners <- sum(sdown$winmain)
  swinners <- sum(sdown$winsplit)
  if ((mwinners > 0) & (swinners > 0)) {
    sdown$potpct <- (sdown$winmain*.5/mwinners) + 
                    (sdown$winsplit*.5/swinners)
  } else { # spot was not split
    if (mwinners > 0) {
      sdown$potpct <- (sdown$winmain/mwinners)
    } else {
      if  (swinners > 0) {
        sdown$potpct <- (sdown$winsplit/mwinners)
      } # end nobody wins logic
    } # end unsplit pot logic
  } # end pot percentage logic
  sdown
} # end internal function iwpshowdown
################################################################################
