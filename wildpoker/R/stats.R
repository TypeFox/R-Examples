##################################################################################
wprules <- function(ngame) {
  #  This function "plays" an entire poker game, from shuffling the deck
  #  to the showdown. It collects extra statistics on the game itself and
  #  stores everything in a single $summary" vector for use by wpstats
  #
  #  Args:
  #    ngame:   any value from rownames(wpsupportedgames) 
  #
  #  Returns: prints the following output to the display device 
  #    (wsg, below = wpsupportedgames)
  #    GAME :  ngame (wsg$Game.Type, max players = wsg$Max.Players)}
  #    DEAL :  wsg$Deal
  #    BETS :  wsg$Bet
  #    WILD :  wsg$Wildcards
  #    MAIN :  wsg$Main.Hand
  #    SPLIT:  wsg$Split.Hand
  #    NOTES:  wsg$Notes
  #
  wsg <- iwpsupportedgames[ngame, ]
  if (rownames(wsg) == "NA") {
    print(c("Supported Games:"))
    print(sort(rownames(iwpsupportedgames)))
    stop(c("ngame parameter: ", ngame, " is not supported"))
  } # end game parameter validation check
  rv1 <- paste("GAME : ", ngame, " (", wsg$Game.Type, ", max players = ",
              wsg$Max.Players, ")", " \n", sep="")
  rv2 <- paste("DEAL : ", wsg$Deal, " \n", sep="")
  rv3 <- paste("BETS : ", wsg$Bet,  " \n", sep="")
  rv4 <- paste("WILD : ", wsg$Wildcards, " \n", sep="")
  rv5 <- paste("MAIN : ", wsg$Main.Hand, " \n", sep="")
  rv6 <- paste("SPLIT: ", wsg$Split.Hand, " \n", sep="")
  rv7 <- paste("NOTES: ", wsg$Notes, " \n", sep="")
  cat(rv1, "\n", rv2, "\n", rv3, "\n", rv4, "\n", 
      rv5, "\n", rv6, "\n", rv7, sep="")
} # end of external function wprules

##################################################################################
wpgame <- function(ngame, players, wcards = NULL) {
  #  This function "plays" an entire poker game, from shuffling the deck
  #  to the showdown. It collects extra statistics on the game itself and
  #  stores everything in a single $summary" vector for use by wpstats
  #
  #  Args:
  #    ngame:   any value from rownames(wpsupportedgames) 
  #    players: Number greater than zero, less than max players for a game
  #    wcard:   a vector with string values either matching the 52 cards
  #             in a normal deck, or several common designations.  Supported
  #             wildcard designations are:
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
  #
  #  Returns: a "pgame" object (list) with the following coefficients:
  #        $ntable    = list of cards showing the deal visually
  #          wildcards are tagged with "w" (eg, "2Dw")
  #          hole cards are prefixed with "#" (eg, "#2Dw"
  #          $Extra     = Extra cards (currently only in "Good Bad Ugly"
  #          $Community = community cards (deck$hand = PC or PF)
  #          $Players   = player cards (deck$hand of form "P1"..."Pn")
  #        $showdown = list showing how each player used their hand
  #  Returns: showdown data frame with one row per unique player in the deck
  #           rowname     = player name (P1, P2 etc)
  #          $maintype    = type of hand (eg 3-Kind or FullHouse)
  #          $mainhand    = text of hand used in $ntable format
  #          $mainscore   = raw score, used to break ties
  #          $splittype   = type of hand (eg 3-Kind or FullHouse)
  #          $splithand   = text of hand used in $ntable format
  #          $splitscore  = raw score, used to break ties
  #          $potpct      = percent of pot won by an individual player
  #        $summary  = character vector with statistics for game as whole
  #          $game     = Name of poker game
  #          $wild     = Additional Wildcards
  #          $mwtype   = Winning hand type of main hand
  #          $swtype   = Winning hand type of split hand
  #        $detail   = character vector with statistics for game as whole
  #          $mwscore  = score of winning main hand
  #          $msscore  = # of players with score = winning main hand
  #          $mstype   = # of players with type = winning main hand
  #          $swscore  = score of winning split hand
  #          $ssscore  = # of players with wmscore = winning split hand
  #          $sstype   = # of players with wmtype = winning split hand
  #          $wcdeck   = Number of Wildcards in deck before deal
  #          $wcdeal   = Number of Wildcards in deck after deal
  #          $pnum     = players
  #          $mppct    = % of main pot earned if main hand won (0-100)
  #          $sppct    = % of main pot earned if split hand won (0-100)
  #          $bppct    = % of entire pot earned if both hands won (0-100)
  #
  
  #### 
  # shuffle and deal the cards. Wildcards not matching legal values are 
  # ignored by iwpshuffle
  deck <- iwpshuffle(wcards)
  wcn <- sum(deck$suit == "W")
  #  the iwpdeal program validates all of the parameters, so just pass
  #  them directly. 
  deck <- iwpdeal(ngame, players, deck = deck)
  #  add the count of fixed wildcards based on game variant to wcn
  if (iwpgames[ngame, ]$wcfix != "None") {
    wcval <- unlist(strsplit(iwpgames[ngame, ]$wcfix, split = " "))
    wcidx <- substr(wcval, nchar(wcval), nchar(wcval)) %in% 
             c("D", "C", "S", "H")    
    wcn <- wcn + sum(wcidx) + 4 * sum(!wcidx) 
  } # end add in fixed wild cards
  # NoWin is for cases where no hand even qualifies to win, and betting
  # is therefore assumed to be nonexistant - only the ante remains in the
  # pot.  Redeal is for cases where hands strong enough to bet on exist,
  # but game variant restrictions can cause everybody to lose with a large
  # pot moving forward to the next game.  Redeal is considered to be both
  # the largest possible hand and a tie with all players, which graphs
  # differently, as such variants greatly increase risk during betting.
  if (iwpgames[ngame, ]$redeal != "None") {
    nowin <- "Redeal"
  } else {
    nowin <- "NoWin"
  } # end of NoWin vs Redeal logic
  cross <- iwpgames[ngame, ]$mkhand == "Cross"
  pgame <- list(table = iwptable(deck, cross),
                showdown = iwpshowdown(deck, ngame),
                summary  = c(game = ngame, 
                             wild = paste(wcards, collapse = ", "),
                             winmainhand = "", winsplithand = ""),
                detail   = c(mwscore = 0, msscore = 0, mstype = 0,
                             swscore = 0, ssscore = 0, sstype = 0, 
                             wcdeck = wcn, wcdeal = 0, pnum = players,
                             mppct = 0, sppct = 0, bppct = 0))
  ######
  # Set wcdeck (the # of wildcards in the final, after-deal deck.
  #  only include cards that don't bucket as "extra"
  pgame[["detail"]]["wcdeal"] <- sum(deck$suit == "W" &
                                     (substr(deck$hand,1,1) == "P" | 
                                      deck$hand == "Deck"))
  #######
  # Set main hand values
  midx <- pgame[["showdown"]]["winmain"] == 1
  if (sum(midx) > 0) {
    pgame[["detail"]]["mppct"] <- 100 * round((1/sum(midx)),2)
    tval <- pgame[["showdown"]]["maintype"][midx]
    pgame[["summary"]]["winmainhand"] <- tval[1]
    sval <- pgame[["showdown"]]["mainscore"][midx]
    pgame[["detail"]]["mwscore"] <- sval[1]
    pgame[["detail"]]["msscore"] <- length(sval)
    sidx <- pgame[["showdown"]]["maintype"] == tval
    pgame[["detail"]]["mstype"] <- sum(sidx)
  } else { # nobody won main hand
    pgame[["summary"]]["winmainhand"] <- nowin
  }
  #######
  # Set split hand values
  if (iwpgames[ngame, ]$split != "None")  {
    tidx <- pgame[["showdown"]]["winsplit"] == 1
    if (sum(tidx) > 0) {
      pgame[["detail"]]["sppct"] <- 100 * round((1/sum(tidx)),2)
      tval <- pgame[["showdown"]]["splittype"][tidx]
      pgame[["summary"]]["winsplithand"] <- tval[1]
      sval <- pgame[["showdown"]]["splitscore"][tidx]
      pgame[["detail"]]["swscore"] <- sval[1]
      pgame[["detail"]]["ssscore"] <- length(sval)
      sidx <- pgame[["showdown"]]["splittype"] == tval
      pgame[["detail"]]["sstype"] <- sum(sidx)
      if (pgame[["summary"]]["winmainhand"] == nowin) {
        pgame[["summary"]]["winmainhand"] <- "SplitWon"
        pgame[["detail"]]["msscore"] <- pgame[["detail"]]["ssscore"]
        pgame[["detail"]]["mstype"]  <- pgame[["detail"]]["ssscore"]  
      } # end split hand gets entire pot
    } else { # no split hand won, see if main hand won
      if (pgame[["summary"]]["winmainhand"] != nowin) {
        pgame[["summary"]]["winsplithand"] <- "MainWon"
        pgame[["detail"]]["ssscore"] <- pgame[["detail"]]["msscore"]
        pgame[["detail"]]["sstype"]  <- pgame[["detail"]]["msscore"] 
      } else { # else nobody won, main or split 
        pgame[["summary"]]["winsplithand"] <- nowin
      } # end pot goes to main hand check
    } # end nobody won split hand
   } else {  # no split hand, set value to None
     pgame[["summary"]]["winsplithand"] <- "None"
  } # end of split hand logic
  bidx <- pgame[["showdown"]]["winmain"] == 1 &
          pgame[["showdown"]]["winsplit"] == 1
  if (sum(bidx) > 0) {
    pgame[["detail"]]["bppct"] <- round(100 * 
                                  pgame[["showdown"]][bidx, "potpct"][1])
  } # end both pot calculation
  pgame
} # end of external function wpgame 


################################################################################## 



iwphandstats <- function(rawstat, split = FALSE) {
  # Utility function to create labels for reporting grouping
  #
  # Args:
  #         rawstat - rawstat object from wpgstats
  #         ngame - a legal game from "rownames(wpsupportedgames)" 
  #         split - if FALSE, use game$hand variables, 
  #                 if TRUE use game$split variables 
  #
  # Returns: a list of tables with 3 values
  #         pctwin - table showing % of hands won by a given hand
  #         cntwin - table showing count of winning hands
  #         wcwin  - table showing wildcard count for each winning hand
  #

  ######################
  # load the game type variables, based on rawstat and split parameters
  players <- rawstat$players[1]
  game    <- iwpgames[rawstat$game[1], ]
  ngame   <- rownames(game)
  if (split) { # make the split hand
    htyp  <- game$split
    hmod  <- game$smod
    hmin  <- game$smin 
    whand <- rawstat$swtype
    wpwin <- rawstat$ssscore
    wptyp <- rawstat$sstype
  } else { # make the primary hand
    htyp  <- game$hand
    hmod  <- game$hmod
    hmin  <- game$hmin 
    whand <- rawstat$mwtype
    wpwin <- rawstat$msscore
    wptyp <- rawstat$mstype
  } # end of game variable assignment
  # Community game logic is not part of this function
  hmod <- ifelse(hmod == "PC" | hmod == "PF", "None", hmod)
   
  if (hmod == "Sum") {
    hlabel <- c("9-", "10+", "20+", "30+", "40+", "50+", 
                "60+", "70+", "80+", "90+")
    if (hmin == "Spots") {
      # There is no such thing as a failure to win in spots
      # and the best hand is 70 points (7 cards 10 value each)
       hlabel <- hlabel[!(hlabel %in% c("70+", "80+", "90+"))]
       hlabel <- paste(hlabel, hmin, sep="")
    } else {  # it is possible to have no winning hand, always split hand
      if (hmin == "None" ) {
        hmin <- ""
      } else {
        hlabel <- paste(hlabel, hmin, sep="")
        hlabel <- c("MainWon", hlabel)
      } # can't lose check.
    } # end of spots check
  } else { # uses more normal hand order
    hlabel <- c( "7-High", "8-High", "9-High", "10-High",
                 "J-High", "Q-High", "K-High", "A-High",
                 "Pair", "2-Pair", "3-Kind", "Straight", "Flush", 
                 "FullHouse", "4-Kind", "StrFlush", "5-Kind")
    if (hmod == "Hole") {
       hlabel <- c("2-High", "3-High", "4-High", "5-High", "6-High", hlabel)
       hrpl <- ifelse(hmin == "None", "-inHole", paste(hmin, "-inHole", sep=""))
       hlabel <- gsub("-High", hrpl, hlabel[1:13])
       if (htyp == "LO") {
         hlabel <- c(rev(hlabel[1:12]),hlabel[13])
       } # end of low logic
       hlabel <- c("MainWon", hlabel)
    } else { # normal logic
      if (htyp == "HI") {
        if (sum(game[1:10]) < 5) {
          hlabel <- hlabel[!(hlabel %in%  c("StrFlush", "Straight", "Flush", 
                                            "FullHouse", "5-Kind"))]
          hlabel <- c("5-High", "6-High", hlabel)
          if (sum(game[1:10]) < 4) {
            hlabel <- hlabel[!(hlabel %in%  c("2-Pair", "4-Kind"))]
            hlabel <- c("4-High", hlabel)
            if (sum(game[1:10]) < 3) {
              hlabel <- hlabel[hlabel != "3-Kind"]
              hlabel <- c("3-High", hlabel)
              if (sum(game[1:10]) < 2) {
                hlabel <- hlabel[hlabel != "Pair"]
                hlabel <- c("2-High", hlabel)
              } # end strip 2 card hands
            } # end strip 3 card hands
          } # end strip 4 card hands
        } else { #  normal 5 card poker hands
          if (sum(rawstat$wcdeal) == 0) {
            hlabel <- hlabel[hlabel != "5-Kind"]
          } # end of 5-kind check
          if (hmin != "None") {
            hlabel <- hlabel[grep(hmin, hlabel)[1]:length(hlabel)]
            if (split) {
              hlabel <- c("MainWon", hlabel)
            } else {
              hlabel <- c("SplitWon", hlabel)
            } # end split check
            if (game$smin != "None" & game$hmin != "None") {
              hlabel <- c("NoWin", hlabel)
            } # end check for nobody wins 
          } else { # convert anything lower than a pair to HighCard
            hlabel <- hlabel[grep("Pair", hlabel)[1]:length(hlabel)]
            hlabel <- c("HighCard", hlabel)
            if (sum(rawstat$mwtype == "NoWin") > 0 ) {
                  hlabel <- c("NoWin", hlabel)
            } # end NoWin check
            whand[!(whand %in% c("Redeal", hlabel))] <- "HighCard"
          } # end of hmin check         
        } # end check for less than 5 cards in hand
      } else { # low hand logic
        hlabel <- rev(hlabel)
        if (hmod == ".A-6") {
          hlabel <- c(hlabel, "6-High")
        } # end of A-6 check
        if (hmod == ".A-5") {
          hlabel <- c(hlabel, "6-High", "5-High")
          hlabel <- hlabel[!(hlabel %in%  c("StrFlush", "Straight", "Flush"))]                                     
        } # end of A-5 check
        if (hmin != "None") {
          hml <- paste(hmin, "High", sep="")
          hlabel <- hlabel[grep(hml, hlabel)[1]:length(hlabel)]
          if (split) {
            hlabel <- c("MainWon", hlabel)
          } else {
            hlabel <- c("SplitWon", hlabel)
          } # end split check
          if (game$smin != "None" & game$hmin != "None") {
            hlabel <- c("NoWin", hlabel)
          } # end check for nobody wins 
        } else { # normal low hand logic
          hlabel <- hlabel[grep("A-High", hlabel)[1]:length(hlabel)]
          hlabel <- c("Pair+", hlabel)
          whand[!(whand %in% hlabel)] <- "Pair+"
        } # end of hmin check            
      } # end of HI vs LO check
    }  # end of hole logic check
  }  # end of sum logic check
  if (game$redeal != "None") {
    hlabel <- c(hlabel, "Redeal")
  } # end redeal logic
  wptyp[wpwin > 1] <- 0
  wtable <- data.frame(matrix(replicate(length(hlabel) * (players + 1), 0),
                              byrow = TRUE, nrow = (players + 1)))
  rownames(wtable) <- c(1:players, "tie")
  names(wtable) <- hlabel                      
  ## add data to the table                          
  rtable <- table(whand, wptyp)
  for (i in  rownames(rtable)) {
    for (j in as.character(unique(wptyp))) {
      if (j == "0") {
        wtable["tie", i] <- rtable[i, j]
      } else {
        wtable[j, i] <- rtable[i, j]
      } # end tie check
    } # end of name loop
  } # end of rowname loop

  ###############
  # create ptable
  ptable <- wtable
  csum <- 0
  for (i in  names(ptable)) {
    for (j in rownames(ptable)) {
      if (j == "1" & csum > 0 ) {
        ptable[j, i] <- csum + ptable[j, i]
      } # end of zero check
    } # end of name loop
    csum <- sum(ptable[, i])
  } # end of rowname loop
  ptable <- ptable/sum(wtable)
  ##################
  # create wctable
  wcnum <- length(names(table(rawstat$wcdeal)))
  wctable <- data.frame(matrix(replicate(length(hlabel) * wcnum, 0),
                              byrow = TRUE, nrow = wcnum))
  rownames(wctable) <- names(table(rawstat$wcdeal))
  names(wctable) <- hlabel 
  rtable <- table(rawstat$wcdeal, whand)   
  for (i in unique(whand)) {
    for (j in rownames(wctable)) {
      wctable[as.character(j), i] <- rtable[as.character(j), i]
    } # end of name loop
  } # end of rowname loop

#####
#debug code 
# print(c(rownames(game), players, rawstat$wcdeck[1]))

  list(pctwin <- ptable, cntwin <- wtable, wcwin <- wctable) 
} # end of internal function handstats


################################################################################## 

wpstats <- function(ngame, players, wcards = NULL, numdeal = 1000, seed = 52, 
                     raw = FALSE) {
  #    This function is designed to run many games, and capture the summary
  #    data from the game object for later analysis
  #
  #  Args:
  #    ngame:    any rowname from wpsupportedgames 
  #    players:  Number greater than zero, less than max players for a game
  #    wcards:   If additional wild cards are desired, see ?wpgame for
  #              how to populate this field
  #    numdeal:  number of games to deal (size of output data frame)
  #    seed:     random number seed
  #    raw:      if TRUE, returns a rawstat data frame
  #              if FALSE, returns a gstat list
  #
  #  Returns: one of the following 
  #    rawstat data frame, rows = numdeal
  #          $game     = Name of poker game
  #          $pnum     = players
  #          $wild     = Additional Wildcards
  #          $mwtype   = Winning hand type of main hand
  #          $swtype   = Winning hand type of split hand
  #          $mwscore  = score of winning main hand
  #          $msscore  = # of players with score = winning main hand
  #          $mstype   = # of players with type = winning main hand
  #          $swscore  = score of winning split hand
  #          $ssscore  = # of players with wmscore = winning split hand
  #          $sstype   = # of players with wmtype = winning split hand
  #          $mppct    = % won if main hand is won
  #          $sppct    = % won if split hand is won
  #          $bppct    = % won if both hands are won
  #          $wcdeck   = Number of Wildcards in deck before deal
  #          $wcdeal   = Number of Wildcards in deck after deal
  #    or
  #
  #    gstat list:
  #          $game     = c(pname, players, wcards)
  #          $stats    = c(numdeal, seed, wcdeck)
  #          $pmain    = p values for winning hands by player with hand
  #          $cmain    = count of winning hands  by player with hand
  #          $wmain    = count of wildcards in deck for winning hand
  #          $psplit   = p values for winning hands by player with hand
  #          $csplit   = count of winning hands by player with hand
  #          $wsplit   = count of wildcards in deck for winning hand
  #          $potpct   = percent of pot won by winning players 
  #         
  #
  #
  etime <- proc.time()
  #####
  # validate numdeal
  if (is.na(as.integer(numdeal)) | as.integer(numdeal) < 1) {
    print(c("parameter numdeal is:", numdeal))
    stop ("numdeal must be positive integer")
  } else {
    numdeal <- as.integer(numdeal)
  } # end numdeal validation
  # set.seed does its own validation
  set.seed(seed)


 # initialize the data frames
  rawstat <- data.frame(game = "", players = 0, wild = "", 
                        mwtype = "", swtype = "",
                        mwscore = 0, msscore = 0, mstype = 0,
                        swscore = 0, ssscore = 0, sstype = 0,
                        mppct = 0, sppct = 0, bppct = 0,
                        wcdeck = 0, wcdeal = 0, stringsAsFactors = FALSE)
  gstat   <- list(game  = c(ngame = ngame, players = players,
                            wcards = paste(wcards, collapse = " ")),
                  stats = c(numdeal = numdeal, seed = seed, wcdeck = 0))
                  
  # wpgame validates the other parameters
  for (i in 1:numdeal) {
    pgame <- wpgame(ngame, players, wcards)
    rawstat[i, c("game", "wild", "mwtype", "swtype")] <- pgame[["summary"]]
    rawstat[i, c("mwscore", "msscore", "mstype",
                 "swscore", "ssscore", "sstype", 
                 "wcdeck", "wcdeal", "players",
                 "mppct", "sppct", "bppct" )] <- pgame[["detail"]]
  } # end of numdeal loop
## Debug code ###
# assign("debugraw", rawstat, envir = .GlobalEnv)
#################
  if (raw) {
    rawstat
  } else {
    gstat[["stats"]]["wcdeck"] <- rawstat$wcdeck[1]
    hstat <- iwphandstats(rawstat, split = FALSE)
    gstat$pmain <- hstat[1]
    gstat$cmain <- hstat[2]
    gstat$wmain <- hstat[3]
    if (iwpgames[rawstat[1, "game"], "split"] == "None") {
      hstat <- c("", "", "", "")
      potpct <- table(rawstat$mppct)/numdeal 
    } else { # load the split data & work out pot percent
      hstat <- iwphandstats(rawstat, split = TRUE)
      potval <- rawstat[rawstat$mppct == 0 & 
                        rawstat$sppct == 0, "mppct"]
      rawstat[rawstat$bppct > 0, c("mppct","sppct")] <- 0
      potval <- c(potval, rawstat[rawstat$bppct > 0,"bppct"])
      potval <- c(potval, rawstat[rawstat$mppct > 0 & 
                                  rawstat$sppct == 0, "mppct"])
      potval <- c(potval, rawstat[rawstat$sppct > 0 & 
                                  rawstat$mppct == 0, "sppct"])
      potval <- c(potval, potval)
      potval <- c(potval, rawstat[rawstat$sppct > 0 & 
                                  rawstat$mppct > 0, "mppct"]/2)
      potval <- c(potval, rawstat[rawstat$sppct > 0 & 
                                  rawstat$mppct > 0, "sppct"]/2)
      potpct <- table(potval)/(2*numdeal) 
    } # end of split hand logic
    gstat$psplit <- hstat[1]
    gstat$csplit <- hstat[2]
    gstat$wsplit <- hstat[3]
    gstat$potpct <- potpct
    gstat
  } # end of rawstat vs gstat check
} # end of external function wpgstats
################################################################################## 

iwpgcreatedb <- function(pokerdb = NULL, gamelist = rownames(iwpgames),
                         numdeal = 1000, seed=52 ) {
  #    This app is designed to loop throught the supported games and run every 
  #    combination of game, number of players (from 2-8, up to max players/game)
  #    & wildcards(0-8, combinations of Suicide King, One-Eyed Jack and Deuces)
  #
  # Args:
  #    pokerdb:  a pokerdb list from this function - anything in the db is skipped
  #    gamelist: must be a subset of rownames(wpsupportedgames) - others skipped
  #    numdeal:  number of games to deal (size of output data frame)
  #    seed:     random number seed
  # 
  # Returns: pokerdb list, which is a list of gstat lists, with list names =
  #          "name.of.game".# of players.# of wildcards
  #
  #          eg, Texas Hold-Em with 6 players and wildcards = "Deuces" has its
  #              ggstat information stored as follows:
  #              pokerdb$Texas.Hold.Em.6.4[[1]]
  #
  # 
  gamelist <- gamelist[gamelist %in% rownames(iwpgames)]
  gameref <- gsub(" ", ".", gamelist, fixed = TRUE)
  gameref <- gsub("-", ".", gameref, fixed = TRUE)
  gameref <- gsub("(", ".", gameref, fixed = TRUE)
  gameref <- gsub(")", ".", gameref, fixed = TRUE)

  # right now, only supporting 0-7 wildcards
  wclist <- list("none" = NULL, "one" = "Suicide King", 
                  "two" = "One Eyed Jacks", "three"= "", "four" = "Deuces", 
                  "five" = "", "six" = "" , "seven" = "")
  wclist$three <- c(wclist$one, wclist$two)
  wclist$five  <- c(wclist$one, wclist$four)
  wclist$six   <- c(wclist$two, wclist$four)
  wclist$seven <- c(wclist$one, wclist$two, wclist$four)
  wcmax <- length(wclist)
  
  # initialize pokerdb if needed
  if (is.null(pokerdb)) {
    pokerdb  <- list(stats = c(numdeal = numdeal, seed = seed))              
  } # end initialize pokerdb
  for (i in 1:length(gamelist)) {
    game <- iwpgames[gamelist[i], ]
    pdeck <- game$dsize - sum(game[c(2, 4, 6, 8, 10, 11)])
    pcards <- sum(game[c(1, 3, 5, 7, 9)])
    plmax  <- ifelse(floor(pdeck/pcards) >= 8, 8, floor(pdeck/pcards))
    for (j in 2:plmax) {
      for (k in 1:wcmax) {
         refgame <- paste(gameref[i], ".", j, ".", (k-1), sep="")
         if ((length(grep(refgame, names(pokerdb))) == 0)) {
           gstat <- wpstats(ngame = gamelist[i], players = j, 
                            wcards = unlist(wclist[k]),
                            numdeal = numdeal, seed = seed)
           refcmd <- paste("pokerdb$", refgame, " <- gstat", sep="")
           eval(parse(text=refcmd))
         } # end skip if it has already been calculated
## Debug code ###
# assign("debugdb", pokerdb, envir = .GlobalEnv)
# save(list = c("debugdb", "debugraw", "debughand"), file = "data/debugdb.rda")
#################
      } # end loop through wildcards
    } # end loop through players
  } # end loop through games
  pokerdb
} # end of external function wpgcreatedb

################################################################################## 
wpgraphs <- function(gstat = NULL, ngame = NULL, players = NULL, wcnum = NULL,
                     stats = FALSE, gtype = "Default", split = "Vertical") {
  #  This function generates graphs from a gstat object, or tries to match
  #  ngame/player/wcnum input to a row in the wpgamestats list
  #
  #  Args:
  #     gstat:  a game statistics list generated either from wpstats or 
  #              taken from a row of iwpgamestats.  
  #     if gstat is null, all three of the following must be populated and
  #     in combination must match a row entry of iwpgamestats when parsed
  #     ngame:   used only if gstat is null, name of poker game
  #     players: used only if gstat is null, number of players
  #     wcnum:   used only if gstat is null, number of wildcards added
  #     stats:   if TRUE the function returns the gstat parameter or 
  #              the gstat list from iwpgamestats instead of graphs 
  #     gtype:   "Default" returns a 4x4 grid with rich information
  #              "Confidence" returns either 1x1 (no split) or 1x2 (split)
  #                           with just the confidence graphs
  #              "Hands" returns either 1x1 (no split) or 1x2 (split)
  #                           with only the count of hands graphs
  #     split:   Split graphs stack "Vertical" (default) or "Horizontal"  
  #
  #  Returns:  
  #     if gstat is NULL and stats is TRUE, returns the precalculated 
  #     iwpgamestats list matching the ngame, players, and wcnum parameters
  #
  #     otherwise returns a 4x4 set of graphs including
  #       probability graph for Main Hand 
  #       win count graph for Main Hand  
  #     
  #       If split hands exist, Probability & Win Count for split hand
  #          and a one-line % pot distribution metric for winning hands
  #
  #       if no split hands exist, the a graph showing distribution of
  #          % pot for winning hands and win count/wildcard graph.
  #     

  # Validate parameters
  if (is.null(gstat)) {
    sgame <- iwpsupportedgames[ngame,]$Stats.Game
    sgame <- ifelse(is.na(sgame),ngame, sgame)
    game <- gsub(" ", ".", sgame, fixed = TRUE)
    game <- gsub("-", ".", game, fixed = TRUE)
    game <- gsub("(", ".", game, fixed = TRUE)
    game <- gsub(")", ".", game, fixed = TRUE)
    gameref <- paste(game, ".", players, ".",wcnum, sep="")
    if (length(grep(gameref, names(iwpgamestats))) > 0) {
      gstat <-  iwpgamestats[[gameref]]
    } else { #exit with error
      # ngame not supported
      if (sum(rownames(iwpgames) == sgame) != 1) {
        print(c("Supported Games:"))
        print(sort(rownames(iwpsupportedgames)))
        stop(c("ngame parameter: ", ngame, " is not supported"))
      } # end game parameter validation check
      # outside error trap, should never be called
      print("Supported game chosen but parameters not in pre-calculated list")
      print(paste("    ngame   =", ngame))
      print(paste("    players =", players))
      print(paste("    wcnum   =", wcnum))
      print("Use the following command to attempt manual generation:")
      if (is.null(wcnum)) {
      print("    wpgraph(wpstats(ngame, players))")
      } else {
      print("    wpgraph(wpstats(ngame, players, wcards = c(your wildcard choices)))")
      } # end wildcard check
      stop(paste("game reference", gameref, "not found in iwpgamestats"))     
    } # end name/player/wc check
  } # end gstat null check
  if (!is.list(gstat)) {
    if (is.null(gstat)) {
      print("If gstat parameter is null, the program requires")
      print("ngame, players and wcnum parameters to be not null")
      stop("Usage: wpgraphs(gstat, ngame, players, wcnum, gtype)")
    } else {
      print(gstat)
      stop("incorrect gstat format - must be a list")
    }
  } # not list test
  if (length(gstat$pmain) == 0 | length(gstat$cmain) == 0 |
      length(gstat$wmain) == 0 | length(gstat$psplit) == 0 |
      length(gstat$csplit) == 0 | length(gstat$wsplit) == 0 |
      length(gstat$potpct) == 0 ) {
    print(gstat)
    stop("incorrect gstat format")
  } # end gstat format check
  if (stats) {
    return(gstat)
  } # return the gstat object, not a set of graphs
  if (is.null(ngame)) {
    ngame <- gstat$game["ngame"]
  } # end ngame check
  # the max # of cards for a given hand is used to scale count graphs
  # the wildcard matrix will have 1-2 rows, grab the largest from each
  maxcard <- 0
  for (i in 1:dim(gstat$cmain[[1]])[1]) {
    maxcard <- maxcard + max(gstat$cmain[[1]][i, ])
  } # end of maxcard
  # set the base color for the graph, the default pot is 100% won
  basecol <- "dimgrey"
  # use basecol for first color, then a heatmap showing that there are
  # rival hands likely with 2 or more competitors, red = tied hand
  gcolor <- rev(heat.colors(dim(gstat$pmain[[1]])[1]-1, alpha = 1))
  gcolor <- c(basecol, gcolor)
  if (gtype == "Default") {
    par(las=2, mfcol = c(2,2), mar = c(5,4,2,2), bg = "snow2")
  } else {
    if (as.matrix(gstat$psplit[[1]])[1,1] == "") {
        par(las=2, mfcol = c(1,1), mar = c(6,4,2,2), bg = "snow2")
    } else {
       if (split == "Horizontal" ) {
        par(las=2, mfcol = c(1, 2), mar = c(6,4,2,2), bg = "snow2")
       } else {
        par(las=2, mfcol = c(2, 1), mar = c(6,4,2,2), bg = "snow2")
       } # end check for horizontal/vertical split
    } # end check for split
  } # end par settings by graph type

  stat1  <- paste(gstat[["game"]]["players"], "players &",
            gstat[["stats"]]["wcdeck"], "wildcards")

  stat2  <- paste(gstat[["stats"]]["numdeal"], 
            "hands - seed =", gstat[["stats"]]["seed"])

  stat3 <- c(paste(stat1,"-", stat2))

  if (gtype == "Default" | gtype == "Confidence") {
    if (gtype == "Default") {
      gmain <- ngame
      gsub <- "Confidence: Main Hand"
    } else {
      gmain <- paste(ngame," - Main Hand Confidence")
      gsub <- stat3
    } # end labels check
    barplot(as.matrix(gstat$pmain[[1]]), ylab = "Confidence", col=gcolor,
            main = gmain, ylim = c(0,1.09), border = NA)
    abline(h=.2, col = basecol)
    abline(h=.4, col = basecol)
    abline(h=.6, col = basecol)
    abline(h=.8, col = basecol)
    abline(h=1, col = basecol)
    legend(x="right", legend = rev(rownames(as.matrix(gstat$pmain[[1]]))), 
           fill=rev(gcolor), cex=.7, bg="white")
    text(x=dim(as.matrix(gstat$pmain[[1]]))[2] * .7,
         y=1.055,label=gsub, col = "Black")
  } # end draw main confidence graph

  if (gtype == "Default" | gtype == "Hands") {
    if (gtype == "Default") {
      gmain <- stat1
      gsub <- "Hands Won: Main Hand"
    } else {
      gmain <- paste(ngame," - Hands Won Main")
      gsub <- stat3
    } # end labels check
    barplot(as.matrix(gstat$cmain[[1]]), ylab = "Hands Won", col=gcolor,
            main = gmain, ylim = c(0, maxcard*1.2),
            border = NA)
    legend(x="right", legend = rev(rownames(as.matrix(gstat$cmain[[1]]))), 
           fill=rev(gcolor), cex=.7, bg="white")
    text(x=dim(as.matrix(gstat$cmain[[1]]))[2] * .8,
         y=maxcard*1.1,label=gsub, col = "Black")

  } # end draw main Hands graph 

  if (as.matrix(gstat$psplit[[1]])[1,1] == "") {
    if (gtype == "Default") {
      pcolor <- rev(heat.colors(length(gstat$potpct)-1, alpha = .5))
      pcolor <- rev(c(basecol, pcolor))
      barplot(gstat$potpct, ylab = "Confidence", col = pcolor,
              main = "No Split Hand", ylim = c(0,1.09), border = NA)
      text(x=length(gstat$potpct) * .7,
           y=1.055,label="Percentage of Pot Won", col = "Black")
      wcolor <- rainbow(dim(as.matrix(gstat$wmain[[1]]))[1], start = .6)
      barplot(as.matrix(gstat$wmain[[1]]), ylab = "Hands Won", col=wcolor,
           ylim = c(0, maxcard*1.2), main = stat2, border = NA)
      legend(x="right", legend = rev(rownames(as.matrix(gstat$wmain[[1]]))), 
             fill=rev(wcolor), cex=.8, bg="white")
      text(x=dim(as.matrix(gstat$wmain[[1]]))[2] * .7,
           y=maxcard*1.1,label="Hands Won: # Wildcards", col = "Black")
    } # if not Default type, these graphs aren't drawn
  } else {
    # split hands generate interest in normal % of pot won
    pctpot <- rev(sort(round(gstat$potpct,3)))
    pctlab <- names(pctpot)
    potlab <- paste(pctlab, "(", pctpot, ")", sep="")
    if (length(potlab) > 3) {
      potlab <- potlab[1:3]
    } # end large number of pots check
    potlab <- paste(potlab, collapse=", ")
    # split hand may have entries with no actual values
    scidx <- NULL
    for (i in 1:dim(gstat$psplit[[1]])[1]) {
      if (sum(gstat$psplit[[1]][i,]) > 0) {
       scidx <- c(scidx, rownames(gstat$psplit[[1]][i,]))
      } # end of test rows for values
    } # end of loop through psplit rows
    scolor <- rev(heat.colors(length(scidx)-1, alpha = 1))
    scolor <- c(basecol, scolor)
    if (gtype == "Default" | gtype == "Confidence") {
      maxcard <- 0
      for (i in 1:dim(gstat$csplit[[1]])[1]) {
        maxcard <- maxcard + max(gstat$csplit[[1]][i, ])
      } # end of maxcard
      if (gtype == "Default") {
        gmain <- paste("%Pot", potlab)
        gsub <- "Confidence: Split Hand"
      } else {
        gmain <- paste(ngame," - Split Hand Confidence")
        gsub <- paste("%Pot", potlab)
      } # end labels check
      barplot(as.matrix(gstat$psplit[[1]]), ylab = "Confidence", col=scolor,
              main = gmain, ylim = c(0,1.09), border = NA)
      abline(h=.2, col = basecol)
      abline(h=.4, col = basecol)
      abline(h=.6, col = basecol)
      abline(h=.8, col = basecol)
      abline(h=1, col = basecol)
      legend(x="right", legend = rev(scidx), 
             fill=rev(scolor), cex=.8, bg="white")
      text(x=dim(as.matrix(gstat$psplit[[1]]))[2] * .7,
           y=1.055,label=gsub, col = "Black")
    } # end draw Split Confidence graph   
    if (gtype == "Default" | gtype == "Hands") {
      if (gtype == "Default") {
        gmain <- stat2
        gsub <- "Hands Won: Split Hand"
      } else {
        gmain <- paste(ngame," - Hands Won Split")
        gsub <- paste("%Pot", potlab)
      } # end labels check
      barplot(as.matrix(gstat$csplit[[1]]), ylab = "Hands Won", col=scolor,
           ylim = c(0, maxcard*1.2), main = gmain, border = NA)
      legend(x="right", legend = rev(scidx), 
             fill=rev(scolor), cex=.8, bg="white")
      text(dim(as.matrix(gstat$csplit[[1]]))[2] * .7,
         y=maxcard*1.1,label=gsub, col = "Black")
    } # end draw Split Hands graph   
  } # end of split hand check

} # end of public function wpgraphs





