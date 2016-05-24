################################################################################
#
#  This file contains functions that operate on deck objects, and calculate
#  values for each distinct player (as defined by a "hand" column entry which
#  begins with "P" and is not "PF" or "PS").  These functions should not be
#  called until wpdeal has operated on a deck.
#
#  There is one external function in this file "wpmakehand" which uses all the 
#  other utilities to pick the best hand based on game and split
#  (from game$hand+game$hmod+game$hmin or game$split+game$smod+game$smin) 
#
#  All of these functions return a three character vector, 
#    index1 = value of hand (relative to other hands of the same type)
#    index2 = string showing the cards and wildcards used to make the hand
#    index3 = hand shortname  (eg, "Straight", "3-Kind")
#
#  All high functions will calculate best wildcard hand of the type even if the
#  same number of wildcards could generate a better hand.  Low functions will
#  use wildcards to break up any hand higher than 1 of a kind.  Wildcards do
#  not apply for hand modifier Hole, but do for Sum (14, 10 if hand min= Spots) 
#
#  Functions prefixed with "h" chose the highest value hand if multiple exist
#  iwpmakehands has all of the "split" logic for sum, hole card or low hand
#  and is the only function called by anything else - directly by iwpshowdown
#  and indirectly by the stats functions
#
#

################################################################################
hflush <- function(hand) {
  #  This is a utility function used to calculate if a flush
  #  is possible, and calculate the highest value flush. 
  #
  #  args - a "hand" or subset("deck", hand == player) data frame
  #         low=TRUE will find the LOWEST VALUE flush 
  # 
  #  returns - vector with 3 values:
  #       1   0 if a flush is not possible, or 5000+ sum of all card values
  #       2   list of rownames of the cards used to make up the flush or ""
  #           wildcards are flagged with a "w". Wildcards are all Ace=14 value
  #       3   "Flush"
  #       eg: [1] "0"        ""         "Flush"
  #           [1] "5006"           "AH 2D 3D 4D 6D" "Flush" 
  #           [1] "5014"            "9Dw JH QH KH AH" "Flush"  
  #

  ###########################################################
  # Wildcard preprocessor - all hands will use WC if avialable
  wsum <- sum(hand$value == "W")
  if (wsum > 0) {
    widx <- grep("W", hand$value)
    for (i in 1:ifelse(wsum > 5,5,wsum)) {
      wcard <- paste(rownames(hand[widx[i], ]),"w", sep="")
      if (i == 1) {
        wFlush <- wcard
      } else {
        wFlush <- paste(wcard, wFlush)
      } # end first card check
    } # end build hand of wildcards 
  } else {
    wFlush <- NULL
  } # end no wild card check

  #########
  # verify that there are at least 5 cards in the hand
  fnum <- sort(as.numeric(as.character(hand[hand$value != "W",]$value)))
  vFlush <- 0
  hFlush <- 0
  if ((length(fnum) + wsum) >= 5) {
    #################
    # group the suits
    if (wsum < 4) {
      suits <- table(hand$suit)
      if (!is.na(suits["W"])) {
       suits <- suits[1:(length(suits)-1)] + suits["W"]
      } # end suit creation
      flsuits <- suits[suits >= 5]
      ##################
      # if flushes exist, loop through possibilities and calculate strongest
      if (length(flsuits) > 0) {
        for (i in names(flsuits)) {
          fhand <- hand[hand$suit == i | hand$suit == "W", ]
          xhand <- fhand[fhand$suit == i, ]
          fnum <- sort(as.numeric(as.character(
                  fhand[fhand$value != "W",]$value)), decreasing = TRUE)
          fval <- 5000 + 14*wsum + sum(fnum[1:(5-wsum)])
          if (fval > vFlush) {
            vFlush <- fval
            cidx <- NULL
            for (j in 1:(5-wsum)) {
              cidx <- c(cidx, grep(valface(fnum[j]),rownames(xhand)))
            } # end build cidx
            hFlush <- paste(c(wFlush, rownames(xhand[cidx, ])), 
                                      collapse = " ")
          } # end of check to see if suit has higher flush value
        } # end of loop through suits
      } # end of check for potential flush suits   
    } else {
  # Note, it is unlikely that this function will be called if there are 4 or
  # more wildcards, because even a straight flush is defeated by a 5 of a kind, 
  # but just in case the function is called, it will return the correct result
      if (wsum == 4) {
        vFlush <- 5000 + max(fnum) + 4*14
        cidx <- grep(valface(max(fnum)), rownames(hand))
        hFlush <- paste(rownames(hand[cidx, ]), wFlush)
      } else { # wsum 5 or greater (hFlush is already built)
        vFlush <- 5000 + 14*5
        hFlush <- wFlush
      } # end 4 wildcard check
    } # end wildcard less than 4 check
  } # end of check for at least 5 cards
  c(vFlush, hFlush, "Flush")
} # end function hflush 

################################################################################
hstraight <- function(hand) {
  #  This is a utility function used to calculate if a straight
  #  is possible, and the highest value straight if possible
  #
  #  args - a "hand" or subset("deck", hand == player) data frame 
  # 
  #  returns - vector with 3 values:
  #       1   0 if a straight is not possible, otherwise 4000+
  #            a number from 5 to 14, indicating the top card (as character)
  #       2   list of rownames of the cards used to make up the straight or ""
  #           wildcards are flagged with a "w".  Aces can be a "1" or a "14"
  #       3   "Flush"
  #       eg: [1] "0"     ""      "Flush"
  #           [1] "5023"           "2D 3D 4D 5D 9D" "Flush"   
  #           [1] "5064"            "JH QH KH AH 9Dw" "Flush" 
  #  
  if (dim(hand)[1] < 5) {
    vStraight <- 0
  } else {
    sWild <- sum(hand$value == "W")
    widx <- grep("W", hand$value)
    ##################################################
    # calculate the best straight value
    if (sWild < 5) {
      sNum <- sort(as.numeric(as.character(
                   unique(hand[hand$value != "W",]$value))),
                   decreasing = TRUE)
      vStraight <- 0
      # Ace can be either 1 or 14
      if (sNum[1] == 14) {
        sNum <- c(sNum, 1)
        vace <- -1
      } else {
        vace <- 0
      } # end Ace check 
      if (length(sNum) + sWild + vace > 4) {  
        for (i in 1:(length(sNum) - 4 + sWild)) {
          if ((sNum[i] - sNum[i + 4 - sWild]) <= 4) {
            vStraight <- sNum[i + 4 - sWild] + 4
            break
          } # end test to allow a straight 
        } # end loop through top cards in hand
      } # not enough unique cards for straight
    # Note, it is unlikely that this function will be called if there are 4 or
    # more wildcards, because even a straight flush is defeated by a 5 of a kind, 
    # but just in case the function is called, it will return the correct result
    } else { 
      # Five wildcards = 14 high straight
      vStraight <- 14
    } # end wildcard check
  } # end of cards < 5 check
  ##################################################
  # Build the hand
  #   if higher than 14, return 14
  vStraight <- ifelse(vStraight > 14, 14, vStraight)
  if (vStraight == 0) {
    rval <- c("0","")
  } else {
    hStraight <- NULL
    for (i in (vStraight - 4):vStraight) {
      cmatch <- ifelse(i == 1, valface(14), valface(i))
      cidx <- grep(cmatch, rownames(hand)) 
      if (length(cidx) == 0) {
        cStraight <- paste(rownames(hand[widx[sWild], ]), "w", sep="")
        sWild <- sWild - 1
      } else { # add the first matched card
        cStraight <- rownames(hand[cidx[1], ])
      } # end of match with wildcard
      if (is.null(hStraight)) {
        hStraight <- cStraight
      } else {
        hStraight <- paste(hStraight, cStraight)
      } # end of first card in hand check
    } # end of building straight hand
    rval <- c(vStraight + 4000, hStraight)
  } # end of no straight check
  c(rval, "Straight")
} # End function hstraight 

################################################################################
hstrflush <- function(hand) {
  #  This is a utility function used to calculate if a straight flush
  #  is possible, and calculate the highest value straight flush if possible
  #
  #  args - a "hand" or subset("deck", hand == player) data frame 
  # 
  #  returns - vector with 3 values:
  #       1   0 if a straight flush is not possible, otherwise
  #           a number from 5 to 14, indicating the top card (as character)
  #       2   list of rownames of the cards used to make up the strflush or ""
  #           wildcards are flagged with a "w". Ace can equal 1 or 14
  #       3   "StrFlush"
  #       eg: [1] "0"        ""         "StrFlush"
  #           [1] "8005"           "AD 2D 3D 4D 5D" "StrFlush"
  #           [1] "8014"            "9Dw JH QH KH AH" "StrFlush" 
  #

  rval <- c("0", "")
  suits <- table(hand$suit)
  if (!is.na(suits["W"])) {
    suits <- suits[1:(length(suits)-1)] + suits["W"]
  } # end suit creation
  flsuits <- suits[suits >= 5]
  if (length(flsuits) > 0) {
    for (i in names(flsuits)) {
      fhand <- hand[hand$suit == i | hand$suit == "W", ]
      sfval <- hstraight(fhand)[1:2]
      sfval[1] <- ifelse(sfval[1] == "0", "0", as.numeric(sfval[1])+4000)
      if (as.numeric(sfval[1]) > as.numeric(rval[1])) {
        rval <- sfval
      } # end check for higher hand value
    } # end straight check of flush suits
  } # end no flush suits check
  c(rval, "StrFlush")
} # End function hstrflush

################################################################################
hofkind <- function(hand) {
  #  This is a utility function used to calculate how many cards of a kind
  #  exist in a hand, and choose the best hand including kickers.  This 
  #  function ignores straight, flush and straight/flush combinations.
  #
  #  args - a "hand" or subset("deck", hand == player) data frame 
  # 
  #  returns - vector with 3 values:
  #       1   wildcards are always part of the "ofKind" card group
  #           if any wildcards are present, 2-Pair and HighCard can't occur
  #           value is based on the third vector, as follows:
  #           5-Kind    = 9000 + 50 * ofKind card value (five wc = 5 aces)
  #           4-Kind    = 7000 + 50 * ofKind card value + Kicker value
  #           FullHouse = 6000 + 50 * 3ofKind card value + Pair value
  #           3-Kind    = 3000 + 50 * ofKind card value + sum(Kicker values)
  #           2-Pair    = 2000 + 50 * 1stPair + 15 * 2nd Pair + Kicker value 
  #           Pair      = 1000 + 50 * ofKind card value + sum(Kicker values)
  #           HighCard  =    0 + 50 * ofKind card value + sum(Kicker Values)
  #       2   list of rownames of the cards used to make up the hand
  #           wildcards are flagged with a "w". Ace always equals 14
  #       3   5-Kind, 4-Kind, FullHouse, 3-Kind, 2-Pair, 2-Kind, HighCard
  #       eg: "745" "AH KH QH JH 9D"       "A-High" 
  #           "1736" "AH 9Dw KH QH JH"            "Pair" 
  #           "2739" "KH KS 5D 5S AH"         "2-Pair" 
  #           "3277" "5D 5S KSw AH KH"          "3-Kind" 
  #           "6655" "KH KS 9Dw 5D 5S"       "FullHouse"
  #           "7664" "KH KS 2Dw 9Dw AH"           "4-Kind" 
  #           "3650" "7C 7D 7H 3Sw 3Hw"         "5-Kind" 
  #
  
  if (dim(hand)[1] < 5) {
    hm <- dim(hand)[1]
  } else {
    hm <- 5
  } # end hand length check

  ###############################
  # Build wildcard tracking values
  wnum <- sum(hand$value == "W")
  if (wnum > 0) {
    widx <- grep("W", hand$value)
    wvect <- paste(rownames(hand[widx, ]), "w", sep="")
    for (i in 1:ifelse(wnum > hm, hm, wnum)) {
      if (i == 1) {
        whand <- wvect[i]
      } else {
        whand <- paste(whand, wvect[i])
      } # end first card check
    } # end build wildcard hand
  } else { # no wildcards
    widx <- NULL
    whand <- NULL
  } # end zero wildcard check

  ############################
  # begin rval, rhand, rname calculations
  if (wnum >= hm) {
    # handle "all cards are wild" case
    rhand <- whand
    if (hm == 5) {
      rval <- 9000 + 50 * 14
      rname <- "5-Kind"
    } else {
      if (hm == 4) {
        rval <- 7000 + 50 * 14
        rname <- "4-Kind"
      } else {
        if (hm == 3) {
          rval <- 3000 + 50 * 14
          rname <- "3-Kind"
        } else {
          if (hm == 2) {
            rval <- 1000 + 50 * 14
            rname <- "Pair"           
          } else {
            rval <- 50 * 14
            rname <- "A-High"   
          } # end 2 wildcard hand
        } # end 3 wildcard hand
      } # end 4 wildcard hand
    } # end 5 wildcard hand"
  } else {
    ####### 
    # Build the best "of kind" hand
    knam <- names(table(hand$value))
    kind <- table(hand$value)[knam[knam != "W"]] + wnum
    kind <- kind[kind == max(kind)]
    mxkind <- c(kind[1], max(as.numeric(names(kind))))
    mxvect <- rownames(hand[hand$value == mxkind[2], ])
    mxkind[1] <- ifelse(mxkind[1] > hm, hm, mxkind[1])
    for (i in 1:(mxkind[1]-wnum)) {
      if (i == 1) {
        mxhand <- mxvect[i]
      } else {
        mxhand <- paste(mxhand, mxvect[i])
      } # end first card check
    } # end build mxhand
    if (!is.null(whand)) {
        mxhand <- paste(mxhand, whand)
    } # add wildcards to mxhand
    ################ 
    # Test 5 of Kind
    if (mxkind[1] >= 5){
      rval <- 9000 + 50 * mxkind[2]
      rhand <- mxhand
      rname <- "5-Kind"
    } else {
      ################ 
      # build next portion of hand
      nkind <- table(hand$value)[knam[knam != "W" & knam != mxkind[2]]]
      nkind <- nkind[nkind > 0]
      ################
      # Test 4 of kind
      if (mxkind[1] == 4){
        if (hm > 4 ) {
          kval <- c(max(as.numeric(names(nkind))))
          khand <- rownames(hand[hand$value == kval, ])[1]
          rhand <- paste(mxhand, khand)
        } else {
          kval <- 0
          rhand <- mxhand
        } # end at least 5 cards check
        rval <- 7000 + 50 * mxkind[2] + kval
        rname <- "4-Kind"
      } else {
        ################ 
        # check for pairs in next portion of hand
        if (length(nkind) > 0) {
          nxkind <- nkind[nkind == max(nkind)]
          nxkind <- max(as.numeric(names(nxkind))) 
          if (max(nkind) > 1) { # choose the best pair & kicker
            nxvect <- rownames(hand[hand$value == nxkind, ])
            nxhand <- nxvect[1]
            nxhand <- paste(nxhand, nxvect[2])
          ################ 
          # generate kicker values and hands
            knam <- names(nkind)
            kval <- knam[knam != nxkind]
            if (length(kval) > 0 ) {
              kval   <- sort(as.numeric(kval), decreasing = TRUE)[1]
              khand  <- rownames(hand[hand$value == kval, ])[1]
            } else {
              kval <- 0   
            } # end of no kickers check
          } else { # no pairs in rest of hand, and at least 2 kickers
            nxhand <- NULL
            kvect <- sort(as.numeric(names(nkind)), decreasing = TRUE)
            kval <- sum(kvect[1:(hm-mxkind[1])])
            khand <- rownames(hand[hand$value == kvect[1], ])
            if ((hm-mxkind[1]) > 1) {
              for (i in 2:(hm-mxkind[1])) {
                khand <- paste(khand, 
                               rownames(hand[hand$value == kvect[i], ]))
              } # end of build kicker hand
            } # end kickers > 1 test
          } # end of pairs test in kickers
        } else { # no cards left in hand
            nxhand <- NULL
            kval <- 0       
        } # end of no cards in hand check
        ################
        # Test Full House
        if (mxkind[1] == 3 & !is.null(nxhand)){
          rval <- 6000 + 50 * mxkind[2] + nxkind
          rhand <- paste(mxhand, nxhand) 
          rname <- "FullHouse" 
        } else {
          ################
          # Test 3-Kind
          if (mxkind[1] == 3) {            
            rval <- 3000 + 50 * mxkind[2] + kval
            if (kval == 0) {
              rhand <- mxhand
            } else {
              rhand <- paste(mxhand, khand)
            } # end of null kicker check 
            rname <- "3-Kind" 
          } else {
            ################
            # Test 2-Pair
            if (mxkind[1] == 2 & !is.null(nxhand)) {                
              rval <- 2000 + 50 * mxkind[2] + 15 * nxkind + kval
              if (kval == 0) {
                rhand <- paste(mxhand, nxhand)
              } else { 
                rhand <- paste(mxhand, nxhand, khand)
              } # end of null kicker check  
              rname <- "2-Pair" 
            } else {
              ################
              # Test Pair
              if (mxkind[1] == 2) {
                rval <- 1000 + 50 * mxkind[2] + kval
                if (kval == 0) {
                  rhand <- mxhand
                } else {
                  rhand <- paste(mxhand, khand)
                } # end of null kicker check 
                rname <- "Pair" 
              } else { 
              ################
              # High Card  
                rval <- 50 * mxkind[2] + kval
                if (kval == 0) {
                  rhand <- mxhand
                } else {
                  rhand <- paste(mxhand, khand)
                } # end of null kicker check 
                rname <- paste(valface(mxkind[2]), "-High", sep="") 
              } # end Pair check
            } # end 2-Pair check
          } # end 3-Kind check  
        } # end FullHouse check
      } # end 4-Kind check
    } # end 5-Kind check
  } # end all wildcard check
  c(rval, rhand, rname)
} # end function hofkind

################################################################################
iwpmakehand <- function(hand, ngame, split=FALSE) {
  #  This is a utility function used to calculate how many cards of a kind
  #  exist in a hand, and choose the best hand including kickers.  This 
  #  function ignores straight, flush and straight/flush combinations.
  #
  #  args - a "hand" or subset("deck", hand == player) data frame
  #         ngame - a legal game from "rownames(wpsupportedgames)" 
  #         split - if FALSE, use game$hand variables, 
  #                 if TRUE use game$split variables 
  # 
  #  returns - vector with 3 values:
  #       1   Value of hand
  #       2   list of rownames of the cards used to make up the hand
  #           wildcards are flagged with a "w". Ace varies by game
  #       3   Name of the hand returned    
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
  if (!((sum(names(hand) == names(iwpstddeck)) 
         == dim(iwpstddeck)[2]))) {
    print("Structure of Hand parameter")
    print(str(hand))
    print("Structure of the standard deck")
    print(str(iwpstddeck))
    stop("Hand less than 5 cards or does not match Deck structure")
  } # end hand structure check 

  ######################
  # load the game type variables, based on game and split parameters
  if (split) { # make the split hand
    htyp <- game$split
    hmod <- game$smod
    hmin <- game$smin 
  } else { # make the primary hand
    htyp <- game$hand
    hmod <- game$hmod
    hmin <- game$hmin 
  } # end of game variable assignment
  # Community game logic is not part of this function
  hmod <- ifelse(hmod == "PC" | hmod == "PF", "None", hmod)
 

  rval <- 0
  rhand <- ""
  rname <- "NoHand"
  ######################
  # split by hand type
  if (htyp == "HI" | hmod == "Hole") {
    switch(hmod,
      ###########
      "Hole" = {
        hhand <- hand[hand$holec, ]
        if (length(grep(hmin, c("D","S","C","H")) > 0)) {
          hhand <- subset(hhand, 
                   namevals(rownames(hhand), "Suits") == hmin)
        } else {
          hmin <- ""
        } # end filter by suit check
        if (dim(hhand)[1] > 0) {
          if (htyp == "LO") {
            if (max(as.numeric(as.character(
                namevals(rownames(hhand), "Values")))) == 14) {
              rval <- 14
            } else{
              rval <- min(as.numeric(as.character(
                          namevals(rownames(hhand), "Values"))))
            }
          } else {
            rval <- max(as.numeric(as.character(
                        namevals(rownames(hhand), "Values"))))
          } # end hi-lo check
          rhand <- ifelse(hhand$suit != "W", rownames(hhand),
                          paste(rownames(hhand), "w", sep=""))
          rname <- paste(valface(rval), hmin, "-inHole", sep="")
          if (htyp == "LO" & rval == 14) {
            rval <- 1
          } # end Ace check
        } # end check for any matching hole cards
      }, # end hmod = Hole logic
      ###########
      "Sum" = {
        if (hmin == "Spots") {
          # throw out face cards, set wildcard=10, ace=1
          shand <- hand[!(hand$value %in% 11:13), ]
          if (sum(shand$suit == "W") > 0) {
            shand[shand$suit == "W", "value"] <- 10
          } # end wildcard check
          shand$numval <- as.numeric(as.character(shand$value))
          if (length(shand[shand$numval == 14, ]$numval) > 0) {
            shand[shand$numval == 14, ]$numval <- 1
          } # end ace = 1
          hmin <- "+Spots"
        } else {       
          # sum by suit
          if (length(grep(hmin, c("D","S","C","H")) > 0)) {
            shand <- hand[hand$suit == hmin | hand$suit == "W", ]
            hmin <- paste("+", hmin,  sep="")
          # sum everything
          } else {
            shand <- hand
            hmin <- "+Sum"
          } # end filter by suit check
          # wildcards count as aces
          if (sum(shand$suit == "W") > 0) {
            shand[shand$suit == "W", "value"] <- 14
          } # end wildcard check
          shand$numval <- as.numeric(as.character(shand$value))
        } # end filter by Spot check
        if (dim(shand)[1] > 0) {
          rval <- sum(shand$numval)
          if (rval < 10) {
            rtyp <- as.character(rval)
          } else {
            if (rval == 70 & hmin == "+Spots") {
              rtyp <- "60"
            } else {
              rtyp <- paste(substr(as.character(rval),1,1),"0", sep="")
            } # end maximum value for spots check
          } # end set rtyp
          rhand <- ifelse(shand$suit != "W", rownames(shand),
                          paste(rownames(shand), "w", sep=""))
          rname <- paste(as.character(rtyp), hmin, sep="")
        } # end check hand has at least one card
      }, # end hmod = Sum logic
      ###########
      # Basic normal high hand logic
      "None" = { 
        # Work out the best "of Kind/Full House/2 Pair" type hand
        ofkind <- hofkind(hand)
        rval <- ofkind[1]
        rhand <- ofkind[2] 
        rname <- ofkind[3]
        # check to see if it is a strflush
        if (as.numeric(ofkind[1]) < 9000) {
          strflush <- hstrflush(hand)
          if (as.numeric(strflush[1]) > 0) {
            rval <- strflush[1]
            rhand <- strflush[2] 
            rname <- strflush[3]  
          } else {
            # check to see if it is a flush
            if (as.numeric(ofkind[1]) < 6000) {
              strflush <- hflush(hand)
              if (strflush[1] > 0) {
                rval <- strflush[1]
                rhand <- strflush[2] 
                rname <- strflush[3] 
              } else {
              # check to see if it is a straight
                strflush <- hstraight(hand)
                if (as.numeric(strflush[1]) > 0) {
                  rval <- strflush[1]
                  rhand <- strflush[2] 
                  rname <- strflush[3] 
                } # end Straight check
              } # end Flush check
            } # end 4-Kind or Full house check       
          } # end straight flush check
        } # end 5-Kind check
        ######
        # Ensure hand is good enough to win
        if (hmin == "Pair" & as.numeric(rval) < 1000) {
          rname <- "NoHand"
        } # end minimum = Pair logic
      }, # end hmod = "None" logic
    ) # end switch hmod
  } else {
  #############
  # Begin low hand logic
    # Hands smaller than 5 cards are "No Hand" for low logic
    if (dim(hand)[1] >= 5) {
      ###############################
      # Build wildcard tracking values
      wnum <- sum(hand$value == "W")
      if (wnum > 0) {
        widx <- grep("W", hand$value)
        wvect <- paste(rownames(hand[widx, ]), "w", sep="")
      } else { # no wildcards
        widx <- NULL
        wvect <- NULL
      } # end zero wildcard check
      if (wnum >= 5) {
        rval <- ifelse(hmod == ".A-5", sum(1:4) + 5 * 50,
                  ifelse(hmod == ".A-6", sum(1:4) +6 * 50, 
                                       sum(2:5) + 7* 50)) 
        rhand  <- wvect[1:5]
        rname <- ifelse(hmod == ".A-5", "5-High",
                  ifelse(hmod == ".A-6", "6-High", "7-High"))
      } else {
        ##########
        # Assemble non-wildcard portion of the hand
        vhand <- hand[hand$value != "W", ]
        # handle Aces in the hand that aren't wild
        vmin <- ifelse(hmod == ".A-5" | hmod == ".A-6",1 , 2)
        vcard <- table(vhand$value)
        vcard <- vcard[vcard > 0]
        vval <- sort(as.numeric(names(vcard)))
        lidx <- NULL
        lrtn <- NULL
        lval <- NULL
        vvc <- 1
        vwc <- 1
        lfv <- 15
        lfi <- 0
        ##########
        # Build the best "number-High" hand
        for (i in vmin:14) {
          # Match the lowest nonwildcard
          if ((i == 1 & vval[length(vval)] == 14) |
              (i != 1 & vval[vvc] == i)) {
            if (i == 1) {
              lidx <- c(lidx, grep(valface(14), rownames(hand))[1])
              if (length(vval) == 1) {
                vval <- 1
              } else {
                vval <- c(1, vval[1:(length(vval) - 1)])
              } # end vval check
            } else {
              lidx <- c(lidx, grep(valface(i), rownames(hand))[1])
            } # end of ace=1 or 14 logic
            vvc <- ifelse(vvc == length(vval), vvc, vvc + 1)
            lval <- c(lval, i)
          } else {
          # match a wildcard if it exists
            if (wnum >= vwc) {
              lidx <- c(lidx, widx[vwc]) 
              vwc <- vwc + 1
              lval <- c(lval, i)
            } # end check for wildcard
          } # end check for number match
          ########
          # check to see if hand is complete
          if (length(lidx) == 5) {
            # make sure no straights or flushes exist, if they count as high
            # ".A-5" ignores flushes and straights
            if (hmod != ".A-5") {
            ########## 
            # Flush check - wildcards in hand will cause check to fail
              if (max(table(hand[lidx, ]$suit)) == 5) {
                # if hand = 5, stuck with either flush or strflush
                if (dim(hand)[1] == 5) {
                  lrtn <- hstrflush(hand[lidx, ])
                  if (lrtn[1] == "0") {
                    lrtn <- hflush(hand[lidx, ])
                  } # end not straightflush check
                  break # this is the final card
                } # end 5-card hand check
                # Try to find a card with different suit in current nums
                thand <- hand[hand$value %in% hand[lidx,]$value, ] 
                if (dim(thand)[1] == 5) {
                  # no alternative cards exist with same numbers
                  # remove the top card from the stack & pick new card
                  if (lfv > lval[5]) {
                    lrtn <- hstrflush(hand[lidx, ])
                    if (lrtn[1] == "0") {
                      lfv <- lval[5]
                      lfi <- lidx[5]
                    } # end not straightflush check
                  } # handle multiple flush card case
                  if (i < max(vval)) {
                    lidx <- lidx[1:4]
                    lval <- lval[1:4]
                    lrtn <- NULL
                    next # loop to next card
                  } else { # all cards are flush, use lowest value
                    lidx <- c(lidx[1:4], lfi)
                    lval <- c(lval[1:4], lfv)
                    lrtn <- hflush(hand[lidx, ])
                    break  # this is the final card                 
                  } # end loop check
                } else {
                  # Replace with a card that breaks the flush
                  tval <- thand[thand$suit != hand[lidx, ]$suit[1], ][1, 1]
                  tidx <- grep(valface(tval), rownames(hand))
                  lidx <- replace(lidx,grep(tidx[1],lidx),tidx[2])
                  rval <- i*50 + sum(lval[1:4]) 
                  rhand <- ifelse(hand[lidx,"value"] == "W", 
                             paste(rownames(hand[lidx,]), "w", sep=""), 
                             rownames(hand[lidx, ]))
                  rname <- paste(valface(i),"-High", sep= "")
                  break  # this is the best hand
                } # end of find same card different suit      
              } # end flush check
              ########## 
              # Straight check 
              thand <- hand[lidx, ]
              thand$value <- lval  # <- real values subbed for wc values
              lrtn <- hstraight(thand)
              # the replaced values confuse hstraight hand type logic
              if (lrtn[1] != "0") {
                # special case - A-5 straight, set ace high instead of low
                # if hand is just A-5 or we've cycled through to A=14
                # A-high beats any other hand that can be constructed
                if ((hmod != ".A-5") & lrtn[1] == "4005" & wnum == 0 &
                     (max(vval) == 5 | i == 14)) {
                    rval <- 14*50 + sum(2:5)
                    if (max(vval) == 5) {
                      lidx <- c(lidx[2:5], lidx[1])
                    } # end of Ace low check 
                    rhand <- rownames(hand[lidx, ])
                    rname <- "A-High"
                    lrtn <- NULL
                    break
                } # return A-High hand
                # if no wildcards and hand = 5, stuck with straight
                if ((wnum == 0) & (dim(hand)[1] == 5)) {
                  break
                } # end stuck with straight
                ###########################
                # The best way to break a straight is trim the 5th card
                # and add a new card one higher than the 5th.  If there
                # is a wild card or higher card it will come up next
                if ((i < max(vval) & (wnum == 0)) | (wnum >= vwc)) {
                  lidx <- lidx[1:4]
                  lval <- lval[1:4]
                  lrtn <- NULL
                  next                    
                } else {
                  if ((wnum == 0)) {
                  # if we are here, the hand has nothing in it but a
                  # mix of straights and pairs (or 3/4kind).  We want 
                  # to find the lowest pair we can use to replace the 
                  # highest value nonpair card
                    lowpair <-  min(as.numeric(names(
                                table(hand$value)[table(hand$value) > 1])))
                    lpidx <- grep(valface(lowpair), rownames(hand))[2]
                    if (hand[lidx[5], "value"] != hand[lpidx,"value"]) {
                      lidx <- replace(lidx,5,lpidx)
                    } else {
                      lidx <- replace(lidx,4,lpidx)
                    } # end low pair is last card check
                    lrtn <- hofkind(hand[lidx, ])
                    break  # this is the best hand
                  } else {
                  # If we get here, the straight has a wildcard in it
                  # but there are no other wildcards in the dex and
                  # the next (i+1) value isn't matched in vhand, so
                  # the best option is to shift the highest wildcard
                  # to the i+1 position.  As wildcards are chosen in
                  # order, the one to shift is widx[wnum]
                    tidx <- grep(widx[wnum], lidx)
                    if (tidx == 1) { # swap first to last
                      lidx <- c(lidx[2:5], widx[wnum])
                      lval <- c(lval[2:5], 7)
                    } else {
                      if (tidx < 5) { # shift highest card
                        lidx <- c(lidx[1:(tidx-1)], 
                                  lidx[(tidx + 1):5], widx[wnum]) 
                      } # if tidx is the last card, don't change index
                      lval <- c(lval[1:(tidx-1)], 
                                lval[(tidx + 1):5], lval[5]+1)                  
                    } # end shift positions
                    lrtn <- NULL
                    rval <- (i+1) * 50 + sum(lval[1:4]) 
                    rhand <- ifelse(hand[lidx,"value"] == "W", 
                             paste(rownames(hand[lidx,]), "w", sep=""), 
                             rownames(hand[lidx, ]))
                    rname <- paste(valface(i+1),"-High", sep= "")
                    break
                  } # end wc in straight check
                } # end break out of straight logic
              } else { # no straight
                lrtn <- NULL
              } # end straight check
            } # end .A-5 check
            # this is the lowest possible hand
            if (is.null(lrtn)) {
              rval <- i*50 + sum(lval[1:4]) 
              rhand <- ifelse(hand[lidx, "value"] == "W", 
                         paste(rownames(hand[lidx, ]), "w", sep=""), 
                         rownames(hand[lidx, ]))
              rname <- paste(valface(i),"-High", sep= "")
              break
            } # end check to see if straight/flush/strflush matched
          } else { 
            ######################################
            ## stop the loop if out of wildcards and have run through
            ## all of the unique card values
            if ((vwc > wnum )& max(vval) <= i) {
                break
            } # end of can't make a 5 unique card low hand
          } # end of hand complete check
        } # end loop through card values
        ##############
        ## check for pairs or other "ofkind" hands
        ## all wildcards have already been assigned
        if (length(lidx) < 5) {
          vcard <- vcard[vcard > 1]
          vval <- sort(as.numeric(names(vcard)))
          vidx <- grep(vval[1], hand$value)
          # this is a case where all remaining cards match 1 card
          if (length(vcard) == 1 | length(lidx) == 4) {
            vlen <- 5 - length(lidx[lidx != vidx[1]])
            lval <- c(lval[lidx != vidx[1]], replicate(vlen, vval[1]))
            lidx <- c(lidx[lidx != vidx[1]], vidx[1:vlen])
          } else {
            nidx <- grep(vval[2], hand$value)
            if (length(vval) > 2 & lval[1] == vval[1] & lval[2] == vval[2]) {
            # This is a "3 pair".  Use largest card as kicker
            # and assemble the 2 higher cards
              lval <- c(replicate(2, vval[2]), replicate(2, vval[1]), vval[3])
              lidx <- c(nidx[1:2], vidx[1:2], grep(vval[3], hand$value)[1])
            } else {
              # pull the lowest 2 pair and add to the kicker, if any
              lval <-  c(lval[!(lidx %in% c(vidx[1], nidx[1]))],
                         replicate(2, vval[1]), replicate(2, vval[2]))
              lidx <-  c(lidx[!(lidx %in% c(vidx[1], nidx[1]))],
                          vidx[1:2], nidx[1:2])
              # if we don't have 5 cards yet, it's a natural full house
              if (length(lidx) < 5)  {
                # if the lowest card has 3 cards, use it, else use the next
                if (length(vidx) > 2) {
                  lval <- c(vval[1], lval)
                  lidx <- c(vidx[3], lidx)
                } else {
                  lval <- c(lval, vval[2])  
                  lidx <- c(lidx, nidx[3])
                } # end check lowest for a 3rd card
              } # end 3 unique card check
            } # end 3nd card assignment                    
          } # end match 1 card
          if (lval[5] == 14 & lval[1] == 1) {
              lval[1] <- 14  
          } # end group of aces check
          tval <- as.numeric(names(table(lval)[table(lval) > 1]))
          tidx <- grep(tval[1], hand[lidx,]$value)
          # tval > 1 indicates a full house/2pair that can't be avoided
          if (length(tval) > 1) {
            tidx <- c(tidx, grep(tval[2], hand[lidx,]$value))
          } # end 2 pair/full house check         
          lrtn <- hofkind(hand[lidx, ][tidx, ])
          if (length(tidx) < 5) {                  
            lrtn[1] <- as.numeric(lrtn[1]) + sum(lval[-tidx])
            lrtn[2] <- paste(ifelse(hand[lidx, "value"] == "W", 
                         paste(rownames(hand[lidx, ]), "w", sep=""), 
                         rownames(hand[lidx, ])), collapse = " ")
          } # end less than 5 card check
        } # end of "ofkind" hand logic
        #####
        # split back the return values of high hands
        if (!is.null(lrtn)) {
          rval  <- lrtn[1]
          rhand <- lrtn[2]
          rname <- lrtn[3]       
        } # end splitting of return values for high hands
        #####
        # apply hmin logic to enforce minimum hands
        if ((hmin == "8-" & as.numeric(rval) > 421) |
            (hmin == "9-" & as.numeric(rval) > 475) ) {
          rname <- "NoHand"
        } # end apply hmin logic
      } # end 5 wildcard check
    } # end LO hand without enough cards check
  } # end HI vs LO (htyp) check
  ##############
  # Build rhand if it is a vector
  if (rhand[1] != "" & length(rhand) > 1) {
    vhand <- rhand[1]
    for (i in 2:length(rhand)) {
      vhand <- paste(vhand, rhand[i])
    } # end merge hand
    rhand <- vhand
  } # end of building rhand
  c(rval, rhand, rname)
} # End of internal function iwpmakehand

################################################################################





 
