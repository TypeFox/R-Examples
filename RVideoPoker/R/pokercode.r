rvideopoker <- function(initbet=1, blnce=1000) {

  state <- new.env()
  
  
  samplespaceJACKS <- paste(sort(rep(c("c","d","h","s"),13)),
                            rep(c(1:10,"j","q","k"),4),sep="")
  payoutscheme <- structure(list(Hand = structure(c(5L, 4L, 10L, 9L, 7L, 1L, 3L, 
                                   2L, 8L, 6L),
                                   .Label = c("Flush", "Four of a Kind",
                                     "Full House",  "Jacks or Better",
                                     "Nothing", "Royal Flush",
                                     "Straight", "Straight Flush", "Three of a Kind", "Two Pairs"),
                                   class = "factor"), Returnvalue = 0:9,  
                                 Credits1 = c(0L, 1L, 2L, 3L, 4L, 6L, 9L, 25L, 50L, 250L), 
                                 Credits2 = c(0L, 2L, 4L, 6L, 8L, 12L, 18L, 50L, 100L, 500L),
                                 Credits3 = c(0L, 3L, 6L, 9L, 12L, 18L, 27L, 75L, 150L, 750L),
                                 Credits4 = c(0L, 4L, 8L, 12L, 16L, 24L, 36L, 100L, 200L, 1000L),
                                 Credits5 = c(0L, 5L, 10L, 15L, 20L, 30L, 45L, 125L, 250L, 4000L)),
                            .Names = c("Hand", "Returnvalue", "Credits1", "Credits2", "Credits3", "Credits4", "Credits5"),
                            class = "data.frame", row.names = c(NA, -10L))
  
  plotFrame <- function() {
    par(mai=c(0,0,0.5,0))
    plot(1,1, type="n", xlim=c(0,100), ylim=c(0,20), axes=FALSE, xlab="",
         ylab="")
  }
  
  drawCardsInitially <- function() {
    sample(x=samplespaceJACKS, size=5)
  }
  
  
  showBack <- function() {
    require(pixmap)
    card1 <- card2 <- card3 <- card4 <- card5 <- read.pnm(system.file("pictures/b1fv.pnm", package="RVideoPoker")[1])
    plotFrame()
    addlogo(x=card1, px=c(2,18), py=c(2,18))
    addlogo(x=card2, px=c(22,38), py=c(2,18))
    addlogo(x=card3, px=c(42,58), py=c(2,18))
    addlogo(x=card4, px=c(62,78), py=c(2,18))
    addlogo(x=card5, px=c(82,98), py=c(2,18))
    title(main=paste("Your Balance: ", blnce, sep=""))
  }


  initcards <- function(panel) {
    state$ahand <- drawCardsInitially()
    showCardsInitially()
    panel
  }


  showCardsInitially <- function() {
    plotFrame()
    card1 <- read.pnm(system.file(paste("pictures/", state$ahand[1], ".pnm", sep=""), package="RVideoPoker")[1])
    card2 <- read.pnm(system.file(paste("pictures/", state$ahand[2], ".pnm", sep=""), package="RVideoPoker")[1])
    card3 <- read.pnm(system.file(paste("pictures/", state$ahand[3], ".pnm", sep=""), package="RVideoPoker")[1])
    card4 <- read.pnm(system.file(paste("pictures/", state$ahand[4], ".pnm", sep=""), package="RVideoPoker")[1])
    card5 <- read.pnm(system.file(paste("pictures/", state$ahand[5], ".pnm", sep=""), package="RVideoPoker")[1])
    addlogo(x=card1, px=c(2,18), py=c(2,18))
    addlogo(x=card2, px=c(22,38), py=c(2,18))
    addlogo(x=card3, px=c(42,58), py=c(2,18))
    addlogo(x=card4, px=c(62,78), py=c(2,18))
    addlogo(x=card5, px=c(82,98), py=c(2,18))
    ante <- (initbet*state$initcredits)
    state$blnce <- blnce-ante
    if (blnce < 0) {
      rp.messagebox("Game Over!\n Game Will Start Again at Initial Values.", title="")
      state$blnce <- state$origbalance
    }
    title(main=paste("Your Balance: ", state$blnce, sep=""))
  }
  
  setbet <- function(panel) {
    state$initbet <- as.numeric(state$penny.nickel.quarter.dollar)
    panel
  }
  
  
  creditset <- function(panel) {
    state$initcredits <- as.numeric(state$no.of.credits)
    panel
  }
  
  state$whichchange <- rep(TRUE, 5)
  
  
  card1check <- function(panel) {
    state$whichchange[1] <- !(state$keepC1)
    panel
  }
  card2check <- function(panel) {
    state$whichchange[2] <- !(state$keepC2)
    panel
  }
  card3check <- function(panel) {
    state$whichchange[3] <- !(state$keepC3)
    panel
  }
  card4check <- function(panel) {
    state$whichchange[4] <- !(state$keepC4)
    panel
  }
  card5check <- function(panel) {
    state$whichchange[5] <- !(state$keepC5)
    panel
  }


  showCardsCheck <- function(panel) {
    plotFrame()
    card1 <- read.pnm(system.file(paste("pictures/", state$ahand[1], ".pnm", sep=""), package="RVideoPoker")[1])
    card2 <- read.pnm(system.file(paste("pictures/", state$ahand[2], ".pnm", sep=""), package="RVideoPoker")[1])
    card3 <- read.pnm(system.file(paste("pictures/", state$ahand[3], ".pnm", sep=""), package="RVideoPoker")[1])
    card4 <- read.pnm(system.file(paste("pictures/", state$ahand[4], ".pnm", sep=""), package="RVideoPoker")[1])
    card5 <- read.pnm(system.file(paste("pictures/", state$ahand[5], ".pnm", sep=""), package="RVideoPoker")[1])
    if (state$whichchange[1]) {
      card1 <- read.pnm(system.file("pictures/b1fv.pnm", package="RVideoPoker")[1])
    }
    if (state$whichchange[2]) {
      card2 <- read.pnm(system.file("pictures/b1fv.pnm", package="RVideoPoker")[1])      
    }
    if (state$whichchange[3]) {
      card3 <- read.pnm(system.file("pictures/b1fv.pnm", package="RVideoPoker")[1])      
    }
    if (state$whichchange[4]) {
      card4 <- read.pnm(system.file("pictures/b1fv.pnm", package="RVideoPoker")[1])      
    }
    if (state$whichchange[5]) {
      card5 <- read.pnm(system.file("pictures/b1fv.pnm", package="RVideoPoker")[1])      
    }
    
    addlogo(x=card1, px=c(2,18), py=c(2,18))
    addlogo(x=card2, px=c(22,38), py=c(2,18))
    addlogo(x=card3, px=c(42,58), py=c(2,18))
    addlogo(x=card4, px=c(62,78), py=c(2,18))
    addlogo(x=card5, px=c(82,98), py=c(2,18))
    panel
  }

  evalhand <- function(panel) {
    suits <- substr(state$ahand2,1,1)
    values <- substr(state$ahand2,2,3)
    returnvalue <- 0
    if(is.jacks.or.better(values=values)) returnvalue <- 1  
    if(is.2.pairs(values=values)) returnvalue <- 2  
    if(is.3.of.a.kind(values=values)) returnvalue <- 3  
    if(is.straight(values=values)) returnvalue <- 4  
    if(is.flush(suits=suits)) returnvalue <- 5  
    if(is.full.house(values=values)) returnvalue <- 6  
    if(is.4.of.a.kind(values=values)) returnvalue <- 7  
    if(is.straight.flush(suits=suits, values=values)) returnvalue <- 8  
    if(is.royal.flush(suits=suits, values=values)) returnvalue <- 9
    meaning <- rev(c("Royal Flush", "Straight Flush", "4 of a kind",
                     "Full House", "Flush", "Straight", "3 of a kind",
                     "Two Pairs", "Jacks or Better"))
    whichcol <- 2 + state$initcredits
    yourwin <- initbet * payoutscheme[returnvalue+1,whichcol]
    state$blnce <- blnce + yourwin
    if (returnvalue %in% 1:9) {
      state$returnmessage <- (paste("You've got: ", meaning[returnvalue],
                               ". You win: ", yourwin,
                               ". Your new balance: ", blnce, sep=""))
    } else {
      state$returnmessage <- (paste("You've got nothing! Your new balance: ",
                               blnce, sep=""))
    }
    panel
  }


  drawNewHand <- function(panel) {
    newsamplespace <- samplespaceJACKS[!(samplespaceJACKS %in% state$ahand)]
    newsample <- sample(x=newsamplespace, size=sum(state$whichchange))
    ahand2internal <- state$ahand
    ahand2internal[state$whichchange] <- newsample
    state$ahand2 <- ahand2internal
    plotFrame()
    card1 <- read.pnm(system.file(paste("pictures/", state$ahand2[1], ".pnm", sep=""), package="RVideoPoker")[1])
    card2 <- read.pnm(system.file(paste("pictures/", state$ahand2[2], ".pnm", sep=""), package="RVideoPoker")[1])
    card3 <- read.pnm(system.file(paste("pictures/", state$ahand2[3], ".pnm", sep=""), package="RVideoPoker")[1])
    card4 <- read.pnm(system.file(paste("pictures/", state$ahand2[4], ".pnm", sep=""), package="RVideoPoker")[1])
    card5 <- read.pnm(system.file(paste("pictures/", state$ahand2[5], ".pnm", sep=""), package="RVideoPoker")[1])
    addlogo(x=card1, px=c(2,18), py=c(2,18))
    addlogo(x=card2, px=c(22,38), py=c(2,18))
    addlogo(x=card3, px=c(42,58), py=c(2,18))
    addlogo(x=card4, px=c(62,78), py=c(2,18))
    addlogo(x=card5, px=c(82,98), py=c(2,18))
    evalhand(panel)
    rp.messagebox(state$returnmessage, title="The results are:")
    panel
  }
  
  is.royal.flush <- function(suits, values) {
    is.straight.flush(values=values, suits=suits) & all(c("10", "j", "q","k", "1") %in% values)
  }
  
  is.straight.flush <- function(suits, values) {
    is.straight(values) & is.flush(suits)
  }

  is.4.of.a.kind <- function(values) {
    return(any(table(values)==4))
  }
  

  is.full.house <- function(values) {
    the.table <- table(values)
    all(c(3,2) %in% the.table)
  }
  
  
  is.flush <- function(suits) {
    any(table(suits)==5)
  }

  is.straight <- function(values) {
    newvalues1 <- values
    newvalues1[values=="j"] <- 11
    newvalues1[values=="q"] <- 12
    newvalues1[values=="k"] <- 13
    newvalues1[values=="1"] <- 14
    newvalues1 <- sort(as.numeric(newvalues1))
    returnvalue1 <- all(diff(newvalues1)==1)
    newvalues2 <- values
    newvalues2[values=="j"] <- 11
    newvalues2[values=="q"] <- 12
    newvalues2[values=="k"] <- 13
    newvalues2[values=="1"] <- 1
    newvalues2 <- sort(as.numeric(newvalues2))
    returnvalue2 <- all(diff(newvalues2)==1)
    return(any(c(returnvalue1, returnvalue2)))
  }

  is.3.of.a.kind <- function(values) {
    return(any(table(values)==3))
  }
  
  is.2.pairs <- function(values) {
    the.table <- table(values)
    sum(the.table == 2)==2
  }

  is.jacks.or.better <- function(values) {
    newvalues <- values[values %in% c("j","q","k","1")]
    return(any(table(newvalues)==2))
  }

  state$initbet <- initbet
  state$blnce <- blnce
  state$origbalance <- blnce
  showBack()
  mypanel <- rp.control(title="Video Poker in R")
  rp.radiogroup(mypanel, var=state$penny.nickel.quarter.dollar, initval=1,
              values=c(1,5,25,100), title="Bet (in Cents)", action=setbet)
  rp.slider(mypanel, var=state$no.of.credits, from=1, to=5, initval=1,
            title="How many Credits?", resolution=1, action=creditset,
            showvalue=TRUE) 
  rp.button(mypanel, action=initcards, title="Start Game")
  rp.checkbox(panel=mypanel, var=state$keepC1, action=card1check, title="Keep Card 1?")
  rp.checkbox(panel=mypanel, var=state$keepC2, action=card2check, title="Keep Card 2?")
  rp.checkbox(panel=mypanel, var=state$keepC3, action=card3check, title="Keep Card 3?")
  rp.checkbox(panel=mypanel, var=state$keepC4, action=card4check, title="Keep Card 4?")
  rp.checkbox(panel=mypanel, var=state$keepC5, action=card5check, title="Keep Card 5?")
  
  rp.button(mypanel, action=showCardsCheck, title="Check Selection")
  rp.button(mypanel, action=drawNewHand, title="Draw Again")
  rp.block(mypanel)
}
