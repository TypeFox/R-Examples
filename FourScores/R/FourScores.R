FourScores <-
function(modus = "normal", rows = 6, collumns = 7){


  ###Here are some help-functions
  
  
  ##help-function which gets and returns the players' names
  getPlayerNames <- function(){
    for(i in 1:length(PlayerNames)){
      plot(1, axes = FALSE, ann = FALSE, col = "white")
      mtext(paste(PlayerNames[i], ", please enter your name:", sep = ""))
  
      #initializing a variable needed for the generation of the names
      SingleLetters <- ""
  
      #boolean variable for checking if enter was given
      enter <- FALSE
      while(!enter){
      
        #getting a single letter, bye hitting a key on the keybord
        tmp <- getGraphicsEvent(onKeybd = typing)
       
        #if the player presses [Enter] (encoded by "ctrl-J")...     
        if(tmp == "ctrl-J"){
        
          #... the boolean gets true
          enter = TRUE
      
        #im the case, the player presses [Backspace] (encoded by "ctrl-H")...
        }else if(tmp == "ctrl-H"){
      
          #... the last saved single letter, is removed from the singleletters
          SingleLetters <- substr(SingleLetters, 1, nchar(SingleLetters) -1)
        
        #in all other cases... ###INCLUDING [PageUp], [Home] etc. SO SHOULD USE tmp %in% letters()
        }else{
      
          #... the temporary variable is pasted onto the existing singleletters
          SingleLetters <- paste(SingleLetters, tmp, sep = "")
        } 
      
        #"empty" plot (white dot on white background) 
        plot(1, axes = FALSE, ann = FALSE, col = "white")
      
        #Instruction to the player
        mtext(paste(PlayerNames[i], ", please enter your name:", sep = ""))
      
        #displaying the progress (that means the current -unconfirmed- singleletters)
        text(x = 1, y = 1, SingleLetters) 
      }
    
      #overwriting the default names, by the input of the players
      PlayerNames[i] <- SingleLetters
    }
  
    #returning the names
    return(PlayerNames)
  }
  
  
  ##help-function which generates the playing-field
  FieldGeneration <- function(){

    #the (at the generation empty) matrix will contain e.g. on line 3 collumn 5 the number 2, 
    #if player 2 was the one how placed the third time a stone in collumn 3
    field <- matrix(nrow = rows, ncol = collumns)
    return(field)
  }

  ##help-function which "throws" the stone into the field and returns the new field
  NewField <- function(){

    #the number of the row, where the stone "lands" in the selected collumn,
    #is equal to the number of free rows in that coullum (when the highest row ist labeled with "1", the second highest with "2" and so on
    roww <- sum(is.na(field[, collumn]))  
  
    #if there is at least one free row
    if(roww >= 1){
    
      #the stone with the number of the player on it, 
      #"falls" in the selected collumn on the calculated row
      field[roww, collumn] <- player
    }

    return(field)
  }



  ##help-function which returns, the key on the keybord which is beeing typed
  typing <- function(key){
    return(key)
  }
  
    ##help-function which return the x-axis-value of the mouse when releasing the mousebutton
  clicking <- function(buttons, x, y) {

    #the x-axis-value of the Mouse
    plx <- grconvertX(x, from = "ndc", to = "user")
  
    #rounds the x-axis-value to a natural number and returns it
    return(round(plx))        
  }      


  ##help-function which gives an optical help for the player
  ##he sees which collumn is selected if he releases the mouse button
  preview <- function(buttons, x, y) {
                 
    #the x-axis-value of the mouse
    plx <- grconvertX(x, "ndc", "user")           
  
    #rounds the x-axis-value to a natural number
    collumn <- round(plx)
  
    #getting the heigth of the field and add 0.5 Points, so that the preview is above and not over the actual field 
    displayheigth <- nrow(field) + 0.5
  
    #variables that has to be used in every plotting, so that the old grey previews "disapear" (are overplottet by white previews)
    xvalues <- 1:ncol(field)              
    yvalues <- rep(displayheigth, times = ncol(field))

    color <- rep("white", times = ncol(field))
    
    ##where the mouse is, the color shall be grey
    color[collumn] <- "grey68" 
  
    #Vector for the symbols which destinguishes the player's stones (in addition to the color - which is unwanted in the preview)
    PlayerPch <- c("O", "X")
  
    #plot a grey preview
    points(x = xvalues, y = yvalues, cex = 1.5, pch = PlayerPch[player], lty = 3, col = color)#, col = farben[spieler + 1])
  }





  ##a major-function which plots the current field, and if given a hint, which player has won
  FieldPlot <- function(){
  
    #vector of the colors
    color <- c("green", "blue")
  
    #vector of the symbols (seen in the title obove the field)
    PlayerPch <- c("O", "X")
  
    par(las = 1)
  
    #plotting the field, dependent on the number of rows and collumns
    plot(c(1, ncol(field)), c(1, nrow(field) + 0.5), type = "n", xlab = "", ylab = "",
      main = paste(PlayerNames[player], ": it's your turn", sep = ""))

    #plotting the stones of each row into the field
    for(i in 1:nrow(field)){
  
      #on the positions x=1,2,..7 and y=7,7,...7, than 6,6,...6 
      points(x = 1:ncol(field), y = rep(nrow(field) - i + 1, times = ncol(field)), 
        col = color[field[i, ]],   #dependent on which player has placed which stone, their color is picked from the color-vector
        cex = 2, pch = PlayerPch[field[i, ]], lwd = 2)  #and which symbol they refer to
    }
  
  }




  ##help-function that checks whether the field is correct
  FieldCorrect <- function(){
  
    #turning warning messages off
    options(warn = -1)
    
    #turning warning messages on again
    options(warn = 0) 
    
    #transforming a typed "7" into a 7
    collumn <- as.integer(collumn)
    roww <- 0
    RangeCorrect <- (collumn >= 1 && collumn <= ncol(field)) 
    if(!is.na(RangeCorrect) && RangeCorrect){
      roww <- sum(is.na(field[, collumn]))
    }
    return(roww >= 1 && RangeCorrect)
  }




  ##help-function that checks whether (at least) one of the for possibilities of winning is given
  FieldWinCheck <- function(){

    #initializing boolean variable that will get TRUE if (at least) one of the for possibilities of winning is given
    won <- FALSE
  
    #check on four in a vertical row
    for(j in 1:ncol(field)){
      for(i in 1:(nrow(field) - 3)){
        if(sum(field[i:(i+3), j] == player, na.rm = TRUE) >= 4){
          won <- TRUE
        }
      }
    }
  
    #check on four in a horizontal row
    for(i in 1:nrow(field)){
      for(j in 1:(ncol(field) - 3)){
        if(sum(field[i, j:(j+3)] == player, na.rm = TRUE) >= 4){
          won <- TRUE
        }
      }
    }
  
    #check on four in a diagonal row(slope = +1)
    for(i in 1:(nrow(field) - 3)){
      for(j in 1:(ncol(field) - 3)){
        if(sum(field[i:(i+3), j:(j+3)][c(4, 7, 10, 13)] == player, na.rm = TRUE) >= 4){
          won <- TRUE
        }
      }
    }
  
    #check on four in a diagonal row(slope = -1)
    for(i in 1:(nrow(field)-3)){
      for(j in 1:(ncol(field) - 3)){
        if(sum(diag(field[i:(i+3), j:(j+3)]) == player, na.rm = TRUE) >= 4){
          won <- TRUE
        }
      }
    }
  
    return(won)
  }



  ###Here starts the main part of the game


  #Entering the names
  PlayerNames <- c("Player1", "Player2")
  if(modus != "quick"){
    PlayerNames <- getPlayerNames()
  }
  
  #generating the playing field (per default 6 collumns and 7 rows)
  field <- FieldGeneration()
  
  #initializing a boolean variable which indicates if a player has won
  won <- FALSE
  
  #initializing a boolean variable which indicates if it is a draw
  draw <- FALSE
  
  #initializing a variable to "keep in mind" which player has to move
  i <- 1
  
  #while nobody has won, the game goes on and a next round is started
  while(!won && !draw){
  
    #increasing the number of rounds
    i <- i + 1 
    
    #Variable for remebering whether player1 oder player2 has to move
    player <- (i %% 2) + 1
    
    #initalizing a boolean variable which indicates wheter a valid input is given by the player
    correct <- FALSE
    
    #while no valid input is given, the player is said to choose a valid collumn
    while(!correct){
    
      #...if a correct input is given, the field is plotted
      FieldPlot()
    
      #get an input for the collum by using 'getGraphicsEvent()'
      collumn <- as.integer(getGraphicsEvent(
        paste(PlayerNames[player], ": select collumn", sep = ""),    #appears in the console, can be removed when the text appears in the plot
        onKeybd = typing,
        onMouseUp = clicking, 
        onMouseMove = preview))
              
      #checking by a help-function whether the choosen cullumn an the current field fits to gether
      correct <- FieldCorrect()       
    }
    
    #generating a new field by adding to the old field the information about the selected collumn and which player has chossen this collumn
    field <- NewField()
    
    #checken by a help-function whether the current player has won by this move
    won <- FieldWinCheck()
    
    #checken if there is a draw
    draw <- !won && ((i - 1) == (ncol(field) * nrow(field)))
  }
  
  #plotting the current field  ###MAYBE WITH THE HIGHLITHED WINNING ROW
  # FieldPlot()
  
  plot(1, axes = FALSE, ann = FALSE, col = "white")
  if(won){
    #giving a message to the winner 
    msg <- paste(PlayerNames[player], " has won!", sep = "")   
    text(x = 1, y = 1, msg, cex = 3, col = "lightgoldenrod3") 
  }else{
    msg <- "It's a draw"      
    text(x = 1, y = 1, msg, cex = 3, col = "slategray4")   
  }
 
}

