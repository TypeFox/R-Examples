.onAttach <- function(...) {

  proton.start = " _____ _          _____         _              _____
|_   _| |_ ___   |  _  |___ ___| |_ ___ ___   |   __|___ _____ ___
  | | |   | -_|  |   __|  _| . |  _| . |   |  |  |  | .'|     | -_|
  |_| |_|_|___|  |__|  |_| |___|_| |___|_|_|  |_____|__,|_|_|_|___|

Your goal is to find Slawomir Pietraszko's credentials for the Proton server.
This is the only way for Bit to find the secret plans of Pietraszko's laboratory. \n
Enter the `proton()` command in order to start the adventure.\n
Remember that at any time you may add `hint=TRUE` argument to the executed command in order to get additional suggestions.
"
   packageStartupMessage(proton.start)
}

dcode <- function(tex) {
  tmp1 <- c(LETTERS, letters)
  tmp2 <- setdiff(unique(unlist(strsplit(tex, split=""))), tmp1)
  let <- c(tmp1, tmp2)
  names(let) <- c(rev(tmp1), tmp2)
  sapply(strsplit(tex, split=""), function(x){
    paste(let[x], collapse="")
  })
}
