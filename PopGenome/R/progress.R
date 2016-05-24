################################################
# Functions for a progress bar and similar stuff
#
#  progressDots: Prints one dot per x loops to
#                the standard output.
#     Arguments: x     the number of loops per dot
#                i     the current value of the
#                      loop variable
#                total the maximum nuber of loops
#         Usage: For example in a for loop:
#                for (i in 1:n) {
#                    progressDots(x = 10, i, total = n)
#                    # something time consuming
#                }

#  progressBar:  Prints a progress bar to the
#                standard output.
#    Arguments:  i          the current value of the
#                           loop variable
#                total      the maximum nuber of loops
#                symb_drawn the number of printed
#                           symbols for the bar
#        Usage:  for example in a for loop:
#                progr <- progressBar()
#                for (i in 1:n) {
#                    progr <- progressBar(i, total = n, symb_drawn = progr)
#                    # something time consuming
#                }

progressDots <- function(x, i, total, symb = ".", num = 80, nl = "\n") {

  # x:     one dot per x loops
  # i:     loop variable
  # total: maximum value of i
  # symb:  printed symbol
  # num:   maximum number of dots per line
  # nl:    newline character

  if (i %% x == 0) {              # print symb after x loops
      cat(symb)
      if (i %% (x * num) == 0) {  # print a newline after num symbs
          cat(nl)
      }
  }
  if (i == total) {               # print a newline after the last symb
      cat(nl)
  }

}


progressBar <- function (i, total, symb_drawn = 0, symb = "=", nl = "\n") {

  # i:          loop variable
  # total:      maximum value of i
  # symb_drawn: number of printet symbs
  # symb:       printed symbol
  # nl:         newline character

  if (missing(i) && missing(total)) {
      init = TRUE
  }
  else {
      init = FALSE
  }

  if (init == TRUE) {
      # print a timeline
      cat("|            :            |            :            | 100 %", nl, "|", sep = "")
  }
  else {

      if (symb_drawn < 51) {

          progress <- round((i/total) * 52, digits = 0) - symb_drawn
          symb_drawn <- progress + symb_drawn

          while (progress > 0) {
              cat(symb);
              progress <- progress - 1
          }

          #return(symb_drawn)
      }
      else {
          if (i == total) {
              cat("| ;-) \n")
          }
      }
  }
  
  return(symb_drawn)
}
