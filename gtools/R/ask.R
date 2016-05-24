ask <- function(msg="Press <RETURN> to continue: ")
        {
          cat(msg);
          readLines(con=stdin(),n=1)
        }
