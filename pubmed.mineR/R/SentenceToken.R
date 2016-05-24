SentenceToken = function (x)
   {
        temp2 = gregexpr("[abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890():][.][ ][ABCDEFGHIJKLMNOPQRSTUVWXYZ]",
                                            x)
         token = NULL
        if (temp2[[1]][1] == -1)
             token = x
         else if (length(temp2[[1]]) == 1) {
               tempa = 1
               tempb = temp2[[1]][1] + 1
               token = c(token, substr(x, start = tempa, stop = tempb))
               tempa = tempb + 1
               tempb = nchar(x)
               token = c(token, substr(x, start = tempa, stop = tempb))
           }
         else if (length(temp2[[1]]) == 2) {
               for (i in 1:length(temp2[[1]])) {
                     beginning = i == 1
                     between = (i > 1) & (i <= length(temp2[[1]]))
                     if (beginning) {
                           tempa = 1
                           tempb = temp2[[1]][i] + 1
                           token = c(token, substr(x, start = tempa, stop = tempb))
                       }
                     else if (between) {
                           tempa = temp2[[1]][i - 1] + 2
                           tempb = temp2[[1]][i] + 1
                           token = c(token, substr(x, start = tempa, stop = tempb))
                           tempa = tempb + 1
                           tempb = nchar(x)
                           token = c(token, substr(x, start = tempa, stop = tempb))
                       }
                 }
           }
         else for (i in 1:length(temp2[[1]])) {
               beginning = i == 1
               between = (i > 1) & (i <= length(temp2[[1]]))
               end = i == length(temp2[[1]])
               if (beginning) {
                     tempa = 1
                     tempb = temp2[[1]][i] + 1
                     token = c(token, substr(x, start = tempa, stop = tempb))
                 }
               else if (between) {
                     tempa = temp2[[1]][i - 1] + 2
                     tempb = temp2[[1]][i] + 1
                     token = c(token, substr(x, start = tempa, stop = tempb))
                 }
               if (end) {
                     tempa = temp2[[1]][i] + 2
                     tempb = nchar(x)
                     token = c(token, substr(x, start = tempa, stop = tempb))
                 }
           }
         return(token)
     }
