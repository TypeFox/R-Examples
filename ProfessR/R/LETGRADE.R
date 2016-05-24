`LETGRADE` <-
function(g)
  {

    lett = rep("I", length=length(g))
    
    SCRS = seq(from=100, by=(-4), length=13)
    LETS = c("A+", "A", "A-", "B+", "B", "B-", "C+", "C", "C-", "D+", "D", "D-", "E", "E")
    SCRS[1] = 100.1
    SCRS[length(SCRS)+1]  = 0

    for(i in 2:length(SCRS))
      {
        lett[g>=SCRS[i] & g<SCRS[i-1] ] = LETS[i]

      }
    

    return(lett)

  }

