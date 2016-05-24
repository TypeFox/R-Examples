## ------------------------------------------------------------------------
shakespeare <- paste(
  "Haste still pays haste, and leisure answers leisure;",
  "Like doth quit like, and MEASURE still FOR MEASURE.",
  "Then, Angelo, thy fault's thus manifested;",
  "Which, though thou wouldst deny, denies thee vantage.",
  "We do condemn thee to the very block",
  "Where Claudio stoop'd to death, and with like haste.",
  "Away with him!")
critic <- paste(
  "The play comes to its culmination where Duke Vincentio, quoting from",
  "the words of the Sermon on the Mount, says,",
  "'Haste still goes very quickly , and leisure answers leisure;",
  "Like doth cancel like, and measure still for measure.'",
  "These titular words sum up the meaning of the play.")

## ------------------------------------------------------------------------
library(textreuse)
align_local(shakespeare, critic)

