str(Aspirin)

# Plot 2 X 2 Contingency tables
plot( ~ treatment + cancer, data = Aspirin)
plot(table(Aspirin), main = "")

# Calculate odds
(Pr.asp <- 18496 / (18496 + 1438))
(Odds.asp <- Pr.asp / (1 - Pr.asp))
(Pr.no.asp <- 18515 / (18515 + 1427))
(Odds.no.asp <- Pr.no.asp / (1 - Pr.no.asp))
(Odds <- Odds.asp / Odds.no.asp)
ln.Odds <- log(Odds)

(SE.Odds <- sqrt(sum(1/table(Aspirin))))
Z <- 1.96
(CI.low <- ln.Odds - Z * SE.Odds)
(CI.high <- ln.Odds + Z * SE.Odds)

exp(CI.low)
exp(CI.high)

# Using oddsRatio() from the mosaic package
# First reformat the data so that "No cancer" is in column 1
# and "Aspirin" is in row 2.
x <- matrix(c(18515, 18496, 1427, 1438), nrow = 2)
x
oddsRatio(x)
