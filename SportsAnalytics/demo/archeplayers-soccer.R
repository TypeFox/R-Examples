#' @demo Archetypal soccer players based on skill ratings; analysis
#'   from the manuscript "Archetypal athletes" by Eugster (2011)

library("SportsAnalytics")
library("archetypes")
library("RColorBrewer")

col_pal <- brewer.pal(7, "Set1")
col_black <- rgb(0, 0, 0, 0.2)



### Data: ############################################################

data("EURO4PlayerSkillsSep11")

dat <- subset(EURO4PlayerSkillsSep11,
              Position != "Goalkeeper",
              select = -c(Birthday, Positions))
dat <- subset(dat,  Attack != 6 & TopSpeed > 0)

mat <- as.matrix(subset(dat, select = -c(Team, Name, Number, Nationality,
                                         Age, InjuryTolerance, Foot, Side,
                                         Position, League, KeeperSkills,
                                         Height, Weight, ConditionFitness,
                                         WeakFootAccuracy, WeakFootFrequency)))
rownames(mat) <- NULL


pcplot(mat, col = col_black, las = 2)



### Archetypes: ######################################################

set.seed(1234)
as <- stepArchetypes(mat, k = 1:15)

screeplot(as)

a4 <- bestModel(as[[4]])


### Archetypal soccer players:

parameters(a4)
barplot(a4, mat, percentiles = TRUE)


### Visualization:

pcplot(a4, mat, data.col = col_black, atypes.col = col_pal[1:4])
legend("topleft", legend = sprintf("A%s", 1:4),
       col = col_pal[1:4], lwd = 1, bg = "white")


### Alpha coefficients:

coef <- coef(a4, "alphas")

pcplot(coef, col = c(NA, NA, col_black),
       rx = matrix(c(0, 1), ncol = 4, nrow = 2), var.label = FALSE)


## ... in relation to player position:
pos <- as.character(dat$Position)

cols <- rep("gray", length(pos))
cols[pos == "Defender"] <- col_pal[1]

pcplot(coef, col = c(NA, NA, cols),
       rx = matrix(c(0, 1), ncol = 4, nrow = 2), var.label = FALSE)



### Player interpretation: ###########################################

coef <- coef(a, "alphas")

## The best player is a combination of Archetyp 1 and Archetype 2 with
## Archetype 1 contributing more than Archetype 2:
which <- which(coef[, 3] == 0 & coef[, 4] == 0 &
               coef[, 1] > 0 & coef[, 2] > 0 &
               coef[, 1] > coef[, 2])

cbind(dat[which, c("Name", "Team")],
      coef[which, ])

