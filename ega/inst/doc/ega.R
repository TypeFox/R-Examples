## ------------------------------------------------------------------------
library(ega)
library(ggplot2)

ggplot(glucose_data, aes(ref, test)) + geom_point()

## ------------------------------------------------------------------------
cor(glucose_data$ref, glucose_data$test)

## ------------------------------------------------------------------------
zones <- getClarkeZones(glucose_data$ref, glucose_data$test)

head(zones)

## ------------------------------------------------------------------------
zones <- factor(zones)

# counts
table(zones)

# percentages
table(zones)/length(zones)*100

## ------------------------------------------------------------------------
plotClarkeGrid(glucose_data$ref, glucose_data$test)

## ------------------------------------------------------------------------
plotClarkeGrid(glucose_data$ref, glucose_data$test, 
               pointsize=1.5, 
               pointalpha=0.6, 
               linetype="dashed")

## ------------------------------------------------------------------------
ceg <- plotClarkeGrid(glucose_data$ref, glucose_data$test)

ceg + theme_gray() + 
  theme(plot.title = element_text(size = rel(2), colour = "blue"))

## ------------------------------------------------------------------------
zones <- getParkesZones(glucose_data$ref, glucose_data$test)

zones <- factor(zones)

# counts
table(zones)

# percentages
table(zones)/length(zones)*100

## ------------------------------------------------------------------------
plotParkesGrid(glucose_data$ref, glucose_data$test)

