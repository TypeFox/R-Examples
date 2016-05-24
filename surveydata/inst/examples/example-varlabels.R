# Extract the variable labels from membersurvey

ms <- membersurvey[, c("id", "Q1", "Q2")]

str(ms)
varlabels(ms)
varlabels(ms)["Q2"]

# Assign a new value to the text of question 2

varlabels(ms)["Q2"] <- "When did you join?"
varlabels(ms)
str(ms["Q2"])

