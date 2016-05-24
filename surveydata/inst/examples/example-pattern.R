# Extract the pattern from membersurvey

oldptn <- pattern(membersurvey)
oldptn

# The pattern is used to extract columns

names(membersurvey)
grep("Q20", names(membersurvey), value=TRUE)

head(membersurvey["Q20"])
head(membersurvey["Q20_other"])

# Define a new pattern

pattern(membersurvey) <- list(sep="_", exclude="")
head(membersurvey["Q20"])

# Reset original pattern

pattern(membersurvey) <- oldptn
rm(oldptn)
