
names(membersurvey)
head(membersurvey["Q1"])
head(membersurvey[c("Q1", "Q2")])
head(membersurvey[membersurvey$Q2=="2009", c("Q1", "Q2")])

# The pattern is used to extract columns

pattern(membersurvey)

grep("Q20", names(membersurvey), value=TRUE)
head(membersurvey["Q20"])
head(membersurvey["Q20_other"])

