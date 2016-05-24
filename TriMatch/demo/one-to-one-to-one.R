require(TriMatch)
data(tutoring)

formu <- ~ Gender + Ethnicity + Military + ESL + EdMother + EdFather + Age +
	Employment + Income + Transfer + GPA

table(tutoring$treat, useNA='ifany')

tutoring.tpsa <- trips(tutoring, tutoring$treat, formu)

tutoring.matched.1to1 <- trimatch(tutoring.tpsa, exact=tutoring[,c("Course")], 
								  method=OneToN, M1=1, M2=1)

tutoring.matched.1to1 <- tutoring.matched.1to1[!duplicated(tutoring.matched.1to1[,1]),]
tutoring.matched.1to1 <- tutoring.matched.1to1[!duplicated(tutoring.matched.1to1[,2]),]
tutoring.matched.1to1 <- tutoring.matched.1to1[!duplicated(tutoring.matched.1to1[,3]),]

nrow(tutoring.matched.1to1)
length(unique(tutoring.matched.1to1$Treat1))
length(unique(tutoring.matched.1to1$Treat2))
length(unique(tutoring.matched.1to1$Control))

summary(tutoring.matched.1to1, tutoring$Grade)

plot(tutoring.matched.1to1)
boxdiff.plot(tutoring.matched.1to1, tutoring$Grade)
