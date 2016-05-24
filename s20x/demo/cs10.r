data (course.df)
diffs<-course.df$Assign-course.df$Test
stripchart(diffs,pch=1,method="stack",main="Assignment score - Test score",xlab="differences")
summaryStats(diffs)
normcheck(diffs)
t.test(course.df$Assign,course.df$Test,paired=T)

