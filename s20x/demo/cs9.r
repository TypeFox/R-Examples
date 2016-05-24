data (course.df)
course.df$Exam
stripchart(course.df$Exam,method="stack",pch=1,main="STATS 20x Exam marks",xlab="exam mark")
summaryStats(course.df$Exam)
normcheck(course.df$Exam)
t.test(course.df$Exam)

