data (course.df)
stripchart(Exam~Attend,pch=1,method="stack", ylab="regularly attend",xlab="exam mark",main="Exam by Attend",data=course.df)
summaryStats(Exam~Attend,data=course.df)
exam.fit<-lm(Exam~Attend,data=course.df)
eovcheck(exam.fit)
normcheck(exam.fit)
t.test(Exam~Attend,var.equal=T,data=course.df)

