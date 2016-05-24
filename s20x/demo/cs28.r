data (course.df)
stripchart(Exam~Stage1,pch=1,method="stack",main="Exam mark by Stage 1",ylab="grade in 10x",xlab="20x exam mark",data=course.df)
summaryStats(Exam~Stage1,data=course.df)
exam.fit<-lm(Exam~Stage1,data=course.df)
eovcheck(exam.fit)
normcheck(exam.fit)
summary1way(exam.fit)
multipleComp(exam.fit)

