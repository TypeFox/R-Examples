### R code from vignette source 'ExamMaker.Rnw'

###################################################
### code chunk number 1: ExamMaker.Rnw:132-136
###################################################
library("ProfessR")

QB1 = Get.testbank('climate_change_abbott.txt') 
QB2 = Get.testbank('tsunami_abbot_web.txt') 


###################################################
### code chunk number 2: ExamMaker.Rnw:141-151
###################################################

L1 = length(QB1)
L2 = length(QB2)

isel1 = sample(1:L1)
QBfinal =  list()
for(i in 1:12) { QBfinal[[i]] = QB1[[isel1[i]]] }
isel2 = sample(1:L2)
for(i in 1:9) { QBfinal[[i+12]] = QB2[[isel2[i]]] }



###################################################
### code chunk number 3: ExamMaker.Rnw:156-157
###################################################
QA1 = ran.exam(QBfinal)


###################################################
### code chunk number 4: ExamMaker.Rnw:163-169
###################################################
examdate="TUES OCT 30 2007"

version.exam(QBfinal, "exam2A" , exnumber="Exam 2", seqnum="1", examdate=examdate)

version.exam(QBfinal, "exam2B" , exnumber="Exam 2", seqnum="2", examdate=examdate)



