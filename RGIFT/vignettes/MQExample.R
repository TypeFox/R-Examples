sink(file="MQExample.txt", type="output")

GIFTcomment("Examples of Matching Questions")

GIFTcategory("MQuestions")

#Question 1
GIFTM("Match the following operations to their respective R commands:",
   c("mean", "variance", "standard deviation"), c("mean()", "var()", "sd()"))

#Question 2
GIFTM("Match each vector with its mean:",
   c("c(1,2,3)", "c(4,5,6)", "c(7,8,9)"), 
   as.character(c(mean(1:3), mean(4:6), mean(7:9))) )

sink()

