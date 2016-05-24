sink(file="MCQExample.txt", type="output")

GIFTcomment("Examples of Multiple Choice Questions")

GIFTcategory("MCQuestions")

#Question 1
GIFTMC("What's the mean of 1, 2, and 3?", c("1", "2", "3"), rightans=2, 
   wwrong="-33.333")

#Question 2
GIFTMC("What's the command to compute the mean?", 
   c("mean()", "sd()", "var()", "length()"), wwrong="-33.333")

#Question 3; Good answer plus a not so good answer
GIFTMC("How can we compute the mean of 1, 2, and 3?", 
   c("mean(c(1,2,3))", "sd(c(1,2,3))", "(1+2+3)/3"), wright=c("100", "75"), 
   rightans=c(1,3),  wwrong="-33.333")

#Question 4; Multiple valid answers
GIFTMC("How can we compute the mean of 1, 2, and 3?", 
   c("mean(c(1,2,3))", "sd(c(1,2,3)", "(1+2+3)/3"), wright=c("50", "50"), 
   rightans=c(1,3),  wwrong="-33.333")

sink()

