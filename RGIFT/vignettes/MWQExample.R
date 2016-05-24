sink(file="MWQExample.txt", type="output")

GIFTcomment("Examples of Missing Word Questions")

GIFTcategory("MWQuestions")

#Question 1
GIFTMW("With command ", " the mean of a vector of values can be computed", 
c("mean()", "sd()", "var()"), rightans=1)

#Question 2
GIFTMW("With command ", " the variance of a vector of values can be computed", 
c("mean()", "sd()", "var()"), rightans=3)

#Question 3
GIFTMW("With command ", " the standard deviation  of a vector of values can be computed", 
c("mean()", "sd()", "var()"), rightans=2)

sink()

