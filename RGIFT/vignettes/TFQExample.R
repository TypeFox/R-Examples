sink(file="TFQExample.txt", type="output")

GIFTcomment("Examples of True-False Questions")

GIFTcategory("TFQuestions")

#Question 1
GIFTTF("The mean of 1, 2, and 3 is 3?", TRUE)

#Question 2
GIFTTF("The command to compute the mean is sd()", FALSE)

#Question 3
GIFTTF("R is the best programming language ever!!", TRUE)

sink()

