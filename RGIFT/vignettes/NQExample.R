sink(file="NQExample.txt", type="output")

GIFTcomment("Examples of Numerical Questions")

GIFTcategory("NQuestions")

#Question 1
GIFTNQ("What's the mean of vector c(.4, .4, .5, .3)",
   as.character(mean(c(.4, .4, .5, .3))), .01)

#Question 2
GIFTNQ("What's the mean of vector c(.4, .4, .5, .3)",
   c(0.39, 0.41))

sink()

