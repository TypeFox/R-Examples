#Question 1

id = c(1:6)
age = c(30, 32,28,39,20,25)
edu = c(0,0,0,0,0,0)
class1 = c("poor","poor","poor","middle","middle","middle")
class2 = c("middle","poor")

factor(class1, levels=class2)

#Question 2
indianMother = data.frame(id = id,
                            age=age,
                            edu=edu,
                            class=class1)
