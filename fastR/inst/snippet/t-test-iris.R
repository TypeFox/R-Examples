# for CI; p-value not interesting here
t.test(iris$Sepal.Width[iris$Species=="virginica"])        
# this gives a more interesting p-value
t.test(iris$Sepal.Width[iris$Species=="virginica"],mu=3)   
