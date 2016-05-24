for (species in levels(iris$Species)) { 
    print(t.test(iris$Sepal.Width[iris$Species==species])) 
}
