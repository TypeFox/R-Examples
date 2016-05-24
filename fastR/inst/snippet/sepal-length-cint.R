for (species in levels(iris$Species)) { 
    print(t.test(iris$Sepal.Length[iris$Species==species])$conf.int) 
}
