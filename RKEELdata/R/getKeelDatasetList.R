#Returns the list of keel datasets
getKeelDatasetList <- function(){

  dataList <- c(
                #Classification
                "abalone",
                "abalone19",
                "banana",
                "breast",
                "bupa",
                "car",
                "car_train",
                "car_test",
                "coil2000",
                "dermatology",
                "ecoli",
                "ecoli1",
                "ecoli4",
                "glass",
                "hepatitis",
                "iris",
                "iris_train",
                "iris_test",
                "kr-vs-k",
                "mammographic",
                "marketing",
                "monk-2",
                "mushroom",
                "new-thyroid1",
                "pima",
                "shuttle",
                "thyroid",
                "tic-tac-toe",
                "winequality-red",
                "winequality-red-4",
                "winequality-white",
                "yeast",
                "yeast1",
                "yeast6",
                "zoo",
                
                #Regression
                "autoMPG6",
                "autoMPG6_train",
                "autoMPG6_test",
                "delta_elv",
                "diabetes",
                "forestFires",
                "friedman",
                "mv",
                "plastic"
                )

  return(dataList)
}
