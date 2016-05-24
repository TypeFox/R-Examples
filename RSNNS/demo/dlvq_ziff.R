library(RSNNS)

data(snnsData)

dataset <- snnsData$dlvq_ziff_100.pat

inputs <- dataset[,inputColumns(dataset)]
outputs <- dataset[,outputColumns(dataset)]

model <- dlvq(inputs, outputs)

mean(fitted(model) - outputs)

