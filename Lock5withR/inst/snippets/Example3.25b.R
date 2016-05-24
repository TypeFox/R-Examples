BootE <- do(3000) * diffmean(Exercise~Gender, data = resample(ExerciseHours))
head(BootE, 3)

