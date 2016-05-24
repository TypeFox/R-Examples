BootE <- do(2000) * diff(mean(Exercise~Gender, data = resample(ExerciseHours)))
head(BootE, 3)
histogram(~M, width = 0.5, fit = "normal", data = BootE)

