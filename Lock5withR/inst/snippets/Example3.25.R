head(ExerciseHours)
favstats(~Exercise|Gender, data = ExerciseHours)
stat <- diffmean(Exercise~Gender, data = ExerciseHours); stat

