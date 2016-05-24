Random.Smiles <- do(2000) * diff(mean(Leniency ~ shuffle(Group), data = Smiles))
head(Random.Smiles, 3)
histogram(~smile, n = 24,, fit = "normal", data = Random.Smiles)

