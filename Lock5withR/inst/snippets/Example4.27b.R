Randomization.Smiles <- do(1000) * diffmean(Leniency ~ shuffle(Group), data = Smiles)
head(Randomization.Smiles, 3)

