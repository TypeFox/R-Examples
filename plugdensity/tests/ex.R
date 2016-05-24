library("plugdensity")

options(digits = 6)

data(faithful)
(pd.geys <- plugin.density(faithful$waiting))
pd.geys$y * 1e4
