
data(reinis)
reinis <- as.gmData(reinis)

m1 <- hllm(~.^. , reinis) 
m1 <- fit(m1,engine="loglm")

m1res <- stepwise(m1)
formula(m1res)

m2 <- hllm(~smoke*phys*protein+mental*phys+mental*family+smoke*systol*protein, reinis) 
#dynamic.Graph(m2)

