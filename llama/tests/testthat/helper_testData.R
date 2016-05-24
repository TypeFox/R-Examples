fold = data.frame(a=rep.int(0, 10), b=c(rep.int(1, 5), rep.int(0, 5)),
                  c=c(rep.int(0, 5), rep.int(1, 5)))
d = list(data=rbind(cbind(fold, id=1:10), cbind(fold, id=11:20)),
         train=list(1:nrow(fold)),
         test=list(1:nrow(fold) + nrow(fold)),
         features=c("a"),
         ids=c("id"),
         minimize=T,
         performance=c("b", "c"),
         best=rep.int(c(rep.int("c", 5), rep.int("b", 5)), 2))
class(d) = "llama.data"
attr(d, "hasSplits") = TRUE


folde = data.frame(a=c(rep.int(0, 5), rep.int(1, 5)),
                  b=c(rep.int(1, 5), rep.int(0, 5)),
                  c=rep.int(1, 10))
e = list(data=rbind(cbind(folde, id=1:10), cbind(folde, id=11:20)),
         train=list(1:nrow(folde)),
         test=list(1:nrow(folde) + nrow(folde)),
         features=c("c"), minimize=T,
         performance=c("a", "b"),
         ids=c("id"),
         best=rep.int("b", 20))
class(e) = "llama.data"
attr(e, "hasSplits") = TRUE


foldf = data.frame(a=rep.int(0, 10), b=c(rep.int(1, 5), rep.int(0, 5)),
                  c=c(rep.int(0, 5), rep.int(1, 5)))
dnosplit = list(data=rbind(cbind(foldf, id=1:10), cbind(foldf, id=11:20)),
         features=c("a"),
         ids=c("id"),
         minimize=T,
         performance=c("b", "c"),
         best=rep.int(c(rep.int("c", 5), rep.int("b", 5)), 2))
class(dnosplit) = "llama.data"


foldg = data.frame(a=rep.int(0, 10), b=rep.int(1, 10), c=rep.int(0, 10),
                   d=rep.int(T, 10), e=rep.int(F, 10))
g = list(data=rbind(cbind(foldg, id=1:10), cbind(foldg, id=11:20)),
         train=list(1:nrow(foldg)),
         test=list(1:nrow(foldg) + nrow(foldg)),
         features=c("a"),
         performance=c("b", "c"),
         success=c("d", "e"),
         ids=c("id"),
         minimize=T,
         best=rep.int("b", 20))
class(g) = "llama.data"
attr(g, "hasSplits") = TRUE

bests = c("a", "a", "a", "b", "b")
bestlist = list("a", "b", c("a", "b"))
bestlistlong = list("a", "a", "b", c("a", "b"), "a", "a", c("a", "b"), "a", "a", "a")

foldmeas = data.frame(a=rep.int(1, 5), b=rep.int(0, 5),
        d=rep.int(F, 5), e=rep.int(T, 5))
dmeas = list(data=rbind(cbind(foldmeas, id=1:5), cbind(foldmeas, id=6:10)),
    test=list(1:5, 6:10), performance=c("a", "b"), minimize=T, ids=c("id"),
    success=c("d", "e"))

asmeas = data.frame(algorithm=rep.int("a", 5), score=1, iteration=1)
bsmeas = data.frame(algorithm=rep.int("b", 5), score=1, iteration=1)
modelameas = list(predictions=rbind(cbind(asmeas, id=1:5), cbind(asmeas, id=6:10)))
class(modelameas) = "llama.model"
attr(modelameas, "hasPredictions") = TRUE
modelbmeas = list(predictions=rbind(cbind(bsmeas, id=1:5), cbind(bsmeas, id=6:10)))
class(modelbmeas) = "llama.model"
attr(modelbmeas, "hasPredictions") = TRUE

foldone = data.frame(a=c(rep.int(0, 10)), b=c(rep.int(1, 10)), c=rep.int(1, 10))
one = list(data=rbind(cbind(foldone, id=1:10), cbind(foldone, id=11:20)),
         train=list(1:nrow(foldone)),
         test=list(1:nrow(foldone) + nrow(foldone)),
         features=c("c"), minimize=T,
         performance=c("a", "b"),
         ids=c("id"),
         best=rep.int("a", 20))
class(one) = "llama.data"
attr(one, "hasSplits") = TRUE
