require(FrF2)
set.seed(3456)
planblock1 <- FrF2(32,7,blocks=list(c(1,2,3),c(2,4,5),c(3,6,7)), alias.block.2fi=TRUE)
set.seed(3456)
planblock2 <- FrF2(32,7,blocks=8, alias.block.2fi=TRUE)
set.seed(3456)
planblockrepl <- FrF2(32,7,blocks=8, alias.block.2fi=TRUE, bbreps=2)
set.seed(3456)
planblockrepe <- FrF2(32,7,blocks=8, alias.block.2fi=TRUE, wbreps=2, repeat.only=TRUE)
set.seed(3456)
planblockreplwithin <- FrF2(32,7,blocks=8, alias.block.2fi=TRUE, wbreps=2)
set.seed(3456)
planblockreplboth <- FrF2(32,7,blocks=8, alias.block.2fi=TRUE, wbreps=2, bbreps=2)
set.seed(3456)
planblockreplrepe <- FrF2(32,7,blocks=8, alias.block.2fi=TRUE, wbreps=2, bbreps=2, repeat.only=TRUE)
set.seed(3456)
planest <- FrF2(16,6,estimable=c("AB","AC","AD","BC","CD"),clear=FALSE)

set.seed(3456)
plansimple <- FrF2(16,7)
set.seed(3456)
planrepl <-FrF2(16,7, repl=2)
set.seed(3456)
planrepe <-FrF2(16,7, repl=2, repe=TRUE)

test1 <- add.center(plansimple,3)  ## 3 center points distributed through the plan (3 positions)
test1b <- add.center(plansimple,3, distribute=1) ## 3 center points at the end
test1c <- add.center(plansimple,6, distribute=4)  ## 6 center points distributed through the plan (4 positions)

test2 <- add.center(planest, 3)
test2b <-add.center(planest, 5, distribute=1)
test2c <-add.center(planest, 10, distribute=7)

test4 <- add.center(planrepl,3)  ## 3 center points for proper replications
test4b <- add.center(planrepl,3, distribute=1)  ## 3 center points for proper replications

test5 <- add.center(planrepe,3)  ## 3 center points for repeated measurements
test5b <- add.center(planrepe,3, distribute=1)  ## 3 center points for repeated measurements

test6 <- add.center(planblock1, 3)  ## center points in blocked plan
test6b <- add.center(planblock1, 3, distribute=1)  ## center points in blocked plan

test7 <- add.center(planblock2, 3)  ## center points in blocked plan
test7b <- add.center(planblock2, 3, distribute=1)  ## center points in blocked plan

test8 <- add.center(planblockrepl, 3)  ## center points in replicated blocked plan
test8b <- add.center(planblockrepl, 3, distribute=1)  ## center points in replicated blocked plan

test9 <- add.center(planblockrepe, 3)  ## center points in replicated blocked plan
test9b <- add.center(planblockrepe, 3, distribute=1)  ## center points in replicated blocked plan

test10 <- add.center(planblockreplwithin, 3)  ## center points in replicated blocked plan
test10b <- add.center(planblockreplwithin, 3, distribute=1)  ## center points in replicated blocked plan

test11 <- add.center(planblockreplboth, 3)  ## center points in replicated blocked plan
test11b <- add.center(planblockreplboth, 3, distribute=1)  ## center points in replicated blocked plan

test12 <- add.center(planblockreplrepe, 3)  ## center points in replicated blocked plan
test12b <- add.center(planblockreplrepe, 3, distribute=1)  ## center points in replicated blocked plan

rownames(test1)
rownames(test1b)
rownames(test1c)
rownames(test2)
rownames(test2b)
rownames(test2c)

cbind(run.order(test1),test1)
cbind(run.order(test2),test2)
cbind(run.order(test4b),test4b)
cbind(run.order(test5),test5)
cbind(run.order(test6),test6)
cbind(run.order(test7),test7)
cbind(run.order(test8),test8)
cbind(run.order(test9),test9)
cbind(run.order(test10),test10)
cbind(run.order(test11),test11)
cbind(run.order(test12),test12)
getblock(planblockreplrepe)

rerandomize.design(test1)
rerandomize.design(test2)
rerandomize.design(test4)
rerandomize.design(test5b)
rerandomize.design(test6)
rerandomize.design(test7)
rerandomize.design(test8)
rerandomize.design(test9)
rerandomize.design(test10)
rerandomize.design(test11)
rerandomize.design(test12)
