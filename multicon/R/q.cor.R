q.cor <-
function(x, set, sex, fem=1, male=2, tails=2, sims=1000, seed=2) {
  new.data <- data.frame(cbind(sex, x, set))
  Comb.N <- nrow(subset(new.data, complete.cases(new.data)))
  female <- subset(new.data, sex == fem)
  Fem.N <- nrow(subset(female, complete.cases(female)))
  fem.cols <- female[,-1:-2]
  male <- subset(new.data, sex == male)
  Male.N <- nrow(subset(male, complete.cases(male)))
  male.cols <- male[,-1:-2]
  Ns <- cbind(Comb.N, Fem.N, Male.N)
  rownames(Ns) <- c("N")
  colnames(Ns) <- c("Combined", "Female", "Male")
  comb.cor <- data.frame(t(cor(x, set, use="pairwise")))
  names(comb.cor) <- c("Combined")
  comb.sig <- data.frame(sig.r(abs(comb.cor$Combined), Comb.N, tails))
  names(comb.sig) <- c("p.C")
  comb.out <- data.frame(cbind(comb.cor, comb.sig))
  fem.cor <- data.frame(t(cor(female[,2], fem.cols, use="pairwise")))
  names(fem.cor) <- c("Female")
  fem.sig <- data.frame(sig.r(abs(fem.cor$Female), Fem.N, tails))
  names(fem.sig) <- c("p.F")
  fem.out <- data.frame(cbind(fem.cor, fem.sig))
  male.cor <- data.frame(t(cor(male[,2], male.cols, use="pairwise")))
  names(male.cor) <- c("Male")
  male.sig <- data.frame(sig.r(abs(male.cor$Male), Male.N, tails))
  names(male.sig) <- c("p.M")
  male.out <- data.frame(cbind(male.cor, male.sig))
  cortable <- data.frame(cbind(comb.out, fem.out, male.out))
  vec.corr <- cor(cortable$Female, cortable$Male)
  cortable.order <- cortable[order(cortable$Combined, cortable$Female, cortable$Male, decreasing=T),]
  comb.order.sig <- data.frame(sig.r(abs(cortable.order$Combined), Comb.N, tails))
  fem.order.sig <- data.frame(sig.r(abs(cortable.order$Female), Fem.N, tails))
  male.order.sig <- data.frame(sig.r(abs(cortable.order$Male), Male.N, tails))
  names(comb.order.sig) <- c("p.C")
  names(fem.order.sig) <- c("p.F")
  names(male.order.sig) <- c("p.M")
  sort.out <- data.frame(cbind(cortable.order$Combined, comb.order.sig, cortable.order$Female, fem.order.sig, cortable.order$Male, male.order.sig))
  names(sort.out) <- c("Combined", "p.C", "Female", "p.F", "Male", "p.M")
  rownames(sort.out) <- rownames(cortable.order)

  if (sims != FALSE) {
    all.test <- rand.test(data.frame(x), set, sims=1000, graph=F, seed=seed)
    fem.test <- rand.test(data.frame(female[,2]), fem.cols, sims=sims, graph=F, seed=seed)
    male.test <- rand.test(data.frame(male[,2]), male.cols, sims=sims, graph=F, seed=seed)
    result <- list("N"=Ns, "corrs"=cortable, "sorted"=sort.out, "vector.cor"=vec.corr, "combined"=all.test, "fem.test"=fem.test, "male.test"=male.test)
    }
  if (sims == FALSE) {
    result <- list("N"=Ns, "corrs"=cortable, "sorted"=sort.out, "vector.cor"=vec.corr)
  }
  class(result) <- c("q.cor")
  return(result)
}
