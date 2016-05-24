### R code from vignette source 'relations.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: relations.Rnw:53-56
###################################################
options(width = 80)
library("sets")
library("relations")


###################################################
### code chunk number 2: relationgenerator
###################################################
## A relation created by specifying the graph:
R <- relation(graph = data.frame(A = c(1, 1:3), B = c(2:4, 4)))
## extract domain
relation_domain(R)
## extract graph
relation_graph(R)
## both ("a pair of domain and graph" ...)
as.tuple(R)
## extract incidence
relation_incidence(R)

## (Almost) the same using the set specification
## (the domain labels are missing).
R <- relation(graph = set(tuple(1,2), tuple(1,3), tuple(2,4), tuple(3,4)))
## equivalent to:
## relation(graph = list(1:2, c(1,3), c(2,4), c(3,4)))
relation_incidence(R)

## Domains can be composed of arbitrary R objects:
R <- relation(domain = set(c, "test"),
              graph = set(tuple(c, c), tuple(c, "test")))
relation_incidence(R)

as.relation(1:3)
relation_graph(as.relation(c(TRUE, FALSE, TRUE)))
relation_graph(as.relation(factor(c("A", "B", "A"))))


###################################################
### code chunk number 3: charfun
###################################################
divides <- function(a, b) b %% a == 0
R <- relation(domain = list(1 : 10, 1 : 10), charfun = divides)
R
"%|%" <- relation_charfun(R)

2L %|% 6L
2:4 %|% 6L
2L %|% c(2:3, 6L)

"%|%"(2L, 6L)


###################################################
### code chunk number 4: predicates
###################################################
R <- as.relation(1:5)
relation_is_binary(R)
relation_is_transitive(R)
relation_is_partial_order(R)


###################################################
### code chunk number 5: ops1
###################################################
x <- matrix(0, 3L, 3L)
R1 <- as.relation(row(x) >= col(x))
R2 <- as.relation(row(x) <= col(x))
R3 <- as.relation(row(x) <  col(x))
relation_incidence(max(R1, R2))
relation_incidence(min(R1, R2))
R3 < R2
relation_dissimilarity(min(R1, R2), max(R1, R2))


###################################################
### code chunk number 6: ops2
###################################################
relation_incidence(R1 * R2)
relation_incidence(! R1)
relation_incidence(t(R2))


###################################################
### code chunk number 7: plot
###################################################
ps <- 2 ^ set("a", "b", "c")
inc <- set_outer(ps, "<=")
if (require("Rgraphviz")) plot(relation(incidence = inc))


###################################################
### code chunk number 8: plotfig
###################################################
ps <- 2 ^ set("a", "b", "c")
inc <- set_outer(ps, "<=")
if (require("Rgraphviz")) plot(relation(incidence = inc))


###################################################
### code chunk number 9: relations.Rnw:365-369
###################################################
data("Cetacea")
ind <- sapply(Cetacea, function(s) all(!is.na(s)))
relations <- as.relation_ensemble(Cetacea[, ind])
print(relations)


###################################################
### code chunk number 10: relations.Rnw:374-377
###################################################
any(duplicated(relations))
thrice <- c(rep(relations, 2L), relations)
all.equal(unique(thrice), relations)


###################################################
### code chunk number 11: relations.Rnw:382-383
###################################################
all.equal(thrice[!duplicated(thrice)], relations)


###################################################
### code chunk number 12: relations.Rnw:388-389
###################################################
relation_dissimilarity(relations[1 : 2], relations["CLASS"])


###################################################
### code chunk number 13: relations.Rnw:393-395
###################################################
d <- relation_dissimilarity(relations)
sort(as.matrix(d)[, "CLASS"])[-1L]


###################################################
### code chunk number 14: relations.Rnw:400-402
###################################################
complement <- !relations
complement


###################################################
### code chunk number 15: projection
###################################################
## projection
Person <-
    data.frame(Name = c("Harry", "Sally", "George", "Helena", "Peter"),
               Age = c(34, 28, 29, 54, 34),
               Weight = c(80, 64, 70, 54, 80),
               stringsAsFactors = FALSE)
Person <- as.relation(Person)
relation_table(Person)
relation_table(relation_projection(Person, c("Age", "Weight")))


###################################################
### code chunk number 16: selection
###################################################
## selection
relation_table(R1 <- relation_selection(Person, Age < 29))
relation_table(R2 <- relation_selection(Person, Age >= 34))
relation_table(R3 <- relation_selection(Person, Age == Weight))


###################################################
### code chunk number 17: unioncomplement
###################################################
## union
relation_table(R1 %U% R2)

## works only for the same domains:
relation_table(R2 | R3)

## complement
relation_table(Person - R2)


###################################################
### code chunk number 18: intersectionsymdiff
###################################################
## intersection
relation_table(relation_intersection(R2, R3))

## works only for the same domains:
relation_table(R2 & R3)

## symmetric difference
relation_table(relation_symdiff(R2, R3))


###################################################
### code chunk number 19: cartesian
###################################################
## cartesian product
Employee <-
    data.frame(Name = c("Harry", "Sally", "George", "Harriet", "John"),
               EmpId = c(3415, 2241, 3401, 2202, 3999),
               DeptName = c("Finance", "Sales", "Finance", "Sales", "N.N."),
	       stringsAsFactors = FALSE)
Employee <- as.relation(Employee)
relation_table(Employee)
Dept <- data.frame(DeptName = c("Finance", "Sales", "Production"),
                   Manager = c("George", "Harriet", "Charles"),
                   stringsAsFactors = FALSE)
Dept <- as.relation(Dept)
relation_table(Dept)

relation_table(Employee %><% Dept)


###################################################
### code chunk number 20: division
###################################################
## division
Completed <-
    data.frame(Student = c("Fred", "Fred", "Fred", "Eugene",
                           "Eugene", "Sara", "Sara"),
               Task = c("Database1", "Database2", "Compiler1",
                        "Database1", "Compiler1", "Database1",
                        "Database2"),
               stringsAsFactors = FALSE)
Completed <- as.relation(Completed)
relation_table(Completed)
DBProject <- data.frame(Task = c("Database1", "Database2"),
                        stringsAsFactors = FALSE)
DBProject <- as.relation(DBProject)
relation_table(DBProject)

relation_table(Completed %/% DBProject)

## division remainder
relation_table(Completed %% DBProject)


###################################################
### code chunk number 21: naturaljoin
###################################################
## Natural join
relation_table(Employee %|><|% Dept)

## left (outer) join
relation_table(Employee %=><% Dept)

## right (outer) join
relation_table(Employee %><=% Dept)

## full outer join
relation_table(Employee %=><=% Dept)


###################################################
### code chunk number 22: semijoin
###################################################
## semijoin
relation_table(Employee %|><% Dept)
relation_table(Employee %><|% Dept)


###################################################
### code chunk number 23: antijoin
###################################################
## antijoin
relation_table(Employee %|>% Dept)
relation_table(Employee %<|% Dept)


###################################################
### code chunk number 24: consensus1a
###################################################
data("Felines")
relations <- as.relation_ensemble(Felines)


###################################################
### code chunk number 25: consensus1b
###################################################
E <- relation_consensus(relations, "SD/E")

ids <- relation_class_ids(E)
split(rownames(Felines), ids)


###################################################
### code chunk number 26: consensus2a
###################################################
pm <- matrix(c(0, 1, 0, 1, 1,
               0, 0, 0, 1, 1,
               1, 1, 0, 0, 0,
               0, 0, 1, 0, 0,
               0, 0, 1, 1, 0),
             nrow = 5L,
             byrow = TRUE,
             dimnames = list(letters[1:5], letters[1:5]))
R <- as.relation(t(pm))
relation_incidence(R)
relation_is_tournament(R)


###################################################
### code chunk number 27: consensus2b
###################################################
L <- relation_consensus(R, "SD/L")
relation_incidence(L)


###################################################
### code chunk number 28: relations.Rnw:747-748
###################################################
relation_class_ids(L)


###################################################
### code chunk number 29: consensus2c
###################################################
L <- relation_consensus(R, "SD/L", control = list(all = TRUE))
print(L)
if(require("Rgraphviz")) plot(L)


###################################################
### code chunk number 30: relations.Rnw:759-760
###################################################
lapply(L, relation_class_ids)


###################################################
### code chunk number 31: consensus2d
###################################################
W3 <- relation_consensus(R, "SD/W", control = list(k = 3))
relation_incidence(W3)
relation_class_ids(W3)


###################################################
### code chunk number 32: consensusfig
###################################################
if(require("Rgraphviz")) plot(L)


