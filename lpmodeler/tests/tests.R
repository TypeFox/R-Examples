require(lpmodeler)

p <- newProblem()
print(p)
p

p <- addVariable(p, "C", 5, "x")
p <- addVariable(p, "C", 4, "y")
p
ls(p$ctrs)
ls(p$vars)

p <- addConstraint(p, ">=", 5, c(1, 2), name = "x+2y greater or equal than 5")
p
ls(p$ctrs)
ls(p$vars)

p <- addConstraint(p, "<=", 10, name = "less or equal than 10")
p
ls(p$ctrs)
ls(p$vars)

# Solve a simple LP
# maximize:   2 x1 + 4 x2 + 3 x3
p <- newProblem()

p <- addVariable(p, "C", 2, name="x1")
p <- addVariable(p, "C", 4, name="x2")
p <- addVariable(p, "C", 3, name="x3")

# subject to: 3 x1 + 4 x2 + 2 x3 <= 60
p <- addConstraint(p, "<=", 60, name="c1")
p <- setPoint(p, "x1", "c1", 3)
p <- setPoint(p, "x2", "c1", 4)
p <- setPoint(p, "x3", "c1", 3)

#             2 x1 +   x2 +   x3 <= 40
p <- addConstraint(p, "<=", 40, c(2, 1, 1))

#               x1 + 3 x2 + 2 x3 <= 80
p <- addConstraint(p, "<=", 80, c(1, 3, 2))

#               x1 >= 0
p <- addConstraint(p, ">=", 0, c(1, 0, 0))

#               x2 >= 0
p <- addConstraint(p, ">=", 0, c(0, 1, 0))

#               x3 >= 0
p <- addConstraint(p, ">=", 0, c(0, 0, 1))

if(require(Rsymphony))
  mipSolve(p)

# Simple mixed integer linear program.
# maximize:    3 x1 + 1 x2 + 3 x3
# subject to: -1 x1 + 2 x2 +   x3 <= 4
#                     4 x2 - 3 x3 <= 2
#                x1 - 3 x2 + 2 x3 <= 3
#                x1 >= 0 (integer)
#                x2 >= 0 (real)
#                x3 >= 0 (integer)
p <- newProblem()
p <- addVariable(p, "I", 3)
p <- addVariable(p, "C", 1)
p <- addVariable(p, "I", 3)
p <- addConstraint(p, "<=", 4, c(-1, 2, 1))
p <- addConstraint(p, "<=", 2, c(0, 4, -3))
p <- addConstraint(p, "<=", 3, c(1, -3, 2))
p <- addConstraint(p, ">=", 0, c(1, 0, 0))
p <- addConstraint(p, ">=", 0, c(0, 1, 0))
p <- addConstraint(p, ">=", 0, c(0, 0, 1))

if(require(Rsymphony))
  mipSolve(p)
