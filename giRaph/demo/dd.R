
library(giRaph)

G<-new("incidenceList",
       V=letters[1:12],
       E=list(
              d(6,5,c(2,4),c(1,3)),
              u(2,4,5),
              d(2,4),d(4,2),
              d(1,7),d(3,7),d(4,7),
              d(5,8),d(5,8),d(5,8),
              u(6,9),d(6,9),
              u(9,9),
              d(9,8),d(9,12),
              u(7,8),u(8,12),u(12,11),u(11,7),
              u(11,8),
              d(11,10)
             )
      )
I<-as(G,"incidenceMatrix")
A<-as(I,"adjacencyList")
X<-as(I,"adjacencyMatrix")

ag<-new("anyGraph",incidenceList=G)
gg<-new("generalGraph",incidenceList=G)
mg<-new("multiGraph",incidenceList=G)
sg<-new("simpleGraph",incidenceList=G)

# Subsection "A class for vertex sets"

V<-new("vertexSet",seq(1,6),seq(4,9))
names(V)
str(V)
is(V)

show(W<-v(seq(4,9),seq(1,6)))
card(W)
isEmpty(W)
areTheSame(V,W)

V[seq(1,9,2)]
W[[1]]

v("a","b")+v("b","c")
v("a","b")*v("b","c")
v("a","b")-v("b","c")

# Subsection "Classes for edges"

show(U<-new("vertexSet",letters[1:12]))

u(2,4,5)
u(11,12)
u(9)
d(6,5,c(2,4),c(1,3))
r(10,11)

showRel(u(11,12),code=U)
showRel(d(6,5,c(2,4),c(1,3)),code=U)

maxId(u(1,13))<=card(U)

recode(u(11,12),src=U,dst=U[c(1,2,11,12)])

card(d(6,5,c(2,4),c(1,3)))
length(d(6,5,c(2,4),c(1,3)))
card(u(11,12))==length(u(11,12))

areTheSame(u(11,12),u(12,11))
areTheSame(d(6,5,c(2,4),c(1,3)),d(6,5,c(4,2),c(3,1)))

u(2,4,5)[1:2]
u(2,4,5)[[3]]
d(6,5,c(2,4),c(1,3))[c(1,2,3)]
d(6,5,c(2,4),c(1,3))[[4]]

# Subsection "A class for multi-sets of edges"

show(E<-new("edgeList",u(11,12),d(11,10)))

maxId(E)
recode(E,src=U,dst=U[10:12])

card(new("edgeList",u(11,12),d(11,10),d(11,10)))
isPresent(d(11,12),E)

areTheSame(E,new("edgeList",u(11,12),d(11,10),d(11,10)))
areTheSame(E,new("edgeList",d(11,10),u(11,12)))

E[1]
E[[1]]

new("edgeList")+u(11,12)
new("edgeList",u(11,12),d(11,10),u(11,12))-u(11,12)

# Subsection "Classes for representations"

new("incidenceList",V=c("f","i"),E=list(u(1,2),d(1,2),u(2,2)))

new("incidenceMatrix",matrix(c(c(-1,0,3,2.1,2.1),rep(0,5)),2,5,byrow=T))

new("adjacencyList",id=letters[1:3],ch=list(c(2,3),0,0),ne=list(0,3,0))

new("adjacencyMatrix",matrix(c(c(0,1,1),c(0,-1,0),c(2,0,0)),3,3,byrow=T))


names(G)<-letters[13:24]
G[[1]]
G[1:6]
str(card(G[1:6]))
isPresent(d(2,4),G[1:6])
names(G)<-letters[1:12]
I<-as(G,"incidenceMatrix")
I[1:4]
A<-as(I,"adjacencyList")
areTheSame(as(G[c(5,7,8)],"adjacencyList"),A[c(5,7,8)])
X<-as(A,"adjacencyMatrix")
areTheSame(as(X,"adjacencyList"),A)

names(G[10:12]+v("y","z"))
areTheSame(G[1:6],G-v(letters[7:12]))
areTheSame(G[1:6],G*v(letters[1:6]))
isPresent(u(2,4,5),G-u(2,4,5))
(I[1:4]+u(3,4))-d(2,4)
areTheSame(A[1:6]*v("b","d"),A[c(2,4)])
isPresent(r(1,2),X[10:12])
X[10:12]+d(1,2)
isPresent(u(1,2),X[10:12]+d(1,2))
isPresent(r(1,2),X[10:12]+d(1,2))

# Classes for graphs

isPresent(d(2,4),ag)
isPresent(u(2,4),ag)
isPresent(u(2,4,5),ag)
isPresent(d(2,4),gg)
isPresent(u(2,4),gg)
isPresent(u(2,4,5),gg)
isPresent(d(2,4),mg)
isPresent(u(2,4),mg)
isPresent(u(2,4,5),mg)
isPresent(d(2,4),sg)
isPresent(u(2,4),sg)
isPresent(u(2,4,5),sg)

areTheSame(ag,gg)
areTheSame(gg,mg)
