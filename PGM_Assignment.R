 #1.3 R–packages
 #source("http://bioconductor.org/biocLite.R"); 
 #biocLite(c("graph", "RBGL", "Rgraphviz"))
 #install.packages("gRbase", dependencies=TRUE); 
 #install.packages("gRain", dependencies=TRUE); 
 #install.packages("gRim", dependencies=TRUE)



# 1.4 The practicals: The coronary artery disease data
data(cad1)
head(cad1)

data(cad2)
head(cad2)


# 3. Undirected Graphs
library(gRbase)
g1 <- ug(~a:b:e + a:c:e + b:d:e + c:d:e + c:g + d:f)
class(g1)

as(g1, "matrix")

library(Rgraphviz)
plot(g1)

m <- matrix(c(0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0,
              1, 1, 0, 1, 1, 1, 1, 1, 0), nrow = 5)
rownames(m) <- colnames(m) <- c("a", "b", "c", "d", "e")
m
as(m, "graphNEL")

g1a <- addEdge("a", "d", g1)
g1b <- removeEdge("c", "d", g1)
par(mfrow = c(1, 3))
plot(g1, main = "g1")
plot(g1a, main = "g1a")
plot(g1b, main = "g1b")

g1c <- subGraph(c("b", "c", "d", "e"), g1)
par(mfrow = c(1, 3))
plot(g1, main = "g1")
plot(g1c, main = "g1c")

is.complete(g1,set=c("a","b","e"))
is.complete(g1)

str(maxClique(g1))

g2 <- ug(~a:b:e + a:c:e + b:d:e + c:d:e)
plot(g2)
separates("a", "d", c("b", "c", "e"), g2)

# 3.1 Factorization and dependence graph
plot((g3 <- ug(~ A:B + B:C:D + C:E + D:E)))

# 3.2 Reading conditional independencies – global Markov property
plot(g3)
separates(c("D","E"), "A", "B", g3)

# 4 Directed acyclic graphs (DAGs)
dag0 <- dag(~a, ~b * a, ~c * a * b, ~d * c * e, ~e * a)
dag0 <- dag(~a + b:a + c:a:b + d:c:e + e:a)
dag0 <- dag("a", c("b", "a"), c("c", "a", "b"), c("d", "c", "e"), c("e", "a"))
dag0

plot(dag0)

(m <- matrix(c(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
               1, 0, 0, 0, 0, 0, 1, 1, 0), nrow = 5))
rownames(m) <- colnames(m) <- letters[1:5]
dg <- as(m, "graphNEL")
plot(dg)

parents("d", dag0)
children("c", dag0)

# 4.1 Factorization and dependence graph – DAGs
plot(dag0)

# 4.2 Reading conditional independencies from DAGs (I)
par(mfrow=c(3,1),mar=c(3,1,2,1))
plot(dag(~a+b:a+c:b),"circo")
plot(dag(~c+b:c+a:b),"circo")
plot(dag(~b+a:b+c:b),"circo")

plot(dag(~a+c+b:a:c),"circo")

# 4.3 Moralization
dag0m <- moralize(dag0)
par(mfrow=c(1,2))
plot(dag0)
plot(dag0m)

# 4.4 Ancestral sets and graphs
ancestralSet(c("a", "c", "e"), dag0)
plot(ancestralGraph(c("a", "c", "e"), dag0))

# 4.5 Reading conditional independences from DAG (II)
par(mfrow=c(1,2))
plot(ancestralGraph(c("a", "c", "e"), dag0))
plot(moralize(ancestralGraph(c("a", "c", "e"), dag0)))

# 5 Bayesian Network (BN) basics
# 6 A small worked example BN
plot((FTH<-dag(~ F + T:F + H:T)), "circo")

# 6.1 Specification of conditional probability tables
(p.F <- parray("F", levels=2, values=c(.01,.99)))
(p.TgF <- parray(c("T","F"), levels=c(2,2), values=c(.95,.05, .001,.999)))
(p.HgT <- parray(c("H","T"), levels=c(2,2), values=c(.80,.20, .010,.990)))

# 6.2 Brute force computations
# 1) Calculate joint distribution p(F T H)
p.FT <- tableMult(p.F, p.TgF)
p.FTH <- tableMult(p.FT, p.HgT)
as.data.frame.table(p.FTH)
# 2) Calculate the marginal distribution p(F H)
p.FH <- tableMargin(p.FTH, margin=c('F','H'))
as.data.frame.table(p.FH)
# 3) calculate conditional distribution p(I|H)
p.H<- tableMargin(p.FH, margin='H')
(p.FgH <- tableDiv(p.FH, p.H))

# 7 Decomposable graphs and junction trees
# 7.1 Decomposable graphs
par(mfrow=c(1,3))
plot(ug(~1:2+2:3:4+3:4:5:6+6:7), "circo") # decomposable
plot(ug(~1:2+2:3+3:4:5+4:1),"circo") # not decomposable
plot(ug(~1:2:5+2:3:5+3:4:5+4:1:5),"circo") # not decomposable

# 7.4 Computations by message passing
par(mfrow=c(1,2), oma=c(2,1,2,1))
plot(FTH)
plot(moralize(FTH))

# 7.5 Clique potential representation
(qFT <- tableMult(p.F, p.TgF))
(qTH <- p.HgT)
(qT <- parray("T",levels=2, values=1))

# 7.6 Working inwards in junction tree
(qTs <- tableMargin(qFT, "T"))
(qTHs <- tableMult(qTH, tableDiv(qTs, qT)))

# 7.7 Working outwards in junction tree
(qTss <- tableMargin(qTHs, "T"))
(qFTs <- tableMult(qFT, tableDiv(qTss, qTs)))
qFTs
qTHs
qTs # probability of temperature
tableMargin(qFT, "F") # probability of fever
tableMargin(qTH, "H") # probability of headache

# 8. Propagating findings
qTH
## Set finding H=H1
qTH[c(2,4)] <- 0
qTH
## Repeat everything
(qTs <- tableMargin(qFT, "T"))
(qTHs <- tableMult(qTH, tableDiv(qTs, qT)))
(qTss <- tableMargin(qTHs, "T"))
(qFTs <- tableMult(qFT, tableDiv(qTss, qTs)))
sum(qFTs)
tableMargin(qFTs, "F")/sum(qFTs)
tableMargin(qFTs, "T")
tableMargin(qTHs, "T")
qTss


# 10 An introduction to the gRain package
yn <- c("yes","no")
a <- cptable(~asia, values=c(1,99),levels=yn)
t.a <- cptable(~tub+asia, values=c(5,95,1,99),levels=yn)
s <- cptable(~smoke, values=c(5,5), levels=yn)
l.s <- cptable(~lung+smoke, values=c(1,9,1,99), levels=yn)
b.s <- cptable(~bronc+smoke, values=c(6,4,3,7), levels=yn)
e.lt <- cptable(~either+lung+tub,values=c(1,0,1,0,1,0,0,1),levels=yn)
x.e <- cptable(~xray+either, values=c(98,2,5,95), levels=yn)
d.be <- cptable(~dysp+bronc+either, values=c(9,1,7,3,8,2,1,9), levels=yn)
plist <- compileCPT(list(a, t.a, s, l.s, b.s, e.lt, x.e, d.be))
bnet <- grain(plist)
bnet
plist
plist$tub
plot(bnet)

# 10.1 Queries
querygrain(bnet, nodes=c('lung', 'tub', 'bronc'))

# 10.2 Setting findings and probability of findings
bnet.f <- setFinding(bnet, nodes=c('asia', 'dysp'), state=c('yes','yes'))
bnet.f
pFinding(bnet.f)

# 10.3 Queries – II
querygrain(bnet.f, nodes=c('lung', 'tub', 'bronc'))
querygrain(bnet.f, nodes=c('lung', 'tub', 'bronc'), type='joint')

# 10.4 Dependence graph, moralization and triangulation
par(mfrow=c(1,2))
plot(bnet$dag)
plot(moralize(bnet$dag))

# 10.5 Triangulation
par(mfrow=c(1,2))
plot(moralize(bnet$dag))
plot(triangulate(moralize(bnet$dag)))

# 12 Contingency tables
data(lizardRAW, package="gRbase")
head(lizardRAW)
dim(lizardRAW)
data(lizard, package="gRbase")
lizard

# 12.1 Notation
lizard
## Marginal table
tableMargin(lizard, c("species","height"))

# 12.2 Log–linear models
(ll1 <- loglin(lizard, list(c("species","diam"),c("species","height"))))
(ll2 <- MASS::loglm(~species:diam+species:height, data=lizard))

# 12.3 Graphical models and decomposable models
par(mfrow=c(1,3))
plot(ug(~A:B:C + B:C:D))
plot(ug(~A:B + A:C + B:C:D))
plot(ug(~A:B + A:C + B:D + C:D))

# 12.4 ML estimation in decomposable models
n.ds <- tableMargin(lizard, c("diam", "species"))
n.hs <- tableMargin(lizard, c("height", "species"))
n.s <- tableMargin(lizard, c("species"))
## Expected cell counts
(fv <- tableDiv( tableMult(n.ds, n.hs), n.s))
as.data.frame.table(tablePerm(fv, c("diam","height","species")))
as.data.frame.table(fitted(ll2))

# 13 Testing for conditional independence
args(ciTest_table)
ciTest(lizard, set=c("diam","height","species"))
# alternatively
ciTest(lizard, set=~diam+height+species)
ciTest(lizard, ~di+he+s)
ciTest(lizard, c("di","he","sp"))
ciTest(lizard, c(2,3,1))

# 13.1 CI-test – stratification
cit <- ciTest(lizard, set=~diam+height+species, slice.info=T)
cit
names(cit)
cit$slice

# 14 Log–linear models using the gRim package
data(wine, package="gRbase")
head(wine,4)
dim(wine)

wine <- cbind(Cult=wine[,1], as.data.frame(lapply(wine[-1], cut, 2, labels=c('L','H'))))
head(wine)
dim(xtabs(~.,wine))
wine <- wine[,1:4]
head(wine)
dim(xtabs(~.,wine))

mm <- dmod(~Cult:Alch+Alch:Mlca:Ash, data=wine)
mm <- dmod(list(c("Cult","Alch"), c("Alch","Mlca","Ash")), data=wine)
mm <- dmod(~C:Alc+Alc:M:As, data=wine)
mm
str(terms(mm))
formula(mm)

# 14.1 Plotting the dependence graph
plot(mm)
# Default: a graphNEL object
DG <- ugList(terms(mm))
DG
# Alternative: an adjacency matrix
ugList(terms(mm), result="matrix")

# 14.2 Model specification shortcuts
str(terms(dmod(~.^., data=wine))) ## Saturated model
str(terms(dmod(~.^1, data=wine))) ## Independence model
str(terms(dmod(~.^3, data=wine))) ## All 3-factor model
marg <- c("Cult", "Alch", "Mlca")
str(terms(dmod(~.^., data=wine, margin=marg))) ## Saturated model

# 14.3 Altering graphical models
mm <- dmod(~Cult:Alch+Alch:Mlca:Ash, data=wine)
mm2 <- update(mm, list(dedge=~Alch:Ash, aedge=~Cult:Mlca)) # No abbreviations
par(mfrow=c(1,2)); plot(mm); plot(mm2)

# 14.4 Model comparison
mm <- dmod(~Cult:Alch+Alch:Mlca:Ash, data=wine)
mm2 <- update(mm, list(dedge=~Alch:Ash+Alch:Cult)) # No abbreviations
compareModels(mm,mm2)

# 14.5 Decomposable models – deleting edges
par(mfrow=c(1,3))
plot(ug(~A:B:C+B:C:D))
plot(ug(~A:C+B:C+B:C:D))
plot(ug(~A:B+A:C+B:D+C:D))

# 14.6 Decomposable models – adding edges
plot(ug(~A:B+B:C+C:D), "circo")
UG <- ug(~A:B+B:C+C:D)
mcs(UG)
UG1 <- addEdge("A","D",UG)
mcs(UG1)
UG2 <- addEdge("A","C",UG)
mcs(UG2)

# 14.7 Test for adding and deleting edges
mm <- dmod(~C:Alc+Alc:M:As, data=wine)
plot(mm)
testdelete(mm, edge=c("Mlca","Ash"))
mm <- dmod(~C:Alc+Alc:M:As, data=wine)
plot(mm)
testadd(mm, edge=c("Mlca","Cult"))

# 14.8 Model search in log–linear models using gRim
args(stepwise.iModel)
dm1 <- dmod(~.^., data=wine)
dm2 <- stepwise(dm1, details=1)
formula(dm2)
terms(dm2)
par(mfrow=c(1,2))
plot(dm1, "circo")
plot(dm2, "circo")

# 15 From graph and data to network
uG2 <- ugList(terms(dm2))
uG2 <- ugList(list(c("Cult", "Ash"), c("Cult", "Alch"), c("Cult", "Mlca")))
uG1 <- ugList(terms(dm1))
(wine1 <- compile(grain(uG1, data=wine)))
(wine2 <- compile(grain(uG2, data=wine)))
querygrain(wine1, "Cult")
querygrain(wine2, "Cult")
querygrain(setFinding(wine1, c("Ash","Alch"), c("L","H")), "Cult")
querygrain(setFinding(wine2, c("Ash","Alch"), c("L","H")), "Cult")
# naive Bayesian model for CAD data
par(mfrow=c(1,2))
plot(UG <- ug(~Heartfail:CAD+Smoker:CAD+Hyperchol:CAD+AngPec:CAD))
plot(DAG <- dag(~Heartfail:CAD+Smoker:CAD+Hyperchol:CAD+AngPec:CAD))
cadmod1 <- compile(grain(UG, cad1))
cadmod2 <- compile(grain(DAG, cad1))
querygrain(cadmod1, nodes="CAD")
querygrain(cadmod2, nodes="CAD")

# 15.1 Prediction
data(cad2)
head(cad2,3)
args(predict.grain)
pred1 <- predict(cadmod1, resp="CAD", newdata=cad2, type="class")
str(pred1)
pred2 <- predict(cadmod1, resp="CAD", newdata=cad2, type="dist")
str(pred2)

# 15.2 Classification error
table(cad2$CAD)
table(cad2$CAD)/sum(table(cad2$CAD))
tt <- table(cad2$CAD, pred1$pred$CAD)
tt
sweep(tt, 1, apply(tt,1,sum), FUN="/")

library(gRim)
data(cad1, package = "gRbase")
# specify saturated model
m.sat <- dmod(~.^., data = cad1)
# Stepwise backward model selection among decomposable models using AIC as
# selection criterion
m.back1 <- stepwise(m.sat, details = 1)
m.back1
summary(m.back1)
plot(m.back1)
# Specify independence model
m.ind <- dmod(~.^1, data = cad1)
# Stepwise forward selection among decomposable models using AIC as
# selection criterion.
m.forw1 <- stepwise(m.ind, direction = "forward", details = 1)
m.forw1
summary(m.forw1)
plot(m.forw1)
#risk factors
x <- terms(m.back1)[!is.na(sapply(terms(m.back1), function(g) match("CAD", g)))]
setdiff(unique(unlist(x)), "CAD")
x <- terms(m.forw1)[!is.na(sapply(terms(m.forw1), function(g) match("CAD", g)))]
setdiff(unique(unlist(x)), "CAD")
#To bayesian network
bn.back1 <- compile(grain(ugList(terms(m.back1)), data = cad1, smooth = 0.001))
bn.forw1 <- compile(grain(ugList(terms(m.forw1)), data = cad1, smooth = 0.001))
data(cad2, package = "gRbase")
table(cad2$CAD)
classMat <- function(trueLabel, classLabel) {
  tt <- table(trueLabel, classLabel)
  sweep(tt, 1, apply(tt, 1, sum), FUN = "/")
}
class.back1 <- predict(bn.back1, resp = "CAD", newdata = cad2, type = "class")
class.forw1 <- predict(bn.forw1, resp = "CAD", newdata = cad2, type = "class")
prob.back1 <- predict(bn.back1, resp = "CAD", newdata = cad2, type = "dist")
prob.forw1 <- predict(bn.forw1, resp = "CAD", newdata = cad2, type = "dist")
classMat(cad2$CAD, class.back1$pred$CAD)
classMat(cad2$CAD, class.forw1$pred$CAD)
m.back2 <- stepwise(m.sat, k = log(nrow(cad1)), details = 1)
m.forw2 <- stepwise(m.ind, k = log(nrow(cad1)), direction = "forward", details = 1)
bn.back2 <- compile(grain(ugList(terms(m.back2)), data = cad1, smooth = 0.001))
bn.forw2 <- compile(grain(ugList(terms(m.forw2)), data = cad1, smooth = 0.001))

class.back2 <- predict(bn.back2, resp = "CAD", newdata = cad2, type = "class")
class.forw2 <- predict(bn.forw2, resp = "CAD", newdata = cad2, type = "class")

classMat(cad2$CAD, class.back2$pred$CAD)
g.naive <- lapply(setdiff(names(cad1), "CAD"), function(x) c(x, "CAD"))
str(g.naive)
bn.naive <- compile(grain(ugList(g.naive), data = cad1))
class.naive <- predict(bn.naive, resp = "CAD", newdata = cad2, type = "class")
classMat(cad2$CAD, class.naive$pred$CAD)
library(rpart)
tree1 <- rpart(CAD ~ ., data = cad1)
plot(tree1, uniform = T, margin = 0.2)
text(tree1, use.n = TRUE)
class.tree <- predict(tree1, newdata = cad2, type = "class")
classMat(cad2$CAD, class.tree)
cad3 <- subset(cad2, select = c("CAD", "AngPec", "Smoker"))

class.tree <- predict(tree1, newdata = cad3, type = "class")
class.bn <- predict(bn.back2, resp = "CAD", newdata = cad3, type = "class")



