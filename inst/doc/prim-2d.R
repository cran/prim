### R code from vignette source 'prim-2d.Rnw'

###################################################
### code chunk number 1: prim-2d.Rnw:64-70
###################################################
library(prim)
library(MASS)
data(Boston)
x <- Boston[,5:6]  
y <- Boston[,1]  
boston.prim <- prim.box(x=x, y=y, threshold.type=1)


###################################################
### code chunk number 2: prim-2d.Rnw:100-101
###################################################
summary(boston.prim, print.box=TRUE)


###################################################
### code chunk number 3: prim-2d.Rnw:108-110
###################################################
plot(boston.prim, col="transparent")
points(x[y>3.5,])


###################################################
### code chunk number 4: prim-2d.Rnw:114-116
###################################################
plot(boston.prim, col="transparent")
points(x[y>3.5,])


###################################################
### code chunk number 5: prim-2d.Rnw:125-126
###################################################
boston.prim.med <- prim.box(x=x, y=y, threshold.type=1, y.fun=median)


###################################################
### code chunk number 6: prim-2d.Rnw:129-131
###################################################
plot(boston.prim, col="transparent")
plot(boston.prim.med, col="transparent", border="red", add=TRUE)


###################################################
### code chunk number 7: prim-2d.Rnw:135-138
###################################################
plot(boston.prim, col="transparent")
plot(boston.prim.med, col="transparent", border="red", add=TRUE)
legend("topleft", legend=c("mean", "median"), col=1:2, lty=1, bty="n")


###################################################
### code chunk number 8: prim-2d.Rnw:146-152
###################################################
x2 <- Boston[,c(5,9)]  
y <- Boston[,1]  
boston.cat.prim <- prim.box(x=x2, y=y, threshold.type=1)
summary(boston.cat.prim, print.box=TRUE)
plot(boston.cat.prim, col="transparent")
points(x2[y>3.5,])


###################################################
### code chunk number 9: prim-2d.Rnw:156-158
###################################################
plot(boston.cat.prim, col="transparent")
points(x2[y>3.5,])


