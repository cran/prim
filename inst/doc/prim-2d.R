## ---- echo=FALSE, message=FALSE-----------------------------------------------
knitr::opts_chunk$set(global.par=TRUE, collapse=TRUE, comment="#>", fig.width=5, fig.height=5, fig.align="center")
options(tibble.print_min=4L, tibble.print_max=4L)

## -----------------------------------------------------------------------------
library(prim)
library(MASS)
data(Boston)
x <- Boston[,5:6]
y <- Boston[,1]
boston.prim <- prim.box(x=x, y=y, threshold.type=1)

## -----------------------------------------------------------------------------
summary(boston.prim, print.box=TRUE)

## ---- fig.asp=1---------------------------------------------------------------
plot(boston.prim, col="transparent")
points(x[y>3.5,])

## -----------------------------------------------------------------------------
boston.prim.med <- prim.box(x=x, y=y, threshold.type=1, y.fun=median)

## ---- fig.asp=1---------------------------------------------------------------
plot(boston.prim, col="transparent")
plot(boston.prim.med, col="transparent", border="red", add=TRUE)
legend("topleft", legend=c("mean", "median"), col=1:2, lty=1, bty="n")

## ---- fig.asp=1---------------------------------------------------------------
x2 <- Boston[,c(5,9)]  
y <- Boston[,1]  
boston.cat.prim <- prim.box(x=x2, y=y, threshold.type=1)
summary(boston.cat.prim, print.box=TRUE)
plot(boston.cat.prim, col="transparent")
points(x2[y>3.5,])

