## ---- echo=FALSE, message=FALSE-----------------------------------------------
knitr::opts_chunk$set(collapse=TRUE, comment="#>", fig.width=5, fig.height=5, fig.align="center", global.par=TRUE, dpi=96) 
options(tibble.print_min=4L, tibble.print_max=4L)

## -----------------------------------------------------------------------------
library(prim)
data(quasiflow)
yflow <- quasiflow[,4]
xflow <- quasiflow[,1:3]
xflowp <- quasiflow[yflow==1,1:3]
xflown <- quasiflow[yflow==-1,1:3]

## -----------------------------------------------------------------------------
pairs(xflowp, cex=0.5, pch=16, col=grey(0,0.1), xlim=c(0,1), ylim=c(0,1))
pairs(xflown, cex=0.5, pch=16, col=grey(0,0.1), xlim=c(0,1), ylim=c(0,1))

## -----------------------------------------------------------------------------
qflow.hdr.pos <- prim.box(x=xflow, y=yflow, threshold=0.38, threshold.type=1)
summary(qflow.hdr.pos)

## -----------------------------------------------------------------------------
qflow.neg <- prim.box(x=xflow, y=yflow, threshold.type=-1)
qflow.hdr.neg1 <- prim.hdr(qflow.neg, threshold=-0.23, threshold.type=-1)
qflow.hdr.neg2 <- prim.hdr(qflow.neg, threshold=-0.43, threshold.type=-1)
qflow.hdr.neg3 <- prim.hdr(qflow.neg, threshold=-0.63, threshold.type=-1)

## -----------------------------------------------------------------------------
summary(qflow.hdr.neg1)

## -----------------------------------------------------------------------------
qflow.prim2 <- prim.combine(qflow.hdr.pos, qflow.hdr.neg1)
summary(qflow.prim2)

## ---- fig.asp=1, fig.height=10------------------------------------------------
plot(qflow.prim2, x.pt=xflow, pch=16, cex=0.5, alpha=0.1)

## ---- fig.asp=1, fig.height=10------------------------------------------------
plot(qflow.prim2, x.pt=xflow, pch=16, cex=0.5, alpha=0.1, splom=FALSE, colkey=FALSE, ticktype="detailed")

