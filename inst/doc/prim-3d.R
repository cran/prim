### R code from vignette source 'prim-3d.Rnw'

###################################################
### code chunk number 1: prim-3d.Rnw:70-76
###################################################
library(prim)
data(quasiflow)
yflow <- quasiflow[,4]
xflow <- quasiflow[,1:3]
xflowp <- quasiflow[yflow==1,1:3]
xflown <- quasiflow[yflow==-1,1:3]


###################################################
### code chunk number 2: prim-3d.Rnw:80-82
###################################################
pairs(xflowp[1:300,])
pairs(xflown[1:300,])


###################################################
### code chunk number 3: prim-3d.Rnw:89-90
###################################################
pairs(xflowp[1:300,])


###################################################
### code chunk number 4: prim-3d.Rnw:93-94
###################################################
pairs(xflown[1:300,])


###################################################
### code chunk number 5: prim-3d.Rnw:104-106 (eval = FALSE)
###################################################
## qflow.thr <- c(0.38, -0.23)
## qflow.prim <- prim.box(x=xflow, y=yflow, threshold=qflow.thr, threshold.type=0)


###################################################
### code chunk number 6: prim-3d.Rnw:115-116
###################################################
qflow.hdr.pos <- prim.box(x=xflow, y=yflow, threshold=0.38, threshold.type=1)


###################################################
### code chunk number 7: prim-3d.Rnw:119-123
###################################################
qflow.neg <- prim.box(x=xflow, y=yflow, threshold.type=-1)
qflow.hdr.neg1 <- prim.hdr(qflow.neg, threshold=-0.23, threshold.type=-1)
qflow.hdr.neg2 <- prim.hdr(qflow.neg, threshold=-0.43, threshold.type=-1)
qflow.hdr.neg3 <- prim.hdr(qflow.neg, threshold=-0.63, threshold.type=-1)


###################################################
### code chunk number 8: prim-3d.Rnw:127-129
###################################################
qflow.prim2 <- prim.combine(qflow.hdr.pos, qflow.hdr.neg1)
summary(qflow.prim2)


###################################################
### code chunk number 9: prim-3d.Rnw:135-136
###################################################
plot(qflow.prim2)


