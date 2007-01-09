
if (0)
{
load("/c/tduong/home/research/flowcyt/densdiff/library/densdiff/data/hiv02.RData")
##names(hiv)[4:8] <- c("FSC", "SSC", "CD3", "CD8", "CD4")

n <- 1e4
hiv.stat <- hiv$hiv.stat
hiv.flow <- hiv[,4:8]
names(hiv.flow) <- c("x1", "x2", "x3", "x4", "x5")
hiv.flow.neg <- hiv.flow[hiv.stat==-1,]
hiv.flow.pos <- hiv.flow[hiv.stat==1,]

## fit 4-component mixture normal to data

ncent <- 4  ## 4 obtained from hierarchical clustering
kmeans.pos <- kmeans(hiv.flow.pos, centers=ncent)
kmeans.neg <- kmeans(hiv.flow.neg, centers=ncent)

mus.m <- round(kmeans.neg$centers,0)
mus.p <- round(kmeans.pos$centers,0)

Sigmas.m <- numeric()
for (i in 1:ncent)
  Sigmas.m <- rbind(Sigmas.m, var(hiv.flow.neg[kmeans.neg$cluster==i,]))

Sigmas.p <- numeric()
for (i in 1:ncent)
  Sigmas.p <- rbind(Sigmas.p, var(hiv.flow.pos[kmeans.pos$cluster==i,]))

for (i in 1:ncent)
print(vech(round(var(hiv.flow.neg[kmeans.neg$cluster==i,]),0)))

Sigmas.m <- round(Sigmas.m, 0)
Sigmas.p <- round(Sigmas.p, 0)

props.m <- round(kmeans.neg$size/sum(kmeans.neg$size), 2)
props.p <- round(kmeans.pos$size/sum(kmeans.pos$size), 2)

set.seed(8192)
x.p <- rmvnorm.mixt(n=n, mu=mus.p, Sigmas=Sigmas.p, props=props.p)
colnames(x.p) <- names(hiv.flow)
x.m <- rmvnorm.mixt(n=n, mu=mus.m, Sigmas=Sigmas.m, props=props.m)
colnames(x.m) <- names(hiv.flow)

quasiflow <- round(rbind(x.p[1:1000,], x.m[1:1000,]),1)
y <- c(rep(1, 1000), rep(-1, 1000))
quasiflow <- cbind(quasiflow, y)
}

