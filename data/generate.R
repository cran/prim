
if (0)
{
quasiflow <- read.table("hiv-scaled.dat", header=TRUE)
quasiflow[,4] <- 2*quasiflow[,4]-1
colnames(quasiflow) <- c("x1", "x2", "x3", "y")
save(quasiflow, file="quasiflow.RData", ascii=TRUE, version=2)
}

