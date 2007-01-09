
###############################################################################
#### PRIM (Patient rule induction method) for bump-hunting
###############################################################################

###############################################################################
## PRIM (Patient rule induction method)
##
## Parameters
## x - matrix of explanatory variables
## y - vector of response variable
## box.init - initial box (should cover range of x)
## mass.min - min. size of box mass
## y.mean.min - min. threshold of mean of y within a box
## pasting - TRUE - include pasting step (after peeling)
##         - FALSE - don't include pasting
##
## Returns
## list with k fields, one for each box
## each field is in turn is a list with fields
## x - data inside box
## y - corr. response values
## y.mean - mean of y
## box - limits of box
## box.mass - box mass
## num.boxes - total number of boxes with box mean >= y.mean.min
###############################################################################


prim.box <- function(x, y, box.init=NULL, peel.alpha=0.05, paste.alpha=0.01,
                 mass.min=0.05, threshold, pasting=TRUE, verbose=FALSE,
                 threshold.type=0)
{
  if (threshold.type==1 | threshold.type==-1)
  {
    if (missing(threshold))
      threshold <- mean(y)

    
    prim.reg <- prim.one(x=x, y=threshold.type*y, box.init=box.init,
                         peel.alpha=peel.alpha,
                         paste.alpha=paste.alpha, mass.min=mass.min[1],
                         threshold.type=threshold.type, 
                         threshold=threshold[1], pasting=pasting,
                         verbose=verbose)
  }
  else
  {
    if (missing(threshold))
      threshold <- c(mean(y), mean(y))
    else if (!missing(threshold))
      if (length(threshold)==1)
        stop("Need both upper and lower values for threshold")
      
    prim.plus <- prim.one(x=x, y=y, box.init=box.init, peel.alpha=peel.alpha,
                          paste.alpha=paste.alpha, mass.min=mass.min,
                          threshold.type=1,
                          threshold=threshold[1], pasting=pasting, verbose=verbose)
    prim.minus <- prim.one(x=x, y=-y, box.init=box.init, peel.alpha=peel.alpha,
                           paste.alpha=paste.alpha,mass.min=mass.min,
                           threshold.type=-1,
                           threshold=threshold[2], pasting=pasting, verbose=verbose)
    prim.reg <- prim.combine(prim.plus, prim.minus)
  }

  return(prim.reg)
    
}


prim.one <- function(x, y, box.init=NULL, peel.alpha=0.05, paste.alpha=0.01,
                     mass.min=0.05, threshold, pasting=TRUE, threshold.type=1,
                     verbose=FALSE)
{
  d <- ncol(x)
  n <- nrow(x)
  k.max <- ceiling(1/mass.min)
  num.boxes <- k.max 
  
  y.mean <- mean(y)
  mass.init <- length(y)/n 
  
  if (is.null(box.init))
  {   
    box.init <- apply(x, 2, range)
    box.init[1,] <- box.init[1,] - 0.1*abs(diff(box.init))
    box.init[2,] <- box.init[2,] + 0.1*abs(diff(box.init))
  }

  ## find first box
  k <- 1
  boxk <- find.box(x=x, y=y, box=box.init, peel.alpha=peel.alpha,
                   paste.alpha=paste.alpha, mass.min=mass.min,
                   threshold=mean(y), d=d, n=n,
                   pasting=pasting, verbose=verbose)

  if (is.null(boxk))
  {
    if (verbose)
      warning(paste("Unable to find box", k, "\n"))
    return(list(list(x=x, y=threshold.type*y, y.mean=threshold.type*y.mean, box=box.init, box.mass=mass.init)))
  }
  else
  {
    if (verbose)
      cat(paste("Found find box ", k, ": y.mean=", signif(threshold.type*boxk$y.mean,4), ", mass=",
                signif(boxk$mass,4), "\n\n", sep=""))
    boxes <- list(x=list(boxk$x), y=list(boxk$y), y.mean=list(boxk$y.mean),
                  box=list(boxk$box), mass=list(boxk$mass))       
  }
 
  ## find subsequent boxes
  if (num.boxes > 1)
  {
    boxk <- list(x=boxes$x[[k]], y=boxes$y[[k]], y.mean=boxes$y.mean[[k]],
                 box=boxes$box[[k]], mass=boxes$mass[[k]])
    ## data still under consideration
    
    x.out.ind.mat <-  matrix(TRUE, nrow=nrow(x), ncol=ncol(x))
   
    for (j in 1:d)
      x.out.ind.mat[,j] <- (x[,j] < boxk$box[1,j]) | (x[,j] > boxk$box[2,j])

    x.out.ind <- apply(x.out.ind.mat, 1, sum)!=0
    
    x.out <- x[x.out.ind,]
    y.out <- y[x.out.ind]
    box.out <- apply(x.out, 2, range)
   
    while ((k < num.boxes) & (!is.null(boxk))) 
    {
      k <- k+1

      boxk <- find.box(x=x.out, y=y.out, box=box.out,
                       peel.alpha=peel.alpha, paste.alpha=paste.alpha,
                       mass.min=mass.min, threshold=min(y), d=d, n=n,
                       pasting=pasting, verbose=verbose)
      
      if (is.null(boxk))
      {
        if (verbose)
          cat(paste("Bump", k, "includes all remaining data\n\n"))

        boxes$x[[k]] <- x.out
        boxes$y[[k]] <- y.out
        boxes$y.mean[[k]] <- mean(y.out)
        boxes$box[[k]] <- box.init
        boxes$mass[[k]] <- length(y.out)/n
      }
      else 
      {
        ## update x and y
        if (verbose)
          cat(paste("Found find box ", k, ": y.mean=", signif(threshold.type*boxk$y.mean,4),
                    ", mass=", signif(boxk$mass,4), "\n\n", sep=""))
        
        x.out.ind.mat <- matrix(TRUE, nrow=nrow(x), ncol=ncol(x))
        for (j in 1:d)
          x.out.ind.mat[,j] <- (x[,j] < boxk$box[1,j]) | (x[,j] > boxk$box[2,j])
        
        x.out.ind <- x.out.ind & (apply(x.out.ind.mat, 1, sum)!=0)
        x.out <- x[x.out.ind,]
        y.out <- y[x.out.ind]
     
        boxes$x[[k]] <- boxk$x
        boxes$y[[k]] <- boxk$y
        boxes$y.mean[[k]] <- boxk$y.mean
        boxes$box[[k]] <- boxk$box
        boxes$mass[[k]] <-boxk$mass 
      }
    }   
  }

  ## adjust for negative hdr  
  for (k in 1:length(boxes$y.mean))
  {
    boxes$y[[k]] <- threshold.type*boxes$y[[k]]
    boxes$y.mean[[k]] <- threshold.type*boxes$y.mean[[k]]
  }
  
  ## highest density region

  prim.res <- prim.hdr(prim=boxes, threshold=threshold, threshold.type=threshold.type)
  
  return(prim.res)
         
}

###############################################################################
## Highest density region for PRIM boxes
###############################################################################

prim.hdr <- function(prim, threshold, threshold.type)
{  
  n <- 0
  for (i in 1:length(prim$box))
    n <- n + length(prim$y[[i]])

  excess.ind <- which(unlist(prim$y.mean)*threshold.type >= threshold*threshold.type)

  if (length(excess.ind) > 0)
    excess.ind <- max(excess.ind)
  else
  {
    if (threshold.type==1)
      warning(paste("No prim box found with mean >=", threshold))
    else if (threshold.type==-1)
      warning(paste("No prim box found with mean <=", threshold))
    return()
  }
  
  ## highest density region  
  x.prim.hdr <- list()
  
  for (k in 1:excess.ind)
  {
    x.prim.hdr$x[[k]] <- prim$x[[k]]
    x.prim.hdr$y[[k]] <- prim$y[[k]]
    x.prim.hdr$y.mean[[k]] <- prim$y.mean[[k]]
    x.prim.hdr$box[[k]] <- prim$box[[k]]
    x.prim.hdr$mass[[k]] <-prim$mass[[k]] 
  }
  
  ## combine non-excess into a `dump' box
  if (excess.ind < length(prim$x))
  {
    x.temp <- numeric()
    y.temp <- numeric()
    for (k in (excess.ind+1):length(prim$x))
    {
      x.temp <- rbind(x.temp, prim$x[[k]])
      y.temp <- c(y.temp, prim$y[[k]])
    }
    
    x.prim.hdr$x[[excess.ind+1]] <- x.temp
    x.prim.hdr$y[[excess.ind+1]] <- y.temp
    x.prim.hdr$y.mean[[excess.ind+1]] <- mean(y.temp)
    x.prim.hdr$box[[excess.ind+1]] <- prim$box[[length(prim$x)]]
    x.prim.hdr$mass[[excess.ind+1]] <- length(y.temp)/n  
    }
  
  x.prim.hdr$num.class <- length(x.prim.hdr$x)
  x.prim.hdr$num.hdr.class <- excess.ind
  x.prim.hdr$threshold <- threshold
  
  x.prim.hdr$ind <- rep(threshold.type, x.prim.hdr$num.hdr.class)
  
  class(x.prim.hdr) <- "prim"
  
  return(x.prim.hdr)
    
}
 
                        
###############################################################################
## Combine (disjoint) PRIM box sequences - useful for joining
## positive and negative estimates
##
## Parameters
## prim1 - 1st PRIM box sequence
## prim2 - 2nd PRIM box sequence
##
## Returns
## same as for prim()
###############################################################################

prim.combine <- function(prim1, prim2)
{
  M1 <- prim1$num.hdr.class
  M2 <- prim2$num.hdr.class

  if (is.null(M1) & !is.null(M2))
    return (prim2)
  if (!is.null(M1) & is.null(M2))
    return(prim1)
  if (is.null(M1) & is.null(M2))
    return(NULL)
    
  overlap <- overlap.box.seq(prim1, prim2)

  x <- numeric()
  y <- vector()
  for (i in 1:prim1$num.class)
  {
    x <- rbind(x, prim1$x[[i]])
    y <- c(y, prim1$y[[i]])
  }  

  if (any(overlap[1:M1,1:M2]))
  {
    warning("Class boundaries overlap - will return NULL")
    return(NULL)
  }
  else
  {
    prim.temp <- list()
    for (i in 1:M1)
    {
      prim.temp$x[[i]] <- prim1$x[[i]]
      prim.temp$y[[i]] <- prim1$y[[i]]
      prim.temp$y.mean[[i]] <- prim1$y.mean[[i]]
      prim.temp$box[[i]] <- prim1$box[[i]]
      prim.temp$mass[[i]] <- prim1$mass[[i]]
      prim.temp$ind[[i]] <- 1
    }
    for (i in 1:M2)
    {
      prim.temp$x[[i+M1]] <- prim2$x[[i]]
      prim.temp$y[[i+M1]] <- prim2$y[[i]]
      prim.temp$y.mean[[i+M1]] <- prim2$y.mean[[i]]
      prim.temp$box[[i+M1]] <- prim2$box[[i]]
      prim.temp$mass[[i+M1]] <- prim2$mass[[i]]
      prim.temp$ind[[i+M1]] <- -1
    }
    
    dumpx.ind <- which.box(x, prim1)==prim1$num.class & which.box(x, prim2)==prim2$num.class
    
    prim.temp$x[[M1+M2+1]] <- x[dumpx.ind,]
    prim.temp$y[[M1+M2+1]] <- y[dumpx.ind]
    prim.temp$y.mean[[M1+M2+1]] <- mean(y[dumpx.ind])
    prim.temp$box[[M1+M2+1]] <- prim1$box[[prim1$num.class]]
    prim.temp$mass[[M1+M2+1]] <- length(y[dumpx.ind])/length(y)
    prim.temp$num.class <- M1+M2+1
    prim.temp$num.hdr.class <- M1+M2
    prim.temp$threshold <- c(prim1$threshold, prim2$threshold) 

    class(prim.temp) <- "prim"
  }

  return (prim.temp)
}
   
###############################################################################
## Finds box
##
## Parameters
## x - matrix of explanatory variables
## y - vector of response variable
## box.init - initial box (should cover range of x)
## mass.min - min box mass
## threshold - min box mean
## pasting - TRUE - include pasting step (after peeling)
##         - FALSE - don't include pasting
##
## Returns
## List with fields
## x - data still inside box after peeling
## y - corresponding response values
## y.mean - mean of y
## box - box limits
## mass - box mass
###############################################################################
find.box <- function(x, y, box, peel.alpha, paste.alpha, mass.min, threshold,
                     d, n, pasting, verbose) 
{
  y.mean <- mean(y)
  mass <- length(y)/n
  
  if ((y.mean >= threshold) & (mass >= mass.min))
    boxk.peel <- peel.one(x=x, y=y, box=box, peel.alpha=peel.alpha,
                           mass.min=mass.min,threshold=threshold, d=d, n=n)   
  else
    boxk.peel <- NULL

  boxk.temp <- NULL
 
  while (!is.null(boxk.peel))
  { 
    boxk.temp <- boxk.peel
    boxk.peel <- peel.one(x=boxk.temp$x, y=boxk.temp$y, box=boxk.temp$box,
                           peel.alpha=peel.alpha,
                           mass.min=mass.min, threshold=threshold, d=d, n=n)
  }
  
  if (verbose)
    cat("Peeling completed \n")
 
  if (pasting)
  {
    boxk.paste <- boxk.temp
    
    while (!is.null(boxk.paste))
    {
      
      boxk.temp <- boxk.paste
      boxk.paste <- paste.one(x=boxk.temp$x, y=boxk.temp$y, box=boxk.temp$box,
                              x.init=x, y.init=y, paste.alpha=paste.alpha,
                              mass.min=mass.min, threshold=threshold, d=d, n=n)      
    }
    if (verbose)
      cat("Pasting completed\n")  
  }
   
  boxk <- boxk.temp

  return(boxk)
}



############################################################################### 
## Peeling stage of PRIM
##
## Parameters
## x - data matrix
## y - vector of response variables
## peel.alpha - peeling quantile
## paste.alpha - peeling proportion
## mass.min - minimum box mass
## threshold - minimum y mean
## d - dimension of data
## n - number of data
## 
## Returns
## List with fields
## x - data still inside box after peeling
## y - corresponding response values
## y.mean - mean of y
## box - box limits
## mass - box mass
###############################################################################

peel.one <- function(x, y, box, peel.alpha, mass.min, threshold, d, n, type=5)
{
  box.new <- box
  mass <- length(y)/n
  y.mean <- mean(y)
  y.mean.peel <- matrix(0, nrow=2, ncol=d)
  box.vol.peel <- matrix(0, nrow=2, ncol=d) 
                     
  for (j in 1:d)
  {
    box.min.new <- quantile(x[,j], peel.alpha, type=type)
    box.max.new <- quantile(x[,j], 1-peel.alpha, type=type)
    
    y.mean.peel[1,j] <- mean(y[x[,j] >= box.min.new])
    y.mean.peel[2,j] <- mean(y[x[,j] <= box.max.new])
    
    box.temp1 <- box
    box.temp2 <- box
    box.temp1[1,j] <- box.min.new
    box.temp2[2,j] <- box.max.new
    box.vol.peel[1,j] <- vol.box(box.temp1)
    box.vol.peel[2,j] <- vol.box(box.temp2)    
  }
  
  y.mean.peel.max.ind <- which(y.mean.peel==max(y.mean.peel), arr.ind=TRUE)

  ## break ties by choosing box with largest volume

  nrr <- nrow(y.mean.peel.max.ind) 
  if (nrr > 1)
  {
    box.vol.peel2 <- rep(0, nrr)
    for (j in 1:nrr)
      box.vol.peel2[j] <- box.vol.peel[y.mean.peel.max.ind[j,1],
                                       y.mean.peel.max.ind[j,2]]
    
    row.ind <- which(max(box.vol.peel2)==box.vol.peel2)
  }
  else
    row.ind <- 1
  
  y.mean.peel.max.ind <- y.mean.peel.max.ind[row.ind,]
  ## peel along dimension j.max
  j.max <- y.mean.peel.max.ind[2]

  ## peel lower 
  if (y.mean.peel.max.ind[1]==1)
  {
    box.new[1,j.max] <- quantile(x[,j.max], peel.alpha, type=type)
    x.index <- x[,j.max] >= box.new[1,j.max] 
  }
  ## peel upper 
  else if (y.mean.peel.max.ind[1]==2)
  {
    box.new[2,j.max] <- quantile(x[,j.max], 1-peel.alpha, type=type)
    x.index <- x[,j.max] <= box.new[2,j.max]
  }
 
  x.new <- x[x.index,]
  y.new <- y[x.index]
  mass.new <- length(y.new)/n
  y.mean.new <- mean(y.new)
  
  ## if min. y mean and min. mass conditions are still true, update
  ## o/w return NULL  

  if ((y.mean.new >= threshold) & (mass.new >= mass.min) & (mass.new < mass))
    return(list(x=x.new, y=y.new, y.mean=y.mean.new, box=box.new,
                mass=mass.new))
}


###############################################################################
## Pasting stage for PRIM
##
## Parameters
## x - data matrix
## y - vector of response variables
## x.init - initial data matrix (superset of x) 
## y.init - initial response vector (superset of y) 
## peel.alpha - peeling quantile
## paste.alpha - peeling proportion
## mass.min - minimum box mass
## threshold - minimum y mean
## d - dimension of data
## n - number of data
## 
## Returns
##
## List with fields
## x - data still inside box after peeling
## y - corresponding response values
## y.mean - mean of y
## box - box limits
## box.mass - box mass
###############################################################################

paste.one <- function(x, y, x.init, y.init, box, paste.alpha,
                      mass.min, threshold, d, n)
{
  box.new <- box
  mass <- length(y)/n
  y.mean <- mean(y)
  n.box <- length(y)
    
  y.mean.paste <- matrix(0, nrow=2, ncol=d)
  box.vol.paste <- matrix(0, nrow=2, ncol=d)
  box.paste <- matrix(0, nrow=2, ncol=d)
  x.paste1.list <- list()
  x.paste2.list <- list()
  y.paste1.list <- list()
  y.paste2.list <- list()

  for (j in 1:d)
  {    
    x.min <- min(x[,j])
    x.max <- max(x[,j])

    ## candidates for pasting
    x.res.ind <- in.box.j(x=x.init, box=box, j=j, d=d, n=nrow(x.init))
    
    x.res <- x.init[x.res.ind,]
    y.res <- y.init[x.res.ind]
    x.cand1 <- x.res[x.res[,j] < x.min,j]
    x.cand2 <- x.res[x.res[,j] > x.max,j]
    y.cand1 <- y.res[x.res[,j] < x.min]
    y.cand2 <- y.res[x.res[,j] > x.max]
    
    ## paste below minimum
    k <- 1
  
    n.paste1 <- ceiling(paste.alpha*n.box)
    x.paste.ind1 <- 1:min(n.paste1, length(x.cand1))
    x.paste1 <- c(x[,j], x.cand1[order(x.cand1, decreasing=TRUE)][x.paste.ind1])
    y.paste1 <- c(y, y.cand1[order(x.cand1, decreasing=TRUE)][x.paste.ind1])    
   
    x.paste1 <- na.omit(x.paste1)
    y.paste1 <- na.omit(y.paste1)
    x.paste1.list[[j]] <- x.paste1
    y.paste1.list[[j]] <- y.paste1
    
    ## paste above maximum
    k <- 1
    
    n.paste2 <- ceiling(paste.alpha*n.box)
    x.paste.ind2 <- 1:min(n.paste2, length(x.cand2))
    x.paste2 <- c(x[,j], x.cand2[order(x.cand2)][x.paste.ind2])
    y.paste2 <- c(y, y.cand2[order(x.cand2)][x.paste.ind2])
    
    x.paste2 <- na.omit(x.paste2)
    y.paste2 <- na.omit(y.paste2)
    x.paste2.list[[j]] <- x.paste2
    y.paste2.list[[j]] <- y.paste2

    ## y means of pasted boxes
    y.mean.paste[1,j] <- mean(y.paste1)
    y.mean.paste[2,j] <- mean(y.paste2)

    box.temp1 <- box
    box.temp2 <- box
    box.temp1[1,j] <- min(x.paste1)
    box.temp2[2,j] <- max(x.paste2)

    ## volume of pasted boxes
    box.vol.paste[1,j] <- vol.box(box.temp1)
    box.vol.paste[2,j] <- vol.box(box.temp2)

    ## limits of pasted boxes
    box.paste[1,j] <- min(x.paste1) 
    box.paste[2,j] <- max(x.paste2)
  }
  
  y.mean.paste.max.ind <- which(y.mean.paste==max(y.mean.paste), arr.ind=TRUE)
  
  ## break ties by choosing box with largest volume

  nrr <- nrow(y.mean.paste.max.ind) 
  if (nrr > 1)
  {
    box.vol.paste2 <- rep(0, nrr)
    for (j in 1:nrr)
      box.vol.paste2[j] <- box.vol.paste[y.mean.paste.max.ind[j,1],
                                         y.mean.paste.max.ind[j,2]]
    row.ind <- which(max(box.vol.paste2)==box.vol.paste2)
  }
  else
    row.ind <- 1
  
  y.mean.paste.max.ind <- y.mean.paste.max.ind[row.ind,]

  ## paste along dimension jmax
  j.max <- y.mean.paste.max.ind[2]

  ## paste lower 
  if (y.mean.paste.max.ind[1]==1)
  { 
    box.new[1,j.max] <- box.paste[1,j.max]
    x.index <- match(x.paste1.list[[j.max]], x.init[,j.max])
  }
  
  ## paste upper
  else if (y.mean.paste.max.ind[1]==2)
  {
    box.new[2,j.max] <- box.paste[2,j.max]
    x.index <- match(x.paste2.list[[j.max]], x.init[,j.max])
  }

  x.new <- x.init[x.index,]
  y.new <- y.init[x.index]
  mass.new <- length(y.new)/n
  y.mean.new <- mean(y.new)

  if ((y.mean.new > threshold) & (mass.new >= mass.min) & (y.mean.new >= y.mean)
      & (mass.new > mass))
    return(list(x=x.new, y=y.new, y.mean=y.mean.new, box=box.new, mass=mass.new))
 
}




###############################################################################
## Output functions for prim objects
###############################################################################

###############################################################################
## Plot function for PRIM objects
###############################################################################


plot.prim <- function(x, ...)
{
  if (ncol(x$x[[1]])==2)
    plotprim.2d(x, ...)
  else if (ncol(x$x[[1]])==3)
    plotprim.3d(x, ...)
  else if (ncol(x$x[[1]])>3)
    plotprim.nd(x, ...)
  
  invisible()
}

plotprim.2d <- function(x, col, xlim, ylim, xlab, ylab, add=FALSE,
   add.legend=FALSE, cex.legend=1, pos.legend, lwd=1, ...)
{ 
  M <- x$num.hdr.class 
  
  ff <- function(x, d) { return (x[,d]) }
   
  if (missing(xlim))
    xlim <- range(sapply(x$box, ff, 1))
  if (missing(ylim))
    ylim <- range(sapply(x$box, ff, 2))

  x.names <- colnames(x$x[[1]])
  if (is.null(x.names)) x.names <- c("x","y")
  
  if (missing(xlab)) xlab <- x.names[1]
  if (missing(ylab)) ylab <- x.names[2] 
  
  if (!add)
    plot(x$box[[1]], type="n", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
   
  text.legend <- paste("box", 1:M, sep="")
  
  if (missing(pos.legend))
  {  
    pos.legend <- c(xlim[2], ylim[2])
    xlim <- c(xlim[1], xlim[2] + 0.1*abs(xlim[1]-xlim[2]))
  }

  if (missing(col))
    col <- topo.colors(M)
  if (length(col) < M)
    col <- rep(col, length=M)
    
  for (i in M:1)
  {
    ## colour i-th box
    box <- x$box[[i]]
    rect(box[1,1], box[1,2], box[2,1], box[2,2], border=TRUE, col=col[i], lwd=lwd)
  }
  
  if (add.legend)
    legend(pos.legend[1], pos.legend[2], legend=text.legend, fill=col, bty="n",
           cex=cex.legend)
  
  invisible()
}


plotprim.3d <- function(x, color, xlim, ylim, zlim, xlab, ylab, zlab, add.axis=TRUE,
                        ...)
{ 
  clear3d()
  rgl.bg(color="white")
  
  M <- x$num.hdr.class
   
  ff <- function(x, d) { return (x[,d]) }
   
  if (missing(xlim))
    xlim <- range(sapply(x$box, ff, 1))
  if (missing(ylim))
    ylim <- range(sapply(x$box, ff, 2))
  if (missing(zlim))
    zlim <- range(sapply(x$box, ff, 3))

  x.names <- colnames(x$x[[1]])
  if (is.null(x.names)) x.names <- c("x","y","z")
  
  if (missing(xlab)) xlab <- x.names[1]
  if (missing(ylab)) ylab <- x.names[2] 
  if (missing(zlab)) zlab <- x.names[3]

  if (add.axis)
  {
    lines3d(xlim[1:2], rep(ylim[1],2), rep(zlim[1],2), size=3, color="black")
    lines3d(rep(xlim[1],2), ylim[1:2], rep(zlim[1],2), size=3, color="black")
    lines3d(rep(xlim[1],2), rep(ylim[1],2), zlim[1:2], size=3, color="black")
  
    texts3d(xlim[2],ylim[1],zlim[1],xlab,size=3,color="black", adj=0)
    texts3d(xlim[1],ylim[2],zlim[1],ylab,size=3,color="black", adj=1)
    texts3d(xlim[1],ylim[1],zlim[2],zlab,size=3,color="black", adj=1)
  }
  
  if (missing(color))
    color <- topo.colors(M)
  if (length(color) < M)
    color <- rep(color, length=M)
  
  for (i in M:1)
  {
    ## colour data in i-th box
    xdata <- x$x[[i]]
    if (is.vector(xdata))
      points3d(xdata[1], xdata[2], xdata[3], color=color[i])
    else
      points3d(xdata[,1], xdata[,2], xdata[,3], color=color[i])
  }
  invisible()
}



plotprim.nd  <- function(x, col, xmin, xmax, xlab, ylab, ...)
{
  M <- x$num.hdr.class
  M2 <- x$num.class
  
  x.names <- colnames(x$x[[1]])
  if (is.null(x.names)) x.names <- c("x","y")
  
  if (missing(xlab)) xlab <- x.names[1]
  if (missing(ylab)) ylab <- x.names[2]
  
  if (missing(col))
    col <- topo.colors(M)
  if (length(col) < M)
    col <- rep(col, length=M)
  
  if (missing(xmin)) xmin <- x$box[[M2]][1,]
  if (missing(xmax)) xmax <- x$box[[M2]][2,]

  xdata <- numeric()
  xcol <- numeric() 
  for (i in M:1)
  {
    xdata <- rbind(xdata, x$x[[i]])
    xcol <- c(xcol, rep(col[i], nrow(as.matrix(x$x[[i]])))) 
  }
 
  pairs(rbind(xmin, xmax, xdata), col=c("transparent", "transparent",xcol), ...)
  invisible()
}


###############################################################################
## Summary function for PRIM objects
##
## Parameters
## x - prim object
##
## Returns
## matrix, 
## with 2 columns y.mean (mean of y) and mass (box mass)
## i-th row corresponds to i-th box of x
## last row is overall y mean and total mass covered by boxes 
###############################################################################

summary.prim <- function(object, ...)
{
  x <- object
  M <- x$num.class

  summ.mat <- vector()
  for (k in 1:M)
    summ.mat <- rbind(summ.mat, c(x$y.mean[[k]], x$mass[[k]], x$ind[k]))

  tot <- c(sum(summ.mat[,1]*summ.mat[,2])/sum(summ.mat[,2]), sum(summ.mat[,2]), NA)
  summ.mat <- rbind(summ.mat, tot)
  
  rownames(summ.mat) <- c(paste("box", 1:(nrow(summ.mat)-1), sep=""), "overall")
  colnames(summ.mat) <- c("box-mean", "box-mass", "box-ind")

  if (x$num.hdr.class < x$num.class)
    for (k in (x$num.hdr.class+1):x$num.class)
      rownames(summ.mat)[k] <- paste(rownames(summ.mat)[k], "*",sep="")

  print(summ.mat)
  
  if (x$num.hdr.class < x$num.class)
    cat("\n* - box not in HDR at level =", x$threshold,"\n\n")

  for (k in 1:M)
  {
    cat(paste("Box limits for box", k, "\n", sep=""))
    box.summ <- x$box[[k]]
    rownames(box.summ) <- c("min", "max")
    print(box.summ)
    cat("\n")
  }
}



prim.thresh.symdiff.mixt <- function(prim.plus, prim.minus, mus1, Sigmas1, props1, mus2, Sigmas2, props2, weight.dd=1/2, weight.mixt=1/2, nmc=1e6, alpha=0.5, xmc.dd, taumc.dd, gmc.mixt.dd, xmc.mixt, threshold.plus.range, threshold.minus.range, verbose=FALSE)
{
  if (missing(xmc.dd))
    xmc.dd <- rmvnorm.mixt.dd(n=nmc, mus1, Sigmas1, props1, mus2, Sigmas2, props2)

  if (missing(taumc.dd))
    taumc.dd <- thresh.mixt.dd.mc(mus1=mus1, Sigmas1=Sigmas1, props1=props1, mus2=mus2, Sigmas2=Sigmas2, props2=props2, weight.dd=weight.dd, xmc=xmc.dd, alpha=alpha, unit.trans=FALSE)
  
  if (missing(xmc.mixt))
    xmc.mixt <- rmvnorm.mixt.mixt(n=nmc, mus1=mus1, Sigmas1=Sigmas1, props1=props1, mus2=mus2, Sigmas2=Sigmas2, props2=props2, weight.mixt=weight.mixt)
  
  if (missing(gmc.mixt.dd))
    (gmc.mixt.dd) <- dmvnorm.mixt.dd(x=xmc.mixt, mus1=mus1, Sigmas1=Sigmas1, props1=props1, mus2=mus2, Sigmas2=Sigmas2, props2=props2, weight.dd=weight.dd)

  xinR.plus <- gmc.mixt.dd >= taumc.dd[1]
  xinR.minus <- gmc.mixt.dd <= taumc.dd[2]

  ## default threshold ranges
  if (missing(threshold.plus.range))
  {
    if (prim.plus$num.hdr.class >1)
    {
      temp.range <- prim.plus$y.mean[1:prim.plus$num.hdr.class]
      threshold.plus.range <- trunc(seq(min(temp.range), max(temp.range), length=11)*100)/100
    }
    else
      threshold.plus.range <- trunc(prim.plus$y.mean[1]*100)/100
        
  }
  
  if (missing(threshold.minus.range))
  {
    if (prim.minus$num.hdr.class >1)
    {
      temp.range <- prim.minus$y.mean[1:prim.minus$num.hdr.class] 
      threshold.minus.range <- trunc(seq(min(temp.range), max(temp.range), length=11)*100)/100
    }
    else
      threshold.minus.range <- trunc(prim.minus$y.mean[1]*100)/100
  }

  
  ## step through threshold ranges to find min. error
  err.plus <- 0
  if (!missing(prim.plus))
  {  
    err.plus <- vector()
    i <- 0
    for (thp in threshold.plus.range)
    {
      prim.plus.hdr <- prim.hdr(prim=prim.plus, threshold=thp, threshold.type=1)
      
      xmc.mixt.class <- which.box(x=xmc.mixt, box.seq=prim.plus.hdr)
      xmc.mixt.class[xmc.mixt.class<=length(prim.plus.hdr$ind)] <- 1
      xmc.mixt.class[xmc.mixt.class>length(prim.plus.hdr$ind)] <- 0
      
      xinRhat <- xmc.mixt.class
      xinRhat.plus <- xinRhat==1
      symdiff.plus <- (xinR.plus & !xinRhat.plus) | (!xinR.plus & xinRhat.plus)
      i <- i+1
      err.plus[i] <- mean(symdiff.plus)    
    }
    
    thresh.plus <- rev(threshold.plus.range)[which.min(rev(err.plus))]
  }
    
 
  err.minus <- 0
  if (!missing(prim.minus))
  {
    err.minus <- vector()
    i <- 0
    
    for (thm in threshold.minus.range)
    {
      prim.minus.hdr <- prim.hdr(prim=prim.minus, threshold=thm, threshold.type=-1)

      if (!is.null(prim.minus.hdr))
      {
        xmc.mixt.class <- which.box(x=xmc.mixt, box.seq=prim.minus.hdr)
        xmc.mixt.class[xmc.mixt.class<=length(prim.minus.hdr$ind)] <- -1
        xmc.mixt.class[xmc.mixt.class>length(prim.minus.hdr$ind)] <- 0
        
        xinRhat <- xmc.mixt.class
        xinRhat.minus <- xinRhat==-1
        symdiff.minus <- (xinR.minus & !xinRhat.minus) | (!xinR.minus & xinRhat.minus)
        i <- i+1
        err.minus[i] <- mean(symdiff.minus)}
      else
      {
        i <- i+1
        err.minus[i] <- 1
      }
    }
  
    thresh.minus <- threshold.minus.range[which.min(err.minus)]
  }

  err <- min(err.plus) + min(err.minus)

  if (verbose)
  {
    tab.plus <- rbind(threshold.plus.range, err.plus)
    rownames(tab.plus) <- c("thresh.plus", "err.plus")
    tab.minus <- rbind(threshold.minus.range, err.minus)
    rownames(tab.minus) <- c("thresh.minus", "err.minus")
    print(tab.plus)
    print(tab.minus)
  }
    
  if (missing(prim.minus) & !missing(prim.plus))
    thresh.val <- thresh.plus
  if (missing(prim.plus) & !missing(prim.minus))
    thresh.val <- thresh.minus
  if (!missing(prim.plus) & !missing(prim.minus))
    thresh.val <- c(thresh.plus, thresh.minus) 
  if (missing(prim.plus) & missing(prim.minus))
    thresh.val <- NULL

  return(thresh.val)
}

