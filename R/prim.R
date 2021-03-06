
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
                 threshold.type=0, y.fun=mean)
{
  if (threshold.type==1 | threshold.type==-1)
  {
    if (missing(threshold)) threshold <- do.call(y.fun, list(x=y))

    prim.temp <- prim.one(x=x, y=threshold.type*y, box.init=box.init,
                         peel.alpha=peel.alpha,
                         paste.alpha=paste.alpha, mass.min=mass.min,
                         threshold.type=threshold.type, 
                         threshold=threshold[1], pasting=pasting,
                         verbose=verbose, y.fun=y.fun)
  }
  else
  {
    if (missing(threshold))
      threshold <- rep(do.call(y.fun, list(x=y)), 2) ##c(mean(y), mean(y))
    else if (!missing(threshold))
      if (length(threshold)==1)
        stop("Need both upper and lower values for threshold")	
	
    prim.pos <- prim.one(x=x, y=y, box.init=box.init, peel.alpha=peel.alpha,
                         paste.alpha=paste.alpha, mass.min=mass.min,
                         threshold.type=1, threshold=threshold[1],
                         pasting=pasting, verbose=verbose, y.fun=y.fun)
    prim.neg <- prim.one(x=x, y=-y, box.init=box.init, peel.alpha=peel.alpha,
                         paste.alpha=paste.alpha,mass.min=mass.min,
                         threshold.type=-1, threshold=threshold[2],
                         pasting=pasting, verbose=verbose, y.fun=y.fun)
    prim.temp <- prim.combine(prim.pos, prim.neg, y.fun=y.fun)
    
  }
  
  ## re-do prim to ensure that no data points are missed from the `dump' box 
  prim.reg <- prim.temp
  prim.labels <- prim.which.box(x=x, box.seq=prim.reg)
  for (k in 1:prim.reg$num.class)
  {
    primk.ind <- which(prim.labels==k)
    prim.reg$x[[k]] <- x[primk.ind,]
    prim.reg$y[[k]] <- y[primk.ind]
    prim.reg$y.fun[k] <- do.call(y.fun, list(x=prim.reg$y[[k]]))
    prim.reg$mass[k] <- length(prim.reg$y[[k]])/nrow(x)
  }
  
  return(prim.reg)
}


prim.one <- function(x, y, box.init=NULL, peel.alpha=0.05, paste.alpha=0.01,
                     mass.min=0.05, threshold, pasting=FALSE, threshold.type=1,
                     verbose=FALSE, y.fun=mean)
{
  d <- ncol(x)
  n <- nrow(x)
  k.max <- ceiling(1/mass.min)
  num.boxes <- k.max 

  ##if (is.vector(x)) x <- as.matrix(t(x))
  y.fun.val <- do.call(y.fun, list(x=y))
  mass.init <- length(y)/n 
  
  if (is.null(box.init))
  {
    box.init <- apply(x, 2, range)
    box.diff <- box.init[2,] - box.init[1,]
    box.init[1,] <- box.init[1,] - 10*paste.alpha*box.diff
    box.init[2,] <- box.init[2,] + 10*paste.alpha*box.diff
  }
  
  ## find first box
  k <- 1
  
  boxk <- find.box(x=x, y=y, box=box.init, peel.alpha=peel.alpha,
                   paste.alpha=paste.alpha, mass.min=mass.min,
                   threshold=min(y)-0.1*abs(min(y)), d=d, n=n, pasting=pasting, verbose=verbose, y.fun=y.fun)
		 
  if (is.null(boxk))
  {
    if (verbose)
      warning(paste("Unable to find box", k, "\n"))

    x.prim <- list(x=list(x), y=list(threshold.type*y), y.fun=threshold.type*y.fun.val, box=list(box.init), box.mass=mass.init, num.class=1, num.hdr.class=1, threshold=do.call(y.fun, list(x=y)))
    class(x.prim) <- "prim"
    
    return(x.prim)
  }
  else
  {
    if (verbose)
      cat(paste("Found box ", k, ": y.fun=", signif(threshold.type*boxk$y.fun,4), ", mass=", signif(boxk$mass,4), "\n\n", sep=""))
    boxes <- list(x=list(boxk$x), y=list(boxk$y), y.fun=list(boxk$y.fun),
                  box=list(boxk$box), mass=list(boxk$mass))       
  }
    
  ## find subsequent boxes
  if (num.boxes > 1)
  {
    boxk <- list(x=boxes$x[[k]], y=boxes$y[[k]], y.fun=boxes$y.fun[[k]],
                 box=boxes$box[[k]], mass=boxes$mass[[k]])

    ## data still under consideration
    x.out.ind.mat <-  matrix(TRUE, nrow=nrow(x), ncol=ncol(x))
    for (j in 1:d)
      x.out.ind.mat[,j] <- (x[,j] < boxk$box[1,j]) | (x[,j] > boxk$box[2,j])

    x.out.ind <- apply(x.out.ind.mat, 1, sum)!=0
    
    x.out <- x[x.out.ind,]
    if (is.vector(x.out)) x.out <- as.matrix(t(x.out)) 
    y.out <- y[x.out.ind]
      
    ##box.out <- apply(x.out, 2, range)
    while ((length(y.out)>0) & (k < num.boxes) & (!is.null(boxk))) 
    {
      k <- k+1
      
      boxk <- find.box(x=x.out, y=y.out, box=box.init,
                       peel.alpha=peel.alpha, paste.alpha=paste.alpha,
                       mass.min=mass.min, threshold=min(y)-0.1*abs(min(y)), d=d, n=n,
                       pasting=pasting, verbose=verbose, y.fun=y.fun)

      if (is.null(boxk))
      {
        if (verbose)
          cat(paste("Bump", k, "includes all remaining data\n\n"))

        boxes$x[[k]] <- x.out
        boxes$y[[k]] <- y.out
        boxes$y.fun[[k]] <- do.call(y.fun, list(x=y.out))
        boxes$box[[k]] <- box.init
        boxes$mass[[k]] <- length(y.out)/n
      }
      else 
      {
        ## update x and y
        if (verbose)
          cat(paste("Found box ", k, ": y.fun=", signif(threshold.type*boxk$y.fun,4),
                    ", mass=", signif(boxk$mass,4), "\n\n", sep=""))
        
        x.out.ind.mat <- matrix(TRUE, nrow=nrow(x), ncol=ncol(x))
        for (j in 1:d)
          x.out.ind.mat[,j] <- (x[,j] < boxk$box[1,j]) | (x[,j] > boxk$box[2,j])
        
        x.out.ind <- x.out.ind & (apply(x.out.ind.mat, 1, sum)!=0)
        x.out <- x[x.out.ind,]
        if (is.vector(x.out)) x.out <- as.matrix(t(x.out))
        y.out <- y[x.out.ind]
     
        boxes$x[[k]] <- boxk$x
        boxes$y[[k]] <- boxk$y
        boxes$y.fun[[k]] <- boxk$y.fun
        boxes$box[[k]] <- boxk$box
        boxes$mass[[k]] <-boxk$mass 
      }
    }   
  }

  ## adjust for negative hdr  
  for (k in 1:length(boxes$y.fun))
  {
    boxes$y[[k]] <- threshold.type*boxes$y[[k]]
    boxes$y.fun[[k]] <- threshold.type*boxes$y.fun[[k]]
  }
  
  ## highest density region

  prim.res <- prim.hdr(prim=boxes, threshold=threshold, threshold.type=threshold.type, y.fun=y.fun)
  
  return(prim.res)
         
}

###############################################################################
## Highest density region for PRIM boxes
###############################################################################

prim.hdr <- function(prim, threshold, threshold.type, y.fun=mean)
{  
  n <- 0
  for (i in 1:length(prim$box))
    n <- n + length(prim$y[[i]])

  hdr.ind <- which(unlist(prim$y.fun)*threshold.type >= threshold*threshold.type)

  if (length(hdr.ind) > 0)
    hdr.ind <- max(hdr.ind)
  else
  {
    if (threshold.type==1)
      warning(paste("No prim box found with y.fun >=", threshold))
    else if (threshold.type==-1)
      warning(paste("No prim box found with y.fun <=", threshold))
    return()
  }
  
  ## highest density region  
  x.prim.hdr <- list()
  
  for (k in 1:hdr.ind)
  {
    x.prim.hdr$x[[k]] <- prim$x[[k]]
    x.prim.hdr$y[[k]] <- prim$y[[k]]
    x.prim.hdr$y.fun[[k]] <- prim$y.fun[[k]]
    x.prim.hdr$box[[k]] <- prim$box[[k]]
    x.prim.hdr$mass[[k]] <-prim$mass[[k]] 
  }
  
  ## combine non-hdr into a `dump' box
  if (hdr.ind < length(prim$x))
  {
    x.temp <- numeric()
    y.temp <- numeric()
    for (k in (hdr.ind+1):length(prim$x))
    {
      x.temp <- rbind(x.temp, prim$x[[k]])
      y.temp <- c(y.temp, prim$y[[k]])
    }
    
    x.prim.hdr$x[[hdr.ind+1]] <- x.temp
    x.prim.hdr$y[[hdr.ind+1]] <- y.temp
    x.prim.hdr$y.fun[[hdr.ind+1]] <- do.call(y.fun, list(x=y.temp))
    x.prim.hdr$box[[hdr.ind+1]] <- prim$box[[length(prim$x)]]
    x.prim.hdr$mass[[hdr.ind+1]] <- length(y.temp)/n  
  }
  
  x.prim.hdr$num.class <- length(x.prim.hdr$x)
  x.prim.hdr$num.hdr.class <- hdr.ind
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

prim.combine <- function(prim1, prim2, y.fun=mean)
{
  M1 <- prim1$num.hdr.class
  M2 <- prim2$num.hdr.class

  if (is.null(M1) & !is.null(M2))
    return (prim2)
  if (!is.null(M1) & is.null(M2))
    return(prim1)
  if (is.null(M1) & is.null(M2))
    return(NULL)

  overlap <- overlap.box.seq(prim1, prim2, rel.tol=0.01)

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
      prim.temp$y.fun[[i]] <- prim1$y.fun[[i]]
      prim.temp$box[[i]] <- prim1$box[[i]]
      prim.temp$mass[[i]] <- prim1$mass[[i]]
      prim.temp$ind[[i]] <- 1
    }
    for (i in 1:M2)
    {
      prim.temp$x[[i+M1]] <- prim2$x[[i]]
      prim.temp$y[[i+M1]] <- prim2$y[[i]]
      prim.temp$y.fun[[i+M1]] <- prim2$y.fun[[i]]
      prim.temp$box[[i+M1]] <- prim2$box[[i]]
      prim.temp$mass[[i+M1]] <- prim2$mass[[i]]
      prim.temp$ind[[i+M1]] <- -1
    }
    
    dumpx.ind <- prim.which.box(x, prim1)==prim1$num.class & prim.which.box(x, prim2)==prim2$num.class
    
    prim.temp$x[[M1+M2+1]] <- x[dumpx.ind,]
    prim.temp$y[[M1+M2+1]] <- y[dumpx.ind]
    prim.temp$y.fun[[M1+M2+1]] <- do.call(y.fun, list(x=y[dumpx.ind]))
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
##         - FALSE - dont include pasting
##
## Returns
## List with fields
## x - data still inside box after peeling
## y - corresponding response values
## y.mean - mean of y
## box - box limits
## mass - box mass
###############################################################################

find.box <- function(x, y, box, peel.alpha, paste.alpha, mass.min, threshold, d, n, pasting, verbose, y.fun=mean) 
{
  y.fun.val <- do.call(y.fun, list(x=y))
  mass <- length(y)/n
  
  if ((y.fun.val >= threshold) & (mass >= mass.min))
    boxk.peel <- peel.one(x=x, y=y, box=box, peel.alpha=peel.alpha,
                          mass.min=mass.min,threshold=threshold, d=d, n=n, y.fun=y.fun, verbose=verbose)   
  else
    boxk.peel <- NULL

  boxk.temp <- NULL
 
  while (!is.null(boxk.peel))
  { 
    boxk.temp <- boxk.peel
    boxk.peel <- peel.one(x=boxk.temp$x, y=boxk.temp$y, box=boxk.temp$box,
                          peel.alpha=peel.alpha,
                          mass.min=mass.min, threshold=threshold, d=d, n=n, y.fun=y.fun, verbose=verbose)
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
                              mass.min=mass.min, threshold=threshold, d=d, n=n, y.fun=y.fun, verbose=verbose)      
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

peel.one <- function(x, y, box, peel.alpha, mass.min, threshold, d, n, type=8, y.fun=mean, verbose)
{
  box.new <- box
  mass <- length(y)/n

  if (is.vector(x)) return(NULL)
  
  ##y.mean <- mean(y)
  y.fun.val <- do.call(y.fun, list(x=y))
  
  y.fun.peel <- matrix(0, nrow=2, ncol=d)
  box.vol.peel <- matrix(0, nrow=2, ncol=d) 
                     
  for (j in 1:d)
  {
    box.min.new <- quantile(x[,j], peel.alpha, type=type)
    box.max.new <- quantile(x[,j], 1-peel.alpha, type=type)

    y.fun.peel[1,j] <- do.call(y.fun, list(x=y[x[,j] >= box.min.new]))
    y.fun.peel[2,j] <- do.call(y.fun, list(x=y[x[,j] <= box.max.new]))

    box.temp1 <- box
    box.temp2 <- box
    box.temp1[1,j] <- box.min.new
    box.temp2[2,j] <- box.max.new
    box.vol.peel[1,j] <- vol.box(box.temp1)
    box.vol.peel[2,j] <- vol.box(box.temp2)    
  }
  
  y.fun.peel.max.ind <- which(y.fun.peel==max(y.fun.peel, na.rm=TRUE), arr.ind=TRUE)

  ## break ties by choosing box with largest volume

  nrr <- nrow(y.fun.peel.max.ind) 
  if (nrr > 1)
  {
    box.vol.peel2 <- rep(0, nrr)
    for (j in 1:nrr)
      box.vol.peel2[j] <- box.vol.peel[y.fun.peel.max.ind[j,1],
                                       y.fun.peel.max.ind[j,2]]
    
    row.ind <- which(max(box.vol.peel2)==box.vol.peel2)
  }
  else
    row.ind <- 1
  
  y.fun.peel.max.ind <- y.fun.peel.max.ind[row.ind,]
  ## peel along dimension j.max
  j.max <- y.fun.peel.max.ind[2]

  ## peel lower 
  if (y.fun.peel.max.ind[1]==1)
  {
    box.new[1,j.max] <- quantile(x[,j.max], peel.alpha, type=type)
    x.index <- x[,j.max] >= box.new[1,j.max] 
  }
  ## peel upper 
  else if (y.fun.peel.max.ind[1]==2)
  {
    box.new[2,j.max] <- quantile(x[,j.max], 1-peel.alpha, type=type)
    x.index <- x[,j.max] <= box.new[2,j.max]
  }

  if (verbose) cat("Peeled in dimension", j.max, ": new limits are",  paste("[", signif(box.new[1,j.max],4), ",", signif(box.new[2,j.max],4) , "]\n", sep=""))
  x.new <- x[x.index,]
  y.new <- y[x.index]
  mass.new <- length(y.new)/n
  y.fun.new <- do.call(y.fun, list(x=y.new))

  ## if min. y.fun and min. mass conditions are still true, update
  ## o/w return NULL  

  if ((y.fun.new >= threshold) & (mass.new >= mass.min) & (mass.new < mass))
    return(list(x=x.new, y=y.new, y.fun=y.fun.new, box=box.new,
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

paste.one <- function(x, y, x.init, y.init, box, paste.alpha, mass.min, threshold, d, n, y.fun=mean, verbose)
{
  box.new <- box
  mass <- length(y)/n
  y.fun.val <- do.call(y.fun, list(x=y))
  n.box <- length(y)

  box.init <- apply(x.init, 2, range)
  
  if (is.vector(x)) x <- as.matrix(t(x))
  
  y.fun.paste <- matrix(0, nrow=2, ncol=d)
  mass.paste <- matrix(0, nrow=2, ncol=d)
  box.paste <- matrix(0, nrow=2, ncol=d)
  x.paste1.list <- list()
  x.paste2.list <- list()
  y.paste1.list <- list()
  y.paste2.list <- list()

  box.paste1 <- box
  box.paste2 <- box
   
  for (j in 1:d)
  {    
    ## candidates for pasting
    box.diff <- (box.init[2,] - box.init[1,])[j]

    box.paste1[1,j] <- box[1,j] - box.diff*paste.alpha
    box.paste2[2,j] <- box[2,j] + box.diff*paste.alpha
    
    x.paste1.ind <- in.box(x=x.init, box=box.paste1, d=d, boolean=TRUE)
    x.paste1 <- x.init[x.paste1.ind,]
    y.paste1 <- y.init[x.paste1.ind]
    
    x.paste2.ind <- in.box(x=x.init, box=box.paste2, d=d, boolean=TRUE)
    x.paste2 <- x.init[x.paste2.ind,]
    y.paste2 <- y.init[x.paste2.ind]
      
    while (length(y.paste1) <= length(y) & box.paste1[1,j] >= box.init[1,j])
    {
      box.paste1[1,j] <- box.paste1[1,j] - box.diff*paste.alpha
      x.paste1.ind <- in.box(x=x.init, box=box.paste1, d=d, boolean=TRUE)
      x.paste1 <- x.init[x.paste1.ind,]
      y.paste1 <- y.init[x.paste1.ind]
    }
    
    while (length(y.paste2) <= length(y) & box.paste2[2,j] <= box.init[2,j])
    {
      box.paste2[2,j] <- box.paste2[2,j] + box.diff*paste.alpha
      x.paste2.ind <- in.box(x=x.init, box=box.paste2, d=d, boolean=TRUE)
      x.paste2 <- x.init[x.paste2.ind,]
      y.paste2 <- y.init[x.paste2.ind]
    }

   
    ## y means of pasted boxes
    y.fun.paste[1,j] <- do.call(y.fun, list(x=y.paste1))
    y.fun.paste[2,j] <- do.call(y.fun, list(x=y.paste2))

    ## mass of pasted boxes
    mass.paste[1,j] <- length(y.paste1)/n
    mass.paste[2,j] <- length(y.paste2)/n
    
    x.paste1.list[[j]] <- x.paste1
    y.paste1.list[[j]] <- y.paste1
    x.paste2.list[[j]] <- x.paste2
    y.paste2.list[[j]] <- y.paste2
    box.paste[1,j] <- box.paste1[1,j]
    box.paste[2,j] <- box.paste2[2,j]
  }

  ## break ties by choosing box with largest mass
  
  y.fun.paste.max <- which(y.fun.paste==max(y.fun.paste, na.rm=TRUE), arr.ind=TRUE)
  
  if (nrow(y.fun.paste.max)>1)
  {
     y.fun.paste.max <- cbind(y.fun.paste.max, mass.paste[y.fun.paste.max])
     y.fun.paste.max.ind <- y.fun.paste.max[order(y.fun.paste.max[,3], decreasing=TRUE),][1,1:2]
  }
  else
    y.fun.paste.max.ind <- as.vector(y.fun.paste.max)       
  
  ## paste along dimension j.max
  j.max <- y.fun.paste.max.ind[2]

  ## paste lower 
  if (y.fun.paste.max.ind[1]==1)
  {
     x.new <- x.paste1.list[[j.max]] 
     y.new <- y.paste1.list[[j.max]]
     box.new[1,j.max] <- box.paste[1,j.max]
  }   
  ## paste upper
  else if (y.fun.paste.max.ind[1]==2)
  {
     x.new <- x.paste2.list[[j.max]] 
     y.new <- y.paste2.list[[j.max]]
     box.new[2,j.max] <- box.paste[2,j.max]
  }

  if (verbose) cat("Pasted in dimension", j.max, ": new limits are",  paste("[", signif(box.new[1,j.max],4), ",", signif(box.new[2,j.max],4) , "]\n", sep=""))
  mass.new <- length(y.new)/n
  y.fun.new <- do.call(y.fun, list(x=y.new))

  if ((y.fun.new > threshold) & (mass.new >= mass.min) & (y.fun.new >= y.fun.val)
      & (mass.new > mass))
    return(list(x=x.new, y=y.new, y.fun=y.fun.new, box=box.new, mass=mass.new))
 
}




###############################################################################
## Output functions for prim objects
###############################################################################

###############################################################################
## Plot function for PRIM objects
###############################################################################


plot.prim <- function(x, splom=TRUE, ...)
{
  if (ncol(x$x[[1]])==2)
    plotprim.2d(x, ...)
  else if (ncol(x$x[[1]])==3 & !splom)
  {  
     plotprim.3d(x, ...)
  }
  else if (ncol(x$x[[1]])>3 | (ncol(x$x[[1]])==3 & splom))
    plotprim.nd(x, ...)
  
  invisible()
}

plotprim.2d <- function(x, col, xlim, ylim, xlab, ylab, add=FALSE,  add.legend=FALSE, cex.legend=1, pos.legend, lwd=1, border, col.vec=c("blue", "orange"), alpha=1, ...)
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
  {  
    col <- rep("transparent", M)
    col[which(x$ind==1)] <- scales::alpha(col.vec[2], alpha=alpha)
    col[which(x$ind==-1)] <- scales::alpha(col.vec[1], alpha=alpha)
  }

  if (missing(border)) border <- rep(1, M)
  if (length(col) < M)
    col <- rep(col, length=M)
    
  for (i in M:1)
  {
    ## colour i-th box
    box <- x$box[[i]]
    rect(box[1,1], box[1,2], box[2,1], box[2,2], border=border[i], col=col[i], lwd=lwd)
  }
  
  if (add.legend)
    legend(pos.legend[1], pos.legend[2], legend=text.legend, fill=col, bty="n", cex=cex.legend)
  
  invisible()
}


plotprim.3d <- function(x, col, xlim, ylim, zlim, xlab, ylab, zlab, col.vec=c("blue","orange"), alpha=1, theta=30, phi=40, d=4, ...)
{
  M <- x$num.hdr.class
  if (missing(col))
  {
    col <- rep("transparent", M)
    col[which(x$ind==1)] <- scales::alpha(col.vec[2], alpha=alpha)
    col[which(x$ind==-1)] <- scales::alpha(col.vec[1], alpha=alpha)
  }
  
  if (missing(xlab)) xlab <- names(x$x[[1]])[1]
  if (missing(ylab)) ylab <- names(x$x[[1]])[2]
  if (missing(zlab)) zlab <- names(x$x[[1]])[3]
  
  xdata <- numeric() 
  for (i in M:1) xdata <- rbind(xdata, x$x[[i]])
  
  xprim <- prim.which.box(xdata, box.seq=x)
  xprim.ord <- order(xprim)
  xdata <- xdata[xprim.ord,]
  xprim.col <- col[xprim][xprim.ord]
  
  plot3D::points3D(xdata[,1], xdata[,2], xdata[,3], col=xprim.col, theta=theta, phi=phi, d=d, xlab=xlab, ylab=ylab, zlab=zlab, ...)
  
  invisible()
}



plotprim.nd  <- function(x, col, xmin, xmax, xlab, ylab, x.pt, m, col.vec=c("blue","orange"), alpha=1, ...)
{
  M <- x$num.hdr.class
  d <- ncol(x$x)
  if (missing(col))
  {
    col <- rep("transparent", M)
    col[which(x$ind==1)] <- scales::alpha(col.vec[2], alpha=alpha)
    col[which(x$ind==-1)] <- scales::alpha(col.vec[1], alpha=alpha)
  }
  if (missing(m) & !missing(x.pt)) m <- round(nrow(x.pt)/10)
  if (missing(x.pt))
  {
    x.pt <- numeric()
    for (j in 1:length(x$x))
      x.pt <- rbind(x.pt,x$x[[j]])
    if (missing(m)) m <- max(round(nrow(x.pt)/10), nrow(x$x[[1]]))
    
    x.pt <- x.pt[sample(1:nrow(x.pt), size=m),]
  }

  xprim <- prim.which.box(x.pt, box.seq=x)
  xprim.ord <- order(xprim)
  x.pt <- x.pt[xprim.ord,]
  xprim.col <- col[xprim][xprim.ord]
 
  pairs(x.pt, col=xprim.col,  ...)

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

summary.prim <- function(object, ..., print.box=FALSE)
{
  x <- object
  M <- x$num.class

  if (M>1)
  {
      x2 <- lapply(lapply(x[c("y.fun","mass", "ind")], head, n=M-1), unlist)
      summ.mat <- data.frame(y.fun=x2$y.fun, mass=x2$mass, ind=x2$ind) 
      summ.mat <- rbind(summ.mat, c(x$y.fun[[M]], x$mass[[M]], NA))
      
      tot <- c(sum(summ.mat[,1]*summ.mat[,2])/sum(summ.mat[,2]), sum(summ.mat[,2]), NA)
      summ.mat <- rbind(summ.mat, tot)
      
    rownames(summ.mat) <- c(paste("box", 1:(nrow(summ.mat)-1), sep=""), "overall")
      colnames(summ.mat) <- c("box-fun", "box-mass", "threshold.type")
      
      if (x$num.hdr.class < x$num.class)
          for (k in (x$num.hdr.class+1):x$num.class)
              rownames(summ.mat)[k] <- paste(rownames(summ.mat)[k], "*",sep="")
  }
  else
  {
      summ.mat <- c(x$y.fun[[1]], x$mass[[1]])
      tot <- c(x$y.fun[[1]], x$mass[[1]])
      summ.mat <- rbind(summ.mat, tot)
      rownames(summ.mat) <- c(paste("box", 1:(nrow(summ.mat)-1), sep=""), "overall")
      colnames(summ.mat) <- c("box-fun", "box-mass")
  }
  
  print(summ.mat)
  cat("\n")
  
  if (print.box)
  { 
    for (k in 1:M)
    {
      cat(paste("Box limits for box", k, "\n", sep=""))
       box.summ <- x$box[[k]]
       rownames(box.summ) <- c("min", "max")
       print(box.summ)
      cat("\n")
    }
  }
}


predict.prim <- function(object, newdata, y.fun.flag=FALSE, ...)
{
  x <- newdata
  prim.obj <- object 
  which.box <- prim.which.box(x=x, box.seq=prim.obj)

  if (y.fun.flag)
  {
    lab <- rep(tail(prim.obj$y.fun, n=1), length(which.box))
    lab[which.box>0] <- prim.obj$y.fun[which.box[which.box>0]]
  }
  else
  {
    lab <- which.box
  }

  return(lab)
}

