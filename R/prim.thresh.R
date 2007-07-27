
############################################################################
## PRIM threshold tuning parameter
## for difference between two normal mixtures
############################################################################

prim.thresh.mixt.dd <- function(prim.pos, prim.neg, mus1, Sigmas1, props1, mus2, Sigmas2, props2, weight.dd=1/2, weight.mixt=1/2, nmc=1e6, tau, xmc.dd, taumc.dd, gmc.mixt.dd, xmc.mixt, threshold.pos.range, threshold.neg.range, verbose=FALSE)
{

  if (missing(taumc.dd))
  {
    xmc.dd <- rmvnorm.mixt.dd(n=nmc, mus1, Sigmas1, props1, mus2, Sigmas2, props2) 
    taumc.dd <- thresh.mixt.dd(mus1=mus1, Sigmas1=Sigmas1, props1=props1, mus2=mus2, Sigmas2=Sigmas2, props2=props2, weight.dd=weight.dd, xmc=xmc.dd, tau=tau, unit.trans=FALSE)
  }
  
  if (missing(xmc.mixt))
    xmc.mixt <- rmvnorm.mixt.mixt(n=nmc, mus1=mus1, Sigmas1=Sigmas1, props1=props1, mus2=mus2, Sigmas2=Sigmas2, props2=props2, weight.mixt=weight.mixt)
  
  if (missing(gmc.mixt.dd))
    gmc.mixt.dd <- dmvnorm.mixt.dd(x=xmc.mixt, mus1=mus1, Sigmas1=Sigmas1, props1=props1, mus2=mus2, Sigmas2=Sigmas2, props2=props2, weight.dd=weight.dd)

  xinR.pos <- gmc.mixt.dd >= taumc.dd[1]
  xinR.neg <- gmc.mixt.dd <= taumc.dd[2]

  ## default threshold ranges
  if (missing(threshold.pos.range))
    threshold.pos.range <- sort(prim.pos$y.mean[prim.pos$y.mean>=mean(prim.pos$y.mean)])
   
  if (missing(threshold.neg.range))
    threshold.neg.range <- sort(prim.neg$y.mean[prim.neg$y.mean<=mean(prim.neg$y.mean)])

  if (verbose) cat("\n    Threshold   Error\n")
  ## step through threshold ranges to find min. error
  err.pos <- 0
  if (!missing(prim.pos))
  {  
    err.pos <- vector()
    i <- 0
    for (thp in threshold.pos.range)
    {
      prim.pos.hdr <- prim.hdr(prim=prim.pos, threshold=thp, threshold.type=1)
      
      xmc.mixt.class <- which.box(x=xmc.mixt, box.seq=prim.pos.hdr)
      xmc.mixt.class[xmc.mixt.class<=length(prim.pos.hdr$ind)] <- 1
      xmc.mixt.class[xmc.mixt.class>length(prim.pos.hdr$ind)] <- 0
      
      xinRhat <- xmc.mixt.class
      xinRhat.pos <- xinRhat==1
      symdiff.pos <- (xinR.pos & !xinRhat.pos) | (!xinR.pos & xinRhat.pos)
      i <- i+1
      err.pos[i] <- mean(symdiff.pos)
      if (verbose) print(c(thp, err.pos[i]))  
    }
    
    thresh.pos <- rev(threshold.pos.range)[which.min(rev(err.pos))]
  }
    
 
  err.neg <- 0
  if (!missing(prim.neg))
  {
    err.neg <- vector()
    i <- 0
    
    for (thm in threshold.neg.range)
    {
      prim.neg.hdr <- prim.hdr(prim=prim.neg, threshold=thm, threshold.type=-1)
      
      if (!is.null(prim.neg.hdr))
      {
        xmc.mixt.class <- which.box(x=xmc.mixt, box.seq=prim.neg.hdr)
        xmc.mixt.class[xmc.mixt.class<=length(prim.neg.hdr$ind)] <- -1
        xmc.mixt.class[xmc.mixt.class>length(prim.neg.hdr$ind)] <- 0
        
        xinRhat <- xmc.mixt.class
        xinRhat.neg <- xinRhat==-1
        symdiff.neg <- (xinR.neg & !xinRhat.neg) | (!xinR.neg & xinRhat.neg)
        i <- i+1
        err.neg[i] <- mean(symdiff.neg)}
      else
      {
        i <- i+1
        err.neg[i] <- 1
      }
      if (verbose) print(c(thm, err.neg[i]))
    }
  
    thresh.neg <- threshold.neg.range[which.min(err.neg)]
  }

  err <- min(err.pos) + min(err.neg)

  if (missing(prim.neg) & !missing(prim.pos))
    thresh.val <- thresh.pos
  if (missing(prim.pos) & !missing(prim.neg))
    thresh.val <- thresh.neg
  if (!missing(prim.pos) & !missing(prim.neg))
    thresh.val <- c(thresh.pos, thresh.neg) 
  if (missing(prim.pos) & missing(prim.neg))
    thresh.val <- NULL

  names(thresh.val) <- c("pos", "neg")
  return(thresh.val)
}

