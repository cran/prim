\name{plot.prim}
\alias{plot.prim}

\title{PRIM plot for multivariate data}
\description{
  PRIM plot for multivariate data.
}

\usage{
\method{plot}{prim}(x, splom=TRUE,  ...)
}

\arguments{
  \item{x}{object of class \code{prim}}
  \item{splom}{flag for plotting 3-d data as scatter plot matrix. Default is TRUE.}
  \item{...}{other graphics parameters}
}
  

\value{
  Plot of 2-dim PRIM is a set of nested rectangles. Plot of 3-dim PRIM
  is a scatter point cloud. Plot of d-dim PRIM is a scatter plot
  matrix. The scatter plots indicate which points belong to which box.
}

\details{
The function headers are
\preformatted{
  ## bivariate
  plot(x, col, xlim, ylim, xlab, ylab, add=FALSE, add.legend=FALSE, 
    cex.legend=1, pos.legend, lwd=1, ...)

  ## trivariate
  plot(x, color, xlim, ylim, zlim, xlab, ylab, zlab, add.axis=TRUE, ...)

  ## d-variate
  plot(x, col, xmin, xmax, xlab, ylab, ...)
}
 The arguments are
 \describe{
  \item{\code{add.legend}}{flag for adding legend (2-d plot)}
  \item{\code{pos.legend}}{(x,y) co-ordinates for legend (2-d plot)}
  \item{\code{cex.legend}}{cex graphics parameter for legend (2-d plot)}
  \item{\code{col}}{vector of plotting colours, one for each box}
  \item{\code{xlab,ylab,zlab,xlim,ylim,zlim,add,lwd}}{usual graphics parameters}
  \item{\code{xmin,xmax}}{vector of minimum and maximum axis plotting values for scatter
    plot matrix}
  \item{\code{color}}{vector of colours, one for each box (3-d plot)}
  \item{\code{add.axis}}{flag for plotting axes (3-d plot)}
} 

  Default colours are \code{topo.colors()}, with one colour per box in the
  PRIM box sequence.   
}
 

\seealso{\code{\link{prim.box}}, \code{\link{predict.prim}}}

\examples{
## see ?predict.prim for bivariate example
## trivariate example
data(quasiflow)
qf <- quasiflow[1:1000,1:3]
qf.label <- quasiflow[1:1000,4]
thr <- c(0.25, -0.3)
qf.prim <- prim.box(x=qf, y=qf.label, threshold=thr, threshold.type=0)
plot(qf.prim)
}

\keyword{ hplot}
