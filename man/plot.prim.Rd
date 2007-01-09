\name{plot.prim}
\alias{plot.prim}

\title{PRIM plot for multivariate data}
\description{
  PRIM plot for multivariate data.
}

\synopsis{\method{plot}{prim}(x, ...)}

\usage{
## bivariate
\method{plot}{prim}(x, col, xlim, ylim, xlab, ylab, add=FALSE,
   add.legend=FALSE, cex.legend=1, pos.legend, lwd=1, ...)

## trivariate
\method{plot}{prim}(x, color, xlim, ylim, zlim, xlab, ylab, zlab,
   add.axis=TRUE, ...)

## d-variate
\method{plot}{prim}(x, col, xmin, xmax, xlab, ylab, ...)
}

\arguments{
  \item{x}{an object of class \code{prim}}
  \item{add.legend}{flag for adding legend (2-d plot)}
  \item{pos.legend}{(x,y) co-ordinates for leaned (2-d plot)}
  \item{cex.legend}{cex graphics parameter for legend (2-d plot)}
  \item{col}{vector of plotting colours, one for each box}
  \item{xlab,ylab,zlab,xlim,ylim,zlim,add,lwd}{usual graphics parameters}
  \item{xmin,xmax}{vector of minimum and maximum axis plotting values for scatter
    plot matrix}
  \item{color}{vector of colours, one for each box (3-d plot)}
  \item{add.axis}{flag for plotting axes (3-d plot)}
  \item{...}{other graphics parameters}
}
  

\value{
  Plot of 2-dim PRIM is a set of nested rectangles. Plot of 3-dim PRIM
  is a scatter point cloud. Plot of d-dim PRIM is a scatter plot
  matrix. The scatter plots indicate which points belong to which box.
}

\references{
  Friedman, J.H. \& Fisher, N.I. (1999) Bump-hunting for high
  dimensional data, \emph{Statistics and Computing}, \bold{9}, 123--143. 
}

\details{
  Default colours are \code{topo.colors()}, with one colour per box in the
  PRIM box sequence.   
}
 

\seealso{\code{\link{prim.box}}}

\examples{
## see ?prim.box
}

\keyword{ hplot}
