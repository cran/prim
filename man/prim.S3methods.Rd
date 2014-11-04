\name{prim S3 methods}
\alias{predict.prim}
\alias{summary.prim}

\title{S3 methods for PRIM for multivariate data}
\description{
  S3 methods PRIM for multivariate data.
}
\usage{
\method{predict}{prim}(object, newdata, y.fun.flag=FALSE, ...)
\method{summary}{prim}(object, ..., print.box=FALSE)
}

\arguments{
  \item{object}{object of class \code{prim}}
  \item{newdata}{data matrix}
  \item{y.fun.flag}{flag to return y value of PRIM box rather than 
    box label. Default is FALSE.}
  \item{print.box}{flag to print out limits of all PRIM boxes. Default
  is FALSE.} 
  \item{...}{other parameters}
}
  


\details{
  --The \code{predict} method returns the value of PRIM box number in
  which \code{newdata} are located.

  --The \code{summary} method displays a table with three columns:
    \code{box-fun} is the y value, \code{box-mass} is the mass of the
    box,  \code{threshold.type} is the threshold direction indicator: 1
    = ">= threshold", -1 = "<=threshold". Each box corresponds to a
    row. The second last row marked with an asterisk is the box
    which collates the remaining data points not belonging to a specific
    PRIM box. The final row is an overall summary, i.e. \code{box-fun} is the
    overall mean of y and \code{box-mass} is 1.  
}


\examples{
data(quasiflow)
qf <- quasiflow[1:1000,1:2]
qf.label <- quasiflow[1:1000,3]*quasiflow[1:1000,4]

qf.prim <- prim.box(x=qf, y=qf.label, threshold=c(0.3, -0.1), threshold.type=0,
  verbose=TRUE)
## verbose=TRUE prints out extra informaton about peeling and pasting

summary(qf.prim)
predict(qf.prim, newdata=c(0.6,0.2))

## using median insted of mean for the response y

qf.prim2 <- prim.box(x=qf, y=qf.label, threshold=c(0.5, -0.2),
  threshold.type=0, y.fun=median)
summary(qf.prim2)
predict(qf.prim2, newdata=c(0.6,0.2))
}

\keyword{multivariate}
