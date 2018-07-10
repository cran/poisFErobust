\name{poisFErobust-package}
\alias{poisFErobust-package}
\alias{poisFErobust}
\docType{package}
\title{
\packageTitle{poisFErobust}
}
\description{
\packageDescription{poisFErobust}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{poisFErobust}
\packageIndices{poisFErobust}
}
\author{
\packageAuthor{poisFErobust}

Maintainer: \packageMaintainer{poisFErobust}
}
\references{
Wooldridge, Jeffrey M. (1999): "Distribution-free estimation of some nonlinear
panel data models," Journal of Econometrics, 90, 77-97.
}

\keyword{ package }
\keyword{ regression }

\examples{
# ex.dt.good satisfies the conditional mean assumption
data("ex.dt.good")
pois.fe.robust(outcome = "y", xvars = c("x1", "x2"), group.name = "id",
               index.name = "id", data = ex.dt.good)
               
# ex.dt.bad violates the conditional mean assumption
data("ex.dt.bad")
pois.fe.robust(outcome = "y", xvars = c("x1", "x2"), group.name = "id",
               index.name = "id", data = ex.dt.bad)
}
