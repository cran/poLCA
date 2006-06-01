\name{poLCA.makeplot.dich}
\alias{poLCA.makeplot.dich}
\title{Plot multi-item dichotomous response probabilities}
\description{This function plots estimated response probabilities to a series of dichotomous outcome variables, conditional on estimated class membership, for an aribitrary number of classes.  It also shows the estimated population proportion in each class.}
\usage{poLCA.makeplot.dich(probs, P, y, ti)}
\arguments{
\item{probs}{a list of length equal to the number of outcome variables, containing estimated conditional item response probabilities for each question and class.}
\item{P}{class membership proportions.}
\item{y}{matrix of observed outcome variables.}
\item{ti}{title of upper plot.}
}
\details{This function is designed to be called internally from \code{\link{poLCA}} only when the attribute \code{graphs} is set to \code{TRUE}, \emph{and} only when every outcome item is dichotomous.  If one or more dependent variables are polytomous, then the \code{\link{poLCA}} function will instead call \code{\link{poLCA.makeplot.poly}} to plot the relevant information.}
\note{This function is simply an elaborate invocation of the \code{\link[scatterplot3d:scatterplot3d]{scatterplot3d}} command in the \pkg{scatterplot3d} package.}
\seealso{\code{\link{poLCA}}, \code{\link{poLCA.makeplot.poly}}}
\keyword{hplot}
