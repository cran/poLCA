\name{poLCA.makeplot.poly}
\alias{poLCA.makeplot.poly}
\title{Plot multi-item polytomous response probabilities}
\description{This function plots estimated response probabilities to a series of polytomous outcome variables, conditional on estimated class membership.}
\usage{poLCA.makeplot.poly(probs, r, y, K.j, ti)}
\arguments{
\item{probs}{a list of length equal to the number of outcome variables, containing estimated conditional item response probabilities for each question and class.}
\item{r}{class for which to plot response probabilities.}
\item{y}{matrix of observed outcome variables.}
\item{K.j}{vector containing number of response items for each dependent variable.}
\item{ti}{title of probability plots.}
}
\details{This function is designed to be called internally from \code{\link{poLCA}} only when the attribute \code{graphs} is set to \code{TRUE}, \emph{and} only if at least one outcome item is polytomous.  If all dependent variables are dichotomous, then the \code{\link{poLCA}} function will instead call \code{\link{poLCA.makeplot.dich}} to plot the relevant information.}
\note{This function is simply an elaborate invocation of the \code{\link[scatterplot3d:scatterplot3d]{scatterplot3d}} command in the \pkg{scatterplot3d} package.}
\seealso{\code{\link{poLCA}}, \code{\link{poLCA.makeplot.dich}}}
\keyword{hplot}

