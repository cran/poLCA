\name{poLCA-internal}
\alias{poLCA.updatePrior}
\alias{poLCA.postClass.C}
\alias{poLCA.probHat.C}
\alias{poLCA.ylik.C}
\alias{poLCA.dLL2dBeta.C}
\alias{poLCA.vectorize}
\alias{poLCA.unvectorize}
\alias{poLCA.compress}
\title{Internal functions for poLCA}
\description{Internal poLCA functions}
\usage{
poLCA.updatePrior(b, x, R)
poLCA.postClass.C(prior, vp, y)
poLCA.probHat.C(rgivy, y, vp)
poLCA.ylik.C(vp, y)
poLCA.dLL2dBeta.C(rgivy, prior, x)
poLCA.vectorize(probs)
poLCA.unvectorize(vp)
poLCA.compress(y)
}
\arguments{
\item{b}{vector of multinomial logit coefficients on covariates.}
\item{R}{number of classes assumed in fitting model.}
\item{prior}{matrix of prior probabilities.}
\item{probs}{a list of length equal to the number of outcome variables, containing estimated conditional item response probabilities for each question and class.}
\item{vp}{vectorized format of \code{probs} for use in C code.}
\item{rgivy}{matrix of posterior probabilities of observation being in class \emph{r} given data \emph{y}.}
\item{x}{matrix of covariates.}
\item{y}{a matrix or data frame of manifest variables.}
}
\value{
\item{poLCA.updatePrior}{When estimating latent class regression models (with covariates), function uses multinomial logit link to update prior probabilities of class membership based upon estimated coefficients and values of the concomitant variables in EM algorithm. Returns matrix of prior probabilities, with one column per class being estimated, and each row representing one observation.}
\item{poLCA.postClass.C}{C call. Calculates posterior probabilities of class membership as prior times likelihood, in EM algorithm.  Returns matrix of posterior probabilities of class memberships, with one column per class being estimated, and each row representing one observation, so that every row sums to 1.}
\item{poLCA.probHat.C}{C call. Calculates updated estimates of manifest variable class-conditional reponse probabilities, in EM algorithm.  Returns list of length equal to number of manifest variables, with each item containing a matrix of conditional response probabilties.}
\item{poLCA.ylik.C}{C call. Calculates for each observation the probabilities of observing the full set of responses to all manifest variables, given conditional response probabilities \code{probs}, for each class.  Probabilities are calculated assuming local independence, meaning that for each observation and each class, the item response probabilities are simply multiplied across the manifest variables.  Returns a matrix with rows equal to number of observations and columns equal to number of classes specified by \code{nclass} in \code{\link{poLCA}}}
\item{poLCA.dLL2dBeta.C}{C call. Calculates gradient vector and Hessian matrix when estimating model with covariates, for updating multinomial logit coefficients \code{b} in EM algorithm.}
\item{poLCA.vectorize}{Converts conditional response probabilties from list format \code{probs} to vector format \code{vp} for use in compiled C code.}
\item{poLCA.unvectorize}{Converts conditional response probabilties from vector format \code{vp} to list format \code{probs} for use in R.}
\item{poLCA.compress}{Compress a matrix or data frame by removing redundant rows, so that each unique row appears only once.  Also counts and returns the frequency of each row pattern.  The resulting matrix is sorted in ascending order, by column, starting with the first column.  Returns a list containing the "compressed" matrix \code{datamat}, in the same format (data frame or matrix) as what was inputted, along with the frequency count \code{freq} of each unique row in that matrix.}
}
\details{
Functions used in updating EM algorithm in main poLCA function, including calls to compiled C code to accelerate estimation.  
  
These functions are masked from the user and are only meant to be called internally from \code{\link{poLCA}}.
}
\seealso{\code{\link{poLCA}}}
\keyword{methods}
