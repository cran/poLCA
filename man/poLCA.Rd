\name{poLCA}
\alias{poLCA}
\title{Latent class analysis of polytomous outcome variables}
\description{Estimates latent class and latent class regression models for polytomous outcome variables.}
\usage{
poLCA(formula, data, nclass = 2, maxiter = 1000, graphs = FALSE, 
      tol = 1e-10, na.rm = TRUE)
}
\arguments{
\item{formula}{a formula expression of the form \code{response ~ predictors}. The details of model specification are given below.}
\item{data}{a data frame containing variables in \code{formula}. Manifest variables must contain \emph{only} integer values, and must be coded with consecutive values from 1 to the maximum number of outcomes for each variable. All missing values should be entered as \code{NA}.}
\item{nclass}{The number of latent classes to assume in the model. Setting \code{nclass=1} results in \code{poLCA} estimating the loglinear independence model. The default is two.}
\item{maxiter}{The maximum number of iterations through which the estimation algorithm will cycle.}
\item{graphs}{Logical, for whether \code{poLCA} should graphically display the parameter estimates at each stage of the updating algorithm. The default is \code{FALSE}.}
\item{tol}{A tolerance value for judging when convergence has been reached.  When the one-iteration change in the estimated log-likelihood is less than \code{tol}, the estimation algorithm stops updating and considers the maximum log-likelihood to have been found.}
\item{na.rm}{Logical, for how \code{poLCA} handles cases with missing values on the manifest variables.  If \code{TRUE}, those cases are removed (listwise deleted) before estimating the model. If \code{FALSE}, cases with missing values are retained. Cases with missing covariates are always removed. The default is \code{TRUE}.}
}
\details{
Latent class analysis, also known as latent structure analysis, is a technique for the analysis of clustering among observations in multi-way tables of qualitative/categorical variables.  The central idea is to fit a model in which any confounding between the manifest variables can be explained by a single unobserved "latent" categorical variable.  \code{poLCA} uses the assumption of local independence to estimate a mixture model of latent multi-way tables, the number of which (\code{nclass}) is specified by the user.  Estimated parameters include the class-conditional response probabilities for each manifest variable, the "mixing" proportions denoting population share of observations corresponding to each latent multi-way table, and coefficients on any class-predictor covariates, if specified in the model.
  
Model specification: Latent class models have more than one manifest variable, so the response variables are \code{cbind(dv1,dv2,dv3...)} where \code{dv#} refer to variable names in the data frame.  For models with no covariates, the formula is \code{cbind(dv1,dv2,dv3)~1}.  For models with covariates, replace the \code{~1} with the desired function of predictors \code{iv1,iv2,iv3...} as, for example, \code{cbind(dv1,dv2,dv3)~iv1+iv2*iv3}.
  
\code{poLCA} treats all manifest variables as qualitative/categorical/nominal -- NOT as ordinal.
}
\value{
\code{poLCA} returns a list containing the following components:
\item{y}{data frame of manifest variables.}
\item{x}{data frame of covariates, if specified.}
\item{N}{number of cases used in model.}
\item{Nobs}{number of fully observed cases (less than or equal to \code{N}).}
\item{probs}{estimated class-conditional response probabilities.}
\item{P}{sizes of each latent class; equal to the mixing proportions in the basic latent class model, or the mean of the posteriors in the latent class regression model.}
\item{posterior}{matrix of posterior class membership probabilities.}
\item{predclass}{vector of predicted class memberships, by modal assignment.}
\item{predcell}{table of observed vs. predicted cell counts.}
\item{llik}{maximum value of the log-likelihood.}
\item{numiter}{number of iterations until reaching convergence.}
\item{coeff}{multinomial logit coefficient estimates on covariates (when estimated). \code{coeff} is a matrix with \code{nclass-1} columns, and one row for each covariate.  All logit coefficients are calculated for classes with respect to class 1.}
\item{coeff.se}{standard errors of coefficient estimates on covariates (when estimated), in the same format as \code{coeff}.}
\item{aic}{Akaike Information Criterion.}
\item{bic}{Bayesian Information Criterion.}
\item{Gsq}{Likelihood ratio/deviance statistic.}
\item{Chisq}{Pearson Chi-square goodness of fit statistic for fitted vs. observed multiway tables.}
\item{time}{length of time it took to run the model.}
\item{npar}{number of degrees of freedom used by the model (estimated parameters).}
\item{resid.df}{number of residual degrees of freedom.}
\item{eflag}{logical: error flag. True if estimation algorithm needed to automatically restart with new initial parameters.}
}
\references{
Agresti, Alan. 2002. \emph{Categorical Data Analysis, second edition}. Hoboken: John Wiley \& Sons.

Bandeen-Roche, Karen, Diana L. Miglioretti, Scott L. Zeger, and Paul J. Rathouz. 1997. "Latent Variable Regression for Multiple Discrete Outcomes." \emph{Journal of the American Statistical Association}. 92(440): 1375-1386.
 
Hagenaars, Jacques A. and Allan L. McCutcheon, eds. 2002. \emph{Applied Latent Class Analysis}. Cambridge: Cambridge University Press.
 
McLachlan, Geoffrey J. and Thriyambakam Krishnan. 1997. \emph{The EM Algorithm and Extensions}. New York: John Wiley \& Sons.
}
\note{
\code{poLCA} uses EM and Newton-Raphson algorithms to maximize the latent class model log-likelihood function. Depending on the starting parameters, this algorithm may only locate a local, rather than global, maximum. This becomes more and more of a problem as \code{nclass} increases. It is therefore highly advisable to run \code{poLCA} multiple times until you are certain you have located the global maximum log-likelihood. Each function call will use different random initial starting parameters.

The term "Latent class regression" (LCR) can have two meanings.  In this package, LCR models refer to latent class models in which the probability of class membership is predicted by one or more covariates.  However, in other contexts, LCR is also used to refer to regression models in which the manifest variable is partitioned into some specified number of latent classes as part of estimating the regression model. It is a way to simultaneously fit more than one regression to the data when the latent data partition is unknown. The \code{\link[fpc:regmix]{regmix}} function in package \pkg{fpc} will estimate this other type of LCR model.  Because of these terminology issues, the LCR models this package estimates are sometimes termed "latent class models with covariates" or "concomitant-variable latent class analysis," both of which are accurate descriptions of this model.
 
Latent class models are closely related to loglinear models; see package \pkg{gllm}.
}
\seealso{\code{\link{poLCA.simdata}}, \code{\link[e1071:lca]{lca}}, \code{\link[fpc:regmix]{regmix}}, \code{\link[gllm:gllm]{gllm}}}
\examples{
##
## Three models without covariates:
## M0: Loglinear independence model.
## M1: Two-class latent class model.
## M2: Three-class latent class model.
##
data(values)
f <- cbind(A,B,C,D)~1
M0 <- poLCA(f,values,nclass=1)			# log-likelihood: -543.6498
M1 <- poLCA(f,values,nclass=2)			# log-likelihood: -504.4677
M2 <- poLCA(f,values,nclass=3,maxiter=8000)	# log-likelihood: -503.3011

##
## Three-class model with a single covariate.
##
data(election)
f2a <- cbind(MORALG,CARESG,KNOWG,LEADG,DISHONG,INTELG,MORALB,CARESB,KNOWB,LEADB,DISHONB,INTELB)~PARTY
nes2a <- poLCA(f2a,election,nclass=3)		# log-likelihood: -16222.32 
pidmat <- cbind(1,c(1:7))
exb <- exp(pidmat \%*\% nes2a$coeff)
matplot(c(1:7),(cbind(1,exb)/(1+rowSums(exb))),ylim=c(0,1),type="l",
	main="Party ID as a predictor of candidate affinity class",
	xlab="Party ID: strong Democratic (1) to strong Republican (7)",
	ylab="Probability of latent class membership")
text(5.9,0.35,"Other")
text(5.4,0.7,"Bush affinity")
text(1.8,0.6,"Gore affinity")

##
## Create a sample dataset with missing values and estimate the
## latent class model including and excluding the missing values.
## Then plot the estimated class-conditional outcome response 
## probabilities against each other for the two methods.
##
simdat3 <- poLCA.simdata(N=5000,niv=1,ndv=5,nclass=2,b=c(-1,2),missval=TRUE,pctmiss=0.2)
f3 <- cbind(Y1,Y2,Y3,Y4,Y5)~X1
lc3.miss <- poLCA(f3,simdat3$dat,nclass=2)
lc3.nomiss <- poLCA(f3,simdat3$dat,nclass=2,na.rm=FALSE)
ifelse((order(lc3.miss$P)==order(lc3.nomiss$P)),o <- c(1,2),o <- c(2,1))
plot(lc3.miss$probs[[1]],lc3.nomiss$probs[[1]][o,],xlim=c(0,1),ylim=c(0,1),
     xlab="Conditional response probabilities (missing values dropped)",
     ylab="Conditional response probabilities (missing values included)")
for (i in 2:5) { points(lc3.miss$probs[[i]],lc3.nomiss$probs[[i]][o,]) }
abline(0,1,lty=3)
}
\keyword{methods}
