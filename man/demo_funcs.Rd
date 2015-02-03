\name{demo_funcs}
\alias{demo.get.f}
\alias{demo.get.grad}
\alias{demo.get.hess.struct}
\title{Functions for the binary choice example in the vignette.}
\description{
Functions for the binary choice example in the vignette.
}
\usage{
demo.get.f(pars, Y, X, inv.Omega, inv.Sigma,T)
demo.get.grad(pars, Y, X, inv.Omega, inv.Sigma,T)
demo.get.hess.struct(N,k)
}

\arguments{
	\item{pars}{parameter vector}
	\item{Y}{Observed counts}
	\item{X}{Covariate matrix}
	\item{inv.Omega}{Hyperparameter}
	\item{inv.Sigma}{Hyperparameter}
	\item{N}{Number of agents (households)}
	\item{k}{Number of covariates}
	\item{T}{Number of purchase opportunities per household}
}
\value{
See vignette for details.
}






 


