\name{sample.GDS}
\alias{sample.GDS}
\title{Collect draws from the target posterior distribution.}
\description{
Runs the accept-reject phase of the Braun and Damien (2013) Generalized Direct Sampling algorithm.
	}
\usage{
sample.GDS(n.draws, log.phi,
		    post.mode,
		    fn.dens.post,
		    fn.dens.prop,
		    fn.draw.prop,
		    prop.params, ...,
		    max.tries=1000000,
		    report.freq=1,
		    announce=FALSE,
		    thread.id=1,
		    seed=.Random.seed
		    )
}
\arguments{
\item{n.draws}{number of draws to take from the target posterior density.}
\item{log.phi}{Vector of log.phi, as computed from the proposal draws.} 
\item{post.mode}{Mode of the target posterior density (numeric vector).}
\item{fn.dens.post}{Function that returns the log posterior density.  Function should take the parameter vector as the first argument.  Additional arguments are passed as ...}
\item{fn.dens.prop}{Function that returns the log density of the proposal distribution. The first argument of the function should take either a vector or a matrix.  If the argument is a matrix, each row is considered a sample.  Additional parameters are passed as a list, prop.params.}
\item{fn.draw.prop}{Function that returns random samples from the proposal density.  This function should return a matrix, with each row being a sample.  Additional parameters are passed as a list, prop.params.}
\item{prop.params}{Object (list or vector) to be passed to both fn.dens.prop and fn.draw.prop.Contains parameters for the proposal distribution.  See details.}
\item{...}{Additional parameters to be passed to fn.dens.post.}
\item{max.tries}{Maximum number of proposal draws to try, without a success.  This prevents the routine from being stuck in an endless loop.}
\item{report.freq}{The frequency that the function will report the current iteration.  For example, if report.freq=5, the function will display a message after every fifth iteration.}
\item{announce}{If TRUE, will print a message when a proposal is accepted as a sample from the target posterior distribution.}
\item{thread.id}{An identifier used in the announce function.  This is useful if running sample.GDS on multiple processors, to collect multiple batches of samples. Defaults to 1.}
\item{seed}{Sets a random seed within the call to sample.GDS.  Useful for assigning different seeds to calls to sample.GDS that are running on different threads or processors.  Defaults to .Random.seed.}
}
\value{
a list with the following elements:
\item{draws}{A matrix with each draw in a row, and each parameter in a column}
\item{counts}{The number of attempts that it took to get an accepted draw.  The accepted draw counts, so the count will always be at least 1.}
\item{gt.1}{A vector that indicates if the phi for that draw was greater than 1.  Available as a diagnostic.  Normally, these should all be FALSE.  Any values of TRUE suggest that a change in proposal density might be warranted.}
\item{log.post.dens}{A numeric vector.  Log posterior density for each draw.}
\item{log.prop.dens}{A numeric vector. Log of the proposal density for each draw.}
\item{log.thresholds}{Vector of threshold draws (log u) from the accept-reject algorithm.  Sorted in ascending order.}
\item{log.phi}{A numeric vector.  Value of log.phi for the accepted draws.}
}







 


