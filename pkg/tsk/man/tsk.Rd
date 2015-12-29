\name{tsk}
\alias{tsk}
\alias{tsk.default}
\alias{tsk.data.frame}
\alias{print.tskresult}
\title{Trimmed Spearman-Karber Method}
\description{Calculates the ED50 (effective dose) from a 
dose-response curve using the Trimmed Spearman-Karber method.}
\usage{
tsk(...)
\method{tsk}{default}(x, n, r, control = 0, trim = 0, 
    conf.level = 0.95, use.log.doses = TRUE,...)
\method{tsk}{data.frame}(input, control = 0, trim = 0, 
    conf.level = 0.95, use.log.doses = TRUE,...)
\method{print}{tskresult}(x,...)
}
\arguments{
  \item{x}{a non-empty numeric vector of doses.}
  \item{n}{a non-empty numeric vector of number of subjects at 
each dose, of the same length as x, or a number if the number of 
subjects are the same at each dose.}
  \item{r}{a non-empty numeric vector of the number of 
responding subjects at each dose.}
  \item{input}{Data frame containing the data to be analyzed. 
Must contain the three components x (doses), n (number of 
subjects at each dose), and r (number of subjects that respond).}
  \item{control}{fraction of subjects responding in a control 
experiment}
  \item{trim}{a number between 0 and .5. 0 gives the untrimmed 
    Spearman-Karber method.}
  \item{conf.level}{confidence level to use for confidence 
intervals.}
  \item{use.log.doses}{if TRUE then the log of the doses in x 
will be taken and used for calculation.}
  \item{...}{further arguments to be passed to or from methods.}
}
\value{A list of type "tskresult" containing the following components:
  \item{use.log.doses}{same as above.}
  \item{trim}{same as above.}
  \item{LD50}{The estimated effective dose for 50 percent of the 
population.}
  \item{gsd}{The geometric standard deviation of the LD50 
estimate.}
  \item{conf.level}{same as above.}
  \item{conf.int}{A confidence interval for the measurement, at 
the level of conf.level.}
If use.log.doses = FALSE, instead of LD50 and gsd, the list will 
contain:
\item{mu}{The estimated LD50, given a different name to prevent 
confusion with calculations on a log scale}
\item{sd}{The standard deviation of the mu estimate.}
}
\note{This code expects the data to have an increasing slope once smoothed. 
However, some studies may result in data with a decreasing trend. For 
instance, one user had data measuring the survival of freshwater fish under 
increasing water salinity: at low salinity all fish survive, while at high 
salinity all fish die. If using data with a decreasing trend, simply use the 
opposite response, e.g., swap the count that survive for the count that do 
not, or replace r with n-r.}
\note{This code was validated using the data given in Hamilton et 
al. (1977) Note that the confidence intervals given in Hamilton 
et al. (1977) are incorrect, as per Hamilton et al. (1978). This 
code replicates the results of the DOS program provided by the US 
EPA to perform the Spearman-Karber method, available at 
\url{http://www.epa.gov/eerd/stat2.htm}.}
\references{Hamilton, M. A.; Russo, R. C.; Thurston, R. V. 
  Trimmed Spearman-Karber Method for Estimating Median Lethal 
  Concentrations in Toxicity Bioassays. Enviro. Sci. Tech. 1977, 
  11 (7), 714-719. 
  \url{http://dx.doi.org/10.1021/es60130a004}

  Ibid, 1978, 12 (4), 417. 
  \url{http://dx.doi.org/10.1021/es60140a017}}
\author{Brenton R. Stone, \email{brstone@ufl.edu}}
\keyword{tsk}
\keyword{Spearman}
\keyword{Karber}
\keyword{ED50}
\keyword{LD50}
          \seealso{
            \code{\link{hamilton}}.
          }
\examples{tsk( c(1, 10, 100, 1000), 20, c(0, 3, 17, 20) )
dr <- data.frame( x=c(1, 10, 100, 1000), n=c(20, 20, 19, 20), 
                  r=c(0, 3, 17, 18) )
tsk( dr, 0, 0.1, 0.99 )

data(hamilton)
tsk(hamilton$dr1a)
}

