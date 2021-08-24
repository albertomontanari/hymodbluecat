\name{hymod.par}
\alias{hymod.par}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to fit the parameter set of the Hymod rainfall-runoff model (Boyle, 2000).}
\description{
It calls a Fortran subroutine to make the actual computation of the simulated river flows. It requires the libraries "devtools", "DescTools" and "DEoptim" if the optimisation routine DEoptim is used (see below).}
\usage{
hymod.par(param.ini,area,tdelta,e,p,qoss,nstep=length(p),qinitial=0,lower=c(10,0.1,0,1,1),upper=c(400,10,0.9,1000,1000),itermax=100,control=list(factr=1e1,fnscale=0.01,parscale=c(100,1,1,1,1)),opt="DEoptim",plot=T))}

\arguments{
\item{param.ini}{Vector of the initial value for the 5 Hymod parameters Cmax, beta, alpha, ks, kq. Suitable values may be Cmax=200 mm, beta=0.5, alpha=0.5, ks=Tc where Tc is the time of concentration of the catchment (expressed in number of simulation time steps) and kq=Tc/2.}
\item{area}{The catchment contributing area given in square kilometers.}
\item{Tdelta}{The model time step given in seconds.}
\item{e}{Vector of potential evapotranspiration. It should have the same length as p.}
\item{p}{Vector of mean areal rainfall. It should have the same length as e.}
\item{qoss}{Vector of observed river flow in cubic meters per second. It should have the same length as e.}
\item{nstep}{Number of simulation time steps. Default is the sample size of e and p. It can be lower but not higher.}
\item{qinitial}{Initial value of river flow. Default is zero.}
\item{lower}{Lower bound for the model parameter vector. Default is c(10,0.1,0,1,1).}
\item{upper}{Upper bound for the model parameter vector. Default is c(400,10,0.9,1000,1000).}
\item{itermax}{Maximum number of trials allowed during optimisation. It used only by DEoptim (see below). Default is 100.}
\item{control}{Control parameters used by DEoptim. See the help of DEoptim for additional details.}
\item{opt}{Specifies the optimisation routine to be used. Allowed options are "DEoptim" and "optim" for the DEoptim routine (requires the library DEoptim) and the optim routine (available in R-base), respectively. See the related helps for additional details. Default is "DEoptim".}
\item{plot}{Specifies if diagnostic is to be returned. If plot=T (the default) a scapperplot of observed versus optimally simulated river flows is returned along with the related Nash efficiency (Nash and Sutcliffe, 1968).}}

\value{The output from DEoptim (or optim) is returned (see the related helps) along with these additional list items:
\item{outputlist$par}{The vector of optimised parameters. It is advisable to check that the resulting parameter values are not close to lower and upper limits. If they are, limits can be adjusted accordingly.}
\item{outputlist$qoss}{The vector of observed calibration data.}
\item{outputlist$qsim}{The vector of simulated calibration data.}
\item{outputlist$eff}{The Nash Sutcliffe efficiency of the best parameter set.}
}

\details{
HyMOD is a five-parameter conceptual rainfall-runoff model that was introduced by Boyle (2000). The model is based on the probability-distributed soil storage capacity principle introduced by Moore (1985).
In HyMOD, the rainfall-runoff process is represented through a nonlinear tank connected with three identical linear tanks in parallel representing the surface flow and a slow-flow tank representing groundwater flow. The model requires the optimization of five parameters: Cmax (the maximum storage capacity within the watershed), beta (the degree of spatial variability of the soil moisture capacity within the watershed), alpha (a factor for partitioning the flow between two series of tanks) and the two residence time parameters of quick-flow and slow-flow tanks, kq and ks respectively. 
The input data consist of precipitation and potential evapotranspiration at the given time scale.}

\author{Demetris Koutsoyiannis and Alberto Montanari, email: alberto.montanari@unibo.it.}

\references{
Boyle, D.P., (2000), Multicriteria calibration of hydrological models, Ph.D. dissertation, Dep. of Hydrol. and Water Resour., Univ of Arizona, Tucson.

Moore, R.J., (1985), The probability-distributed principle and runoff production at point and basin scale, Hydrol. Sci. J., 30(2), 273-297.

Nash, J. E., and Sutcliffe, J. V., (1970), River flow forecasting through conceptual models part I – A discussion of principles, J. Hydrol., 10, 282–290, \doi{https://doi.org/10.1016/0022-1694(70)90255-6}. 
}

\examples{
## Load data
data(arnosubbiano)
e=arnosubbiano[,3]
p=arnosubbiano[,2]
qoss=arnosubbiano[,4]
## Run optimisation
pr1=hymod.par(c(100,1,0.5,200,0.5),area=752,tdelta=86400,e,p,qoss,nstep=length(p),qinitial=15,lower=c(10,0.1,0.1,0.1,0.1),upper=c(800,10,0.9,1000,100),opt="DEoptim")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Hydrology}
\keyword{Rainfall-runoff}
\keyword{Hymod}% __ONLY ONE__ keyword per line
