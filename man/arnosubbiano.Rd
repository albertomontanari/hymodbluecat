\name{arnosubbiano}
\alias{arnosubbiano}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sample dataset - Arno River basin - Daily data
}
\description{
The data set consists of daily time series for the Arno River basin at Subbiano, Italy. The basin has a drainage area of 751 km^2, with an elevation ranging from 259 m to 1622 m above sea level. The data set contains twenty-two years (i.e. from 1992-01-01 to 2013-12-31) of observed discharges, spatially averaged precipitation and potential evapotranspiration.}
\examples{
## Load data
data(arnosubbiano)

## Precipitation [mm/day]
precipit<-arnosubbiano[,2]

## Potential evapotranspiration [mm/day]
evapo<-arnosubbiano[,3]

## Observed discharge [m^3/s]
qoss<-arnosubbiano[,4]

## Plot the data
time<-arnosubbiano[,1]

par(mfrow=c(3,1),mar=c(4,5,2,1)) 
plot(time,precipit,ty="l",lwd=0.7,col="cyan",ylab="Precipitation [mm/day]")
title(main="Arno@Subbiano - Daily data",cex.main=1.2,col.main="darkblue")
plot(time,evapo,ty="l",lwd=0.7,col="orange",
ylab="Potential evapotranspiration [mm/day]")
plot(time,qoss,ty="l",lwd=0.7,col="blue",ylab="Obeserved discharge [m^3 s-1]")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
