\name{sievefornacina}
\alias{sievefornacina}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sample dataset - Sieve River basin - Hourly data
}
\description{
The data set consists of hourly time series for the Sieve River basin, a tributary of
the Arno River in Central Italy. The basin has a drainage area of 830 km^2, with an elevation ranging from 92 m to 1020 m above sea level. The data set contains three years (i.e. from 1992-01-01 to 1994-12-31) of observed discharges, spatially averaged precipitation and potential evapotranspiration.}
\examples{
## Load data
data(sievefornacina)

## Precipitation [mm/hour]
precipit<-sievefornacina[,2]

## Potential evapotranspiration [mm/hour]
evapo<-sievefornacina[,3]

## Observed discharge [m^3/s]
qoss<-sievefornacina[,4]

## Plot the data
time<-sievefornacina[,1]

par(mfrow=c(3,1),mar=c(4,5,2,1)) 
plot(time,precipit,ty="l",lwd=0.7,col="cyan",ylab="Precipitation [mm/hour]")
title(main="Sieve@Fornacina - Hourly data",cex.main=1.2,col.main="darkblue")
plot(time,evapo,ty="l",lwd=0.7,col="orange",ylab="Potential evapotranspiration [mm/hour]")
plot(time,qoss,ty="l",lwd=0.7,col="blue",ylab="Obeserved discharge [m^3 s-1]")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
