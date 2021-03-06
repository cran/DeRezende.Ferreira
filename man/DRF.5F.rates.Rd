% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DeRezende.Ferreira.5F.R
\name{DRF.5F.rates}
\alias{DRF.5F.rates}
\title{Estimation of spot rates with the De Rezende-Ferreira 5 Factor model}
\usage{
DRF.5F.rates(beta, maturity)
}
\arguments{
\item{beta}{Matrix or Vector of class "zoo", which contains the coefficients of the De Rezende-Ferreira 5 Factor model:
\eqn{\left( \beta_{0t},\beta_{1t}, \beta_{2t}, \beta_{3t}, \beta_{4t}, \tau_{1t}, \tau_{2t} \right) }}

\item{maturity}{Vector of class "numeric", wich contains the maturities  \cr}
}
\value{
An object of class "xts" - "zoo", which contains fitted interest rates \cr
}
\description{
The command estimates the spot rates using the De Rezende-Ferreira 5 Factor model
}
\examples{
#
# Fitting the Chinese spot rates using the De Rezende-Ferreira 5F moodel with Variable tau
#

data(ZC_China)
real.rate = ZC_China
ZC_China[["Date"]] = NULL
rate = zoo(ZC_China)
index(rate) = as.POSIXct(paste(real.rate[["Date"]]))
maturity <- c(1,2,3,4,5,6,7,8,9,10,12,15,20,30)

RF.5F.Parameters <- DRF.5F.tVar(rate, maturity)
RF.5F.Rates <- DRF.5F.rates(RF.5F.Parameters, maturity )

plot(maturity,rate[5,],xlab="Maturity",ylab="Yields",ylim=c(3.5,4.7),col="black",lwd = 1)
lines(maturity, RF.5F.Rates[5,], col = "blue", lwd = 1)
grid(nx = 12, ny = 12)

#
#
#
#
# Fitting the South African spot rates using the De Rezende-Ferreira 5F model with fixed tau
#

data(ZC_SouthAfrica)
real.rate = ZC_SouthAfrica
ZC_SouthAfrica[["Date"]] = NULL
rate = zoo(ZC_SouthAfrica)
index(rate) = as.POSIXct(paste(real.rate[["Date"]]))
maturity <- c(0.25, 1,2,3,4,5,6,7,8,9,10,12,15,20,25,30)
fixed_tau1 = (1.07612)
fixed_tau2 = (6.23293)

RF.5F.Parameters <- DRF.5F.tFix(rate, maturity, fixed_tau1, fixed_tau2)
RF.5F.Rates <- DRF.5F.rates(RF.5F.Parameters, maturity )

plot(maturity,rate[5,],xlab="Maturity",ylab="Yields",ylim=c(6.5,10.0),col="black",lwd = 1)
lines(maturity, RF.5F.Rates[5,], col = "blue", lwd = 1)
grid(nx = 12, ny = 12)

}
