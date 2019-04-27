#################################################################################################################
#
# START LIBRARY RAFAEL B. DE REZENDE - MAURO S. FERREIRA 5 FACTORS
#
#################################################################################################################

#' @title Estimation of the De Rezende-Ferreira 5 Factor model's parameters with variable \eqn{\ \tau}
#'
#' @description The command estimates the parameters of the De Rezende-Ferreira 5 Factor model
#'              using variable \eqn{\ \tau_{1}} and \eqn{\ \tau_{2}}
#'
#' @param   rate Vector or matrix of class "zoo", which contains interest rates
#' @param   maturity Vector of class "numeric", wich contains the maturities \cr
#'
#' @importFrom xts try.xts reclass
#' @importFrom stats coef na.omit lm resid time median optimize
#'
#' @details The De Rezende-Ferreira model used to fit the forward rates is: \cr
#'
#' \deqn{f_{t}\left(m\right) = \beta_{0t} +
#' \beta_{1t} e^{-\frac{m}{\tau_{1t}}} + \beta_{2t} e^{-\frac{m}{\tau_{2t}}} +
#' \beta_{3t} \left ( {\frac{m}{\tau_{1t}}} e^{-\frac{m}{\tau_{1t}}}  \right ) +
#' \beta_{4t} \left ( {\frac{m}{\tau_{2t}}} e^{-\frac{m}{\tau_{2t}}}  \right )} \cr
#'
#'
#' The spot rates, derived from the forward rates \eqn{f_{t}\left(m\right)}, are given by: \cr
#'
#' \deqn{y_{t}\left ( m \right ) = \beta_{0t} +
#' \beta_{1t} \left (\frac{1 - e^{-\frac{m}{\tau_{1t}}}}{\frac{m}{\tau_{1t}}}\right) +
#' \beta_{2t} \left (\frac{1 - e^{-\frac{m}{\tau_{2t}}}}{\frac{m}{\tau_{2t}}}\right) +
#' \beta_{3t} \left (\frac{1 - e^{-\frac{m}{\tau_{1t}}}}{\frac{m}{\tau_{1t}}} - e^{-\frac{m}{\tau_{1t}}}\right ) +
#' \beta_{4t} \left (\frac{1 - e^{-\frac{m}{\tau_{2t}}}}{\frac{m}{\tau_{2t}}} - e^{-\frac{m}{\tau_{2t}}}\right )}\cr
#'
#'
#' The set of optimal parameters will be chosen according to the lowest RMSE value:\cr
#'
#'   \deqn{\left (\widehat{\tau}_{1t},\widehat{\tau}_{2t} \right) = argmin\left \{\frac{1}{N}\sum_{t=1}^{N}
#'   \sqrt{\frac{1}{T}\sum_{t=1}^{T}\left [ y_{t}\left (t_{n} \right ) - \widehat{y}_{t}\left  (t_{n},\tau_{1t},
#'   \tau_{2t},\widehat{\beta_{t}}\right )\right ]^{2} } \right \}}\cr
#'
#' @return An object of class "zoo", that contains
#' \eqn{\ \left (\beta_{0t},\ \beta_{1t},\ \beta_{2t},\ \beta_{3t},\ \beta_{4t},\ \tau_{1t},\ \tau_{2t},\ SSR_{t},\ R^{2}_{t} \right)}
#'
#' @examples
#' #
#' # De Rezende-Ferreira 5F model on the Brazilian Data-Set
#' #
#'
#' data(ZC_Brazil)
#' real.rate = ZC_Brazil
#'
#' ZC_Brazil[["Date"]] = NULL
#'
#' rate = zoo(ZC_Brazil)
#' index(rate) = as.POSIXct(paste(real.rate[["Date"]]))
#' maturity <- c(0.5, 0.75, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
#'
#' RF.5F.Parameters <- DRF.5F.tVar(rate, maturity)
#'
#' par(mfrow=c(3,2))
#'  plot(RF.5F.Parameters[,"beta0"],xlab="Date",ylab="BETA0",ylim=c(9.5,12.0),col="blue",lwd=1)
#'  grid(nx=12, ny=12)
#'  plot(RF.5F.Parameters[,"beta1"],xlab="Date",ylab="BETA1",ylim=c(-18.0,2.3),col= "blue",lwd=1)
#'  grid(nx=12, ny=12)
#'  plot(RF.5F.Parameters[,"beta2"],xlab="Date",ylab="BETA2",ylim=c(-6.0,13.0),col= "blue",lwd=1)
#'  grid(nx=12, ny=12)
#'  plot(RF.5F.Parameters[,"beta3"],xlab="Date",ylab="BETA3",ylim=c(-10.0,0.0),col= "blue",lwd=1)
#'  grid(nx=12, ny=12)
#'  plot(RF.5F.Parameters[,"beta4"],xlab="Date",ylab="BETA4",ylim=c(-5.0,5.0),col="blue",lwd=1)
#'  grid(nx=12, ny=12)
#' par(mfrow=c(1,1))
#'
#' par(mfrow=c(2,1))
#'  plot(RF.5F.Parameters[,"tau1"],xlab="Date",ylab="TAU1",ylim=c(0.2,1.3),col="blue",lwd=1)
#'  grid(nx=12, ny=12)
#'  plot(RF.5F.Parameters[,"tau2"],xlab="Date",ylab="TAU2",ylim=c(2.5,5.5),col="blue",lwd=1)
#'  grid(nx=12, ny=12)
#' par(mfrow=c(1,1))
#'
#' @examples
#' #
#' # De Rezende-Ferreira 5F on the Russian Data-Set
#' #
#'
#' data(ZC_Russia)
#' real.rate = ZC_Russia
#'
#' ZC_Russia[["Date"]] = NULL
#'
#' rate = zoo(ZC_Russia)
#' index(rate) = as.POSIXct(paste(real.rate[["Date"]]))
#' maturity <- c(0.25, 0.5, 0.75, 1,2,3,5,7,10,15,20,30)
#' RF.5F.Parameters <- DRF.5F.tVar(rate, maturity)
#'
#' par(mfrow=c(3,2))
#'  plot(RF.5F.Parameters[,"beta0"],xlab="",ylab="BETA0",ylim=c(10.5,12.5),col="blue",lwd=1)
#'  grid(nx=12, ny=12)
#'  plot(RF.5F.Parameters[,"beta1"],xlab="Date",ylab="BETA1",ylim=c(-1.5,0.5),col="blue",lwd=1)
#'  grid(nx=12, ny=12)
#'  plot(RF.5F.Parameters[,"beta2"],xlab="Date",ylab="BETA2",ylim=c(-7.0,-3.5),col="blue",lwd=1)
#'  grid(nx=12, ny=12)
#'  plot(RF.5F.Parameters[,"beta3"],xlab="Date",ylab="BETA3",ylim=c(-1.5,3.5),col="blue",lwd=1)
#'  grid(nx=12, ny=12)
#'  plot(RF.5F.Parameters[,"beta4"],xlab="Date",ylab="BETA4",ylim=c(-5.5,-0.1),col="blue",lwd=1)
#'  grid(nx=12, ny=12)
#' par(mfrow=c(1,1))
#'
#'
#' par(mfrow=c(2,1))
#'  plot(RF.5F.Parameters[,"tau1"],xlab="Date",ylab="TAU1",ylim=c(0.1,1.9),col="blue",lwd=1)
#'  grid(nx=12, ny=12)
#'  plot(RF.5F.Parameters[,"tau2"],xlab="Date",ylab="TAU2",ylim=c(7.5,16.8),col="blue",lwd=1)
#'  grid(nx=12, ny=12)
#' par(mfrow=c(1,1))
#'
#' @export

DRF.5F.tVar <- function(rate, maturity )
{

  #
  # Setting of the parameters
  #

  start.Maturity1 = maturity[1]
  end.Maturity1 = median(maturity)
  step.Maturity1 = 0.75

  start.Maturity2 = median(maturity)
  end.Maturity2 = max(maturity)
  step.Maturity2 = 1

  ncolumn.tabResults = 9
  ncolumn.pos.SSR = 8

  dim.maturity <- length(maturity)

  #
  # Definition of the two subsets
  #

  Range.Maturity1 <- seq(start.Maturity1, end.Maturity1, by = step.Maturity1)
  Range.Maturity2 <- seq(start.Maturity2, end.Maturity2, by = step.Maturity2)

  #
  # Initialization of the matrix tabResults, including its columns names,
  # which will contain the parameters returned.
  #

  tabResults <- matrix(0, nrow(rate), ncolumn.tabResults)
  colnames(tabResults) <- c("beta0","beta1","beta2","beta3","beta4","tau1","tau2","SSR","R^2" )

  #
  # Initialization of the matrix TabTau which will contain the set of parameters characterized by
  # the lowest SSR obtained by fixing tau1 and varying tau2
  #

  TabTau <- matrix(0, length(Range.Maturity1), ncolumn.tabResults)

  #
  # Conversion of the rates, expressed in an arbitrary class, to xts class.
  # In case of error, a matrix object is used
  #

  rate <- try.xts(rate, error = as.matrix)

  #
  # External "while" cycle used to select the current date and the ZC bonds values to be fitted
  #

  j <- 1
  while (j <= nrow(rate))
  {

    #
    # Initialization of the matrix TempResults containing some partial results concerning the values
    # assumed by the parameters Beta, SSR and RQuadro obtained within the internal cycle "for"
    #

    TempResults <- matrix(0,length(Range.Maturity2), ncolumn.tabResults)

    #
    # Intermediate  cycle "for" used to select a maturity within Range.Maturity1
    #

    for (i in 1:length(Range.Maturity1))
    {

      #
      # The temporary optimal tau1 that maximise the factor loading FactLoad2 for each maturity within Range.Maturity1, are estimated
      #

      Tau1 <- optimize(FactLoad2,interval = c(0.001,max(Range.Maturity1)),maturity = Range.Maturity1[i],maximum = TRUE)$maximum

      #
      # Inner cycle "for" used to select a maturity within Range.Maturity2
      #

      for (a in 1:length(Range.Maturity2))
      {

        #
        # Fixed tau1, the temporary optimal tau2 that maximise the factor loading FactLoad2 for each maturity within Range.Maturity2, are estimated
        #

        Tau2 <- optimize(FactLoad2,interval = c(0.001,maturity[dim.maturity]),maturity = Range.Maturity2[a],maximum = TRUE)$maximum

        #
        # Beta parameters estimation using OLS given tau1 and tau2
        #

        EstParamOLS <- DRF.5F.estimator(as.numeric(rate[j,]), maturity, Tau1, Tau2)

        #
        # Update of the temporary table TempResults with the Beta, SSR and RQuadro values
        #

        Beta <- EstParamOLS$Par
        SSR <- sum(EstParamOLS$Res^2)
        RQuadro <- EstParamOLS$RQuadro
        TempResults[a,] <- c(Beta, Tau1, Tau2, SSR, RQuadro)
      }

      #
      # Fixed tau1, choice of the parameters record characterized by the
      # lowest SSR value with the associated tau2
      #

      BestRowTemp <- which.min( TempResults[,8])
      TabTau[i,] <- TempResults[BestRowTemp,]
    }


    #
    # Choise of the record with the lowst SSR value, associated to the current rate to fill
    # the table tabResults that will be returned to the Main Program
    #

    BestRow <- which.min(TabTau[,8])
    tabResults[j,] <- TabTau[BestRow,1:ncolumn.tabResults]

    #
    # Update of the date and bond yields record to be fitted
    #

    j <- j + 1
  }

  reclass(tabResults, rate)
}



#################################################################################################################

#' @title Estimation of the De Rezende-Ferreira 5 Factor model's parameters with fixed \eqn{\ \tau}
#'
#' @description The command estimates the parameters of the De Rezende-Ferreira 5 Factor model
#'              using fixed \eqn{\ \tau_{1}} and \eqn{\ \tau_{2}}
#'
#' @param   rate Vector or matrix of class "zoo", which contains interest rates
#' @param   maturity Vector of class "numeric", wich contains the maturities
#' @param   fixed_tau1 Decaying parameter of class "numeric" (Slope)
#' @param   fixed_tau2 Decaying parameter of class "numeric" (Curvature)
#'
#' @importFrom xts try.xts reclass
#' @importFrom stats coef na.omit lm resid time median  optimize
#'
#' @return An object of class "zoo", that contains
#' \eqn{\ \left (\beta_{0t},\ \beta_{1t},\ \beta_{2t},\ \beta_{3t},\ \beta_{4t},\ \tau_{1t},\ \tau_{2t},\ SSR_{t},\ R^{2}_{t} \right)}
#'
#' @examples
#' #
#' # De Rezende-Ferreira 5F model on the Indian Data-Set
#' #
#'
#' data(ZC_India)
#' real.rate = ZC_India
#' ZC_India[["Date"]] = NULL
#' rate = zoo(ZC_India)
#' index(rate) = as.POSIXct(paste(real.rate[["Date"]]))
#' maturity <- c(0.25, 0.5, 0.75, 1,2,3,4,5,6,7,8,9,10,12,15,20,25,30)
#' fixed_tau1 = (1.07612)
#' fixed_tau2 = (6.23293)
#'
#' RF.5F.Parameters <- DRF.5F.tFix(rate, maturity, fixed_tau1, fixed_tau2)
#'
#' par(mfrow=c(3,2))
#'  plot(RF.5F.Parameters[,"beta0"],xlab="Date",ylab="BETA0",ylim=c(7.0,9.0),col="blue",lwd=1)
#'  grid(nx=12, ny=12)
#'  plot(RF.5F.Parameters[,"beta1"],xlab="Date",ylab="BETA1",ylim=c(-3.5,0.2),col="blue",lwd=1)
#'  grid(nx=12, ny=12)
#'  plot(RF.5F.Parameters[,"beta2"],xlab="Date",ylab="BETA2",ylim=c(-1.5,1.0),col="blue",lwd=1)
#'  grid(nx=12, ny=12)
#'  plot(RF.5F.Parameters[,"beta3"],xlab="Date",ylab="BETA3",ylim=c(-2.0,0.5),col="blue",lwd=1)
#'  grid(nx=12, ny=12)
#'  plot(RF.5F.Parameters[,"beta4"],xlab="Date",ylab="BETA4",ylim=c(-2.5,5.0),col="blue",lwd=1)
#'  grid(nx=12, ny=12)
#' par(mfrow=c(1,1))
#'
#' @export

DRF.5F.tFix <- function(rate, maturity, fixed_tau1, fixed_tau2)
{

  #
  # Setting of the parameters
  #

  ncolumn.tabResults = 9

  #
  # Initialization of the matrix tabResults, including its columns names,
  # which will contain the parameters returned.
  #

  tabResults <- matrix(0, nrow(rate), ncolumn.tabResults)
  colnames( tabResults ) <- c("beta0","beta1","beta2","beta3","beta4","tau1","tau2","SSR","R^2" )

  #
  # Initialization of the decaying parameters tau1 and tau2
  # passed to the Main Program as arguments by value
  #

  Tau1 = fixed_tau1
  Tau2 = fixed_tau2

  #
  # Conversion of the rates, expressed in an arbitrary class, to xts class.
  # In case of error, a matrix object is used
  #

  rate <- try.xts(rate, error = as.matrix)

  #
  # External "while" cycle used to select the current date and the ZC bonds values to be fitted
  #

  j <- 1
  while (j <= nrow(rate) )
  {

    #
    # Beta parameters estimation using OLS
    #

    EstParamOLS <- DRF.5F.estimator(as.numeric(rate[j,]), maturity, Tau1, Tau2)

    #
    # Building of the object tabResults returned by the Function
    #

    Beta <- EstParamOLS$Par
    SSR <- sum(EstParamOLS$Res^2)
    RQuadro <- EstParamOLS$RQuadro

    tabResults[j,] <- c(Beta, Tau1, Tau2, SSR, RQuadro)


    #
    # Update of the date and bond yields record to be fitted
    #

    j <- j + 1
  }


  reclass(tabResults, rate)
}


#################################################################################################################

#' @title Estimation of spot rates with the De Rezende-Ferreira 5 Factor model
#'
#' @description The command estimates the spot rates using the De Rezende-Ferreira 5 Factor model
#'
#' @param   beta Matrix or Vector of class "zoo", which contains the coefficients of the De Rezende-Ferreira 5 Factor model:
#'               \eqn{\left( \beta_{0t},\beta_{1t}, \beta_{2t}, \beta_{3t}, \beta_{4t}, \tau_{1t}, \tau_{2t} \right) }
#'
#' @param   maturity Vector of class "numeric", wich contains the maturities  \cr
#'
#' @importFrom xts try.xts reclass xts
#' @importFrom stats coef na.omit lm resid time median  optimize
#'
#' @return  An object of class "xts" - "zoo", which contains fitted interest rates \cr
#'
#' @examples
#' #
#' # Fitting the Chinese spot rates using the De Rezende-Ferreira 5F moodel with Variable tau
#' #
#'
#' data(ZC_China)
#' real.rate = ZC_China
#' ZC_China[["Date"]] = NULL
#' rate = zoo(ZC_China)
#' index(rate) = as.POSIXct(paste(real.rate[["Date"]]))
#' maturity <- c(1,2,3,4,5,6,7,8,9,10,12,15,20,30)
#'
#' RF.5F.Parameters <- DRF.5F.tVar(rate, maturity)
#' RF.5F.Rates <- DRF.5F.rates(RF.5F.Parameters, maturity )
#'
#' plot(maturity,rate[5,],xlab="Maturity",ylab="Yields",ylim=c(3.5,4.7),col="black",lwd = 1)
#' lines(maturity, RF.5F.Rates[5,], col = "blue", lwd = 1)
#' grid(nx = 12, ny = 12)
#'
#' @examples
#' #
#' #
#' #
#' #
#' # Fitting the South African spot rates using the De Rezende-Ferreira 5F model with fixed tau
#' #
#'
#' data(ZC_SouthAfrica)
#' real.rate = ZC_SouthAfrica
#' ZC_SouthAfrica[["Date"]] = NULL
#' rate = zoo(ZC_SouthAfrica)
#' index(rate) = as.POSIXct(paste(real.rate[["Date"]]))
#' maturity <- c(0.25, 1,2,3,4,5,6,7,8,9,10,12,15,20,25,30)
#' fixed_tau1 = (1.07612)
#' fixed_tau2 = (6.23293)
#'
#' RF.5F.Parameters <- DRF.5F.tFix(rate, maturity, fixed_tau1, fixed_tau2)
#' RF.5F.Rates <- DRF.5F.rates(RF.5F.Parameters, maturity )
#'
#' plot(maturity,rate[5,],xlab="Maturity",ylab="Yields",ylim=c(6.5,10.0),col="black",lwd = 1)
#' lines(maturity, RF.5F.Rates[5,], col = "blue", lwd = 1)
#' grid(nx = 12, ny = 12)
#'
#' @export

DRF.5F.rates <- function(beta, maturity)
{

  #
  # Initialization of the YC and beta matrixes, which will contain respectively
  # the Yields Curves and the beta coefficients found applying the OLS
  #

  YC <- xts(matrix( 0, nrow(beta), length(maturity) ), order.by = time(beta))
  colnames(YC) <- make.names(maturity)
  beta <- as.matrix(beta)

  #
  # Building the fitted Yield Curves for each date
  #

  for (i in 1:nrow(beta))
  {
    YC[i,] <- beta[i,1] + beta[i,2] * FactLoad1( maturity, beta[i,6] ) + beta[i,3] * FactLoad1( maturity, beta[i,7] ) +
      beta[i,4] * FactLoad2( maturity, beta[i,6] ) + beta[i,5] * FactLoad2( maturity, beta[i,7] )
  }

  finalYC <- YC
  reclass(finalYC, beta)
}



#################################################################################################################
# Title:       Factor Loading 1 of the De Rezende-Ferreira 5 Factor model
#
# Description: The function defines the slope
#
# Parameters:  "maturity" --> Vector of class "numeric", wich contains the maturities
#              "tau"      --> Decaying parameter (Slope)
#
# Return:      FactLoad1
#

FactLoad1 <- function(maturity, tau)
{
  as.numeric( (1 - exp(-maturity/tau))/(maturity/tau))
}



#################################################################################################################

# Title:        Factor Loading 2 of the De Rezende-Ferreira 5 Factor model
#
# Description: The function defines the curvature
#
# Parameters:  "maturity" --> Vector of class "numeric", wich contains the maturities
#              "tau"      --> Decaying parameter (Curvature)
#
# Return:      FactLoad2
#

FactLoad2 <- function(maturity, tau)
{
  as.numeric(  ((1 - exp(-maturity/tau))/(maturity/tau) - exp(-maturity/tau)) )
}



#################################################################################################################

# Title:       Estimation of the parameters of the De Rezende-Ferreira 5 Factor model
#
# Description: The command estimates the parameters with OLS
#
# Parameters:  "rate"     -->  Vector or matrix of class "zoo", which contains interest rates
#              "maturity" -->  Vector of class "numeric", wich contains the maturities
#              "tau1"     -->  Decaying parameter (Slope)
#              "tau2"     -->  Decaying parameter (Curvature)
#
# ImportFrom:  stats --> (coef, na.omit, lm, resid, time, median, optimize)
#
# Return:      EstResults
#

DRF.5F.estimator <- function(rate, maturity, tau1, tau2)
{

  #
  # Beta parameters estimation using OLS
  #

  OLS.Results <- lm( rate ~ 1 + FactLoad1(maturity, tau1) + FactLoad1(maturity, tau2) +
                       FactLoad2(maturity, tau1) + FactLoad2(maturity, tau2) )

  #
  # Extraction of the parameters Beta, SSR and RQuadro from the OLS's Summary,
  # which will be returned to the Main Program
  #

  Beta <- coef(OLS.Results)

  #
  # Check for inconsistencies in the parameters
  #

  NaValues <- na.omit(Beta)
  if (length(NaValues) < 5) Beta <- c(0,0,0,0,0)

  #
  # Filling the fields of the returned object EstResults
  #

  names(Beta) <- c("beta0", "beta1", "beta2","beta3","beta4")
  Resid.Beta <- resid(OLS.Results)
  Temp.RQuadro <- summary(OLS.Results)$r.squared

  EstResults <- list(Par = Beta, Res = Resid.Beta, RQuadro = Temp.RQuadro)

  return(EstResults)
}



#
# END LIBRARY RAFAEL B. DE REZENDE - MAURO S. FERREIRA
#

