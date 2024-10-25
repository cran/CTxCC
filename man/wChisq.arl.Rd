\name{wChisq.arl}
\alias{wChisq.arl}

\title{
Compute ARLs of Weighted Chi-Squared control charts for monitoring multivariate normal mean.
}
\description{
Computation of the Average Run Length (ARL) for a Weighted Chi-Squared control chart for a given mean vector, delta, correlation matrix, R, control limit, h, and the vector of weights, w. The mean vector, delta, is defined in Propositon 4.2 from Paper Criticality Assessment for Enhanced Multivariate Process Monitoring.
}
\usage{
wChisq.arl(delta, R, h,  w)
}

\arguments{
\item{delta }{Vector of values representing the change in the mean for each variable, 1xk}
\item{R }{correlation matrix, kxk}
\item{h }{Control limit of Weighted Chi-Squared Control chart}
\item{w }{vector of weigths, 1xk}
}


\value{
\item{arl }{Average Run Length (ARL) for a Weighted Chi-Squared control chart for a given mean vector}

}
\references{
Paper
}
\author{
Dr. Burcu Aytaçoğlu (burcuaytacoglu@gmail.com)
Dr. Diana Barraza-Barraza (diana.barraza@ujed.mx),
Dr. Víctor G. Tercero-Gómez (victor.tercero@tec.mx),
Dr. A. Eduardo Cordero-Franco (lalo.cordero@gmail.com),

}

\examples{
#Table 1 in the Paper Criticality Assessment for Enhanced Multivariate Process Monitoring.
delta <- c(0.5, 0.5)    # mean vector (change vector)
R <- diag(2)    	# correlation matrix
h <- 2.649506  	# Control limit
w <- c(0.50153, 0.49847)	# vector of weights
wChisq.arl(delta, R, h,  w)

}
