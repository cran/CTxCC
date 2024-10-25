\name{wChisq.CLim}
\alias{wChisq.CLim}

\title{
Compute  control limit of Weighted Chi-Squared control charts for monitoring multivariate normal mean.
}
\description{
Computation of a control limit of the Weighted Chi-Squared control chart for a given vector of weights, w, correlation matrix, R, and the false alarm rate, alpha.
}
\usage{
wChisq.CLim(w,R,alpha)
}

\arguments{
  \item{w}{
vector of weigths, 1xk
}

 \item{R}{
correlation matrix, kxk
}

 \item{alpha}{
false alarm rate
}
}

\value{
\item{ContLim }{control limit of the Weighted Chi-Squared control chart}


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
# Table 1 in the Paper Criticality Assessment for Enhanced Multivariate Process Monitoring.

w <- c(0.29836,0.70164)  #vector of weights
R <- diag(2)
alpha <- 0.005
wChisq.CLim(w,R,alpha)

w <- c(0.23912,0.76088)  #vector of weights
R <- diag(2)
R[1,2] <- R[2,1] <- 0.25
alpha <- 0.005
wChisq.CLim(w,R,alpha)

}
