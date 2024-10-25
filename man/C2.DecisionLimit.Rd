\name{C2.DecisionLimit}
\alias{C2.DecisionLimit}

\title{
Conditional decision limit for z, given x already in model
}
\description{
Calculates the conditional decision limit for z, given x already in model, using the exact distribution for the conditional contribution of z to C_k
}
\usage{
C2.DecisionLimit(z, mu.C, R.C, A, x.var, alpha)
}

\arguments{
  \item{z}{
observation vector, kx1, where z[x.var, ] correspond to variables already in the model
}
  \item{mu.C}{
scalar, conditional mean for z given x
}
  \item{R.C}{
scalar, conditional covariance for z given x
}
  \item{A}{
list containing matrix decomposition of A, preferably, obtained from function decomposeA
}
  \item{x.var}{
vector indicating variables already present in the model. length: k-1
}
  \item{alpha}{
confidence level for decision limit
}
}
\details{
Proposition Distribution of a C^2 contribution from Paper Criticality Assessment for Enhanced Multivariate Process Monitoring
}
\value{
conditionalCL, conditional decision limit for z's contribution to C_k

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
k<-3
sigma0 = matrix(diag(rep(1,k)),ncol = k)
mu0 = matrix(c(0,0,0), ncol = 1)
Weights = diag(c(0.5, 0.25,0.25))

library(mvtnorm)
set.seed(1000)
X = matrix(ncol= 1, data = rmvnorm(n = 1, mean = mu0, sigma = sigma0))
Z = (X - mu0)/sqrt(as.numeric(diag(sigma0)))
Corr<-get.R(Sigma0 = sigma0)

A<-decomposeA(W = Weights, R = Corr, x.var = 1:2, z.var = 3)

Par<-zConditionalParameters(mean0 = mu0, R0 = Corr, z = Z, x.var = 1:2, z.var = 3)
C2.DecisionLimit(z = Z, mu.C = Par$muC, R.C = Par$RC, A = A, x.var = 1:2, alpha = 0.95)


}

