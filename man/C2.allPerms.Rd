\name{C2.allPerms}
\alias{C2.allPerms}

\title{
Contribution to C^2 for all variables
}
\description{
Returns a matrix with values for C^2_1 and C^2_k|C^2_k-1,C^2_k-2, ..., C^2_1, k=2,3, 4... for all possible permutations among k variables
}
\usage{
C2.allPerms(z, W, R)
}

\arguments{
  \item{z}{
observation vector, kx1
}
  \item{W}{
matrix of variables weigths, kxk
}
  \item{R}{
correlation matrix, kxk
}
}

\value{
Data frame where, the first k columns correspond to variable that entred the model first, second... k-th. The following (k+1) to 2*k columns contain the conditional contribution of the variable. The last column contains the sum of all contributions, meaning C^2_k
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

C2.allPerms(z = Z, W = Weights, R = Corr)

}
