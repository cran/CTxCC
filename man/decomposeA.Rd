\name{decomposeA}
\alias{decomposeA}

\title{
Calculation and decomposition of matrix A
}
\description{
Decomposition of matrix A, required in Proposition 4.3. Decomoposition given by equation 41
}
\usage{
decomposeA(W, R, x.var, z.var)
}

\arguments{
  \item{W}{
diagonal matrix containing the corresponding weigth for each monitored variable. Dimensions kxk
}
  \item{R}{
correlation matrix for monitores variables, kxk
}
  \item{x.var}{
vector indicating variables already present in the model. length: k-1.
}
  \item{z.var}{
scalar indicating variables to be included.
}
}
\details{
Note that length(z.var) + length(x.var) = k
}
\value{
Returns decomposition of matrix A according to Equation 41 in paper.
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
k<-6 # variables
B<-matrix(runif(n = k*k),ncol= k)### creating random matrix for sigma0
sigma0 <- B\%*\%t(B)
R<-get.R(sigma0)
Weights = diag(rep(1/k,k))
decomposeA(W = Weights, R = R, x.var = 1:5, z.var = 6)
}
