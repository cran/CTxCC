\name{ms.optim}
\alias{ms.optim}
\title{Multistart Weight Optimization for Design of a CTX Control Chart}
\description{
Computes the optimal weight vector and control limit for a CTX control chart
that minimizes the overall out-of-control Average Run Length (ARL) across
a given set of shift scenarios.
Optimization is performed using a multistart approach.
A quasi-Newton method (\code{"BFGS"}) is used by default, with the Nelder–Mead algorithm
as a backup in multistart optimization.
}
\usage{
ms.optim(delta, R, alpha, n_starts,
method_primary, method_secondary,
scale_sum,jitter_sd, seed)
}
\arguments{
  \item{delta}{List. A set of shift vectors to be evaluated in the optimization.}
  \item{R}{Numeric matrix. Correlation (or covariance) matrix for the CTX control chart.}
  \item{alpha}{Numeric. Significance level for the in-control ARL, where
    \code{alpha = 1/ARL0}.}
  \item{n_starts}{Integer. Number of random starting points for multistart optimization.}
  \item{method_primary}{Character. Primary optimization method.}
  \item{method_secondary}{Character. Secondary (backup) optimization method used if the
    primary method fails.}
  \item{scale_sum}{Numeric. Scaling constant applied to normalize the weight vector.}
  \item{jitter_sd}{Numeric. Standard deviation of Gaussian noise used to jitter
    starting points in multistart optimization.}
  \item{seed}{Integer. Random seed for reproducibility.}
}
\details{
The optimization is performed over log-transformed proxy weights to ensure
positivity, and the resulting weights are normalized to sum to \code{scale_sum}.
The objective function minimizes the sum of out-of-control ARLs across all
specified shift scenarios. The control limit is computed to satisfy the
desired in-control ARL, given by \code{alpha = 1/ARL0}.
}
\value{

  \item{weights}{Optimized weight vector.}
  \item{control_limit}{Computed control limit.}
  \item{objective}{Optimal value of the objective function.}
  \item{method}{Optimization method used.}
  \item{n_starts}{Number of starting points used in the optimization.}

}

\references{
Paper (pending complete bibliographic reference when available).

}
\author{
Dr. Burcu Aytaçoğlu (\email{burcuaytacoglu@gmail.com}),
Dr. Diana Barraza-Barraza (\email{diana.barraza@ujed.mx}),
Dr. Víctor G. Tercero-Gómez (\email{victor.tercero@tec.mx}),
Dr. A. Eduardo Cordero-Franco (\email{lalo.cordero@gmail.com})
}

\examples{
\dontrun{
# Define shift scenarios

R <- diag(3)
delta <- list(c(0.5, 0, 0), c(0, 0.5, 0), c(0, 0, 0.5))

# Run multi-start optimization
res1 <- ms.optim(delta = delta, R = R, alpha = 0.005, n_starts = 10)

# Run multistart optimization with jittered initial points
res2 <- ms.optim(delta = delta, R = R, alpha = 0.005, n_starts = 5,
                 method_primary = "BFGS", method_secondary = "Nelder-Mead")
}

}

