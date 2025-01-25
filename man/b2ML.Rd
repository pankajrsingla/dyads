\name{b2ML}
\alias{b2ML}

\title{
  MCMC estimates for the (multilevel) b2 model
}
\description{
  Estimates the (multilevel) b2 model parameters, which is the bidirectional counterpart of the multilevel p2 model as described in Zijlstra, Van Duijn and Snijders (2006) <doi: 10.1027/1614-2241.2.1.42>.
}
\usage{
b2ML(nets, actor = NULL, density = NULL, adapt = NULL, burnin = NULL, center = NULL, 
separate= NULL, densVar = NULL, seed = NULL) 
}

\arguments{
  \item{nets}{List with n dichotomous symmetric dependent networks.}
  \item{actor}{Optional matrix with a stacked actor covariate, corresponding to the n networks. Multiple actor covariates can be added as a formula object, see example below}
  \item{density}{Optional matrix with symmetric a stacked density covariate, with dimensions similar to the n dependent networks. Multiple density covariates can be added as a formula object, see example below}
  \item{adapt}{Optional number of adaptive sequenses (default is 100).}
  \item{burnin}{Optional specification of number of burn-in iterations (default is 5000).}
  \item{center}{Optional argument for centering predictors (default is TRUE).}
  \item{separate}{Optional argument for estimating separate coefficients for the n dependent networks (default is FALSE).}
  \item{densVar}{Optional argument for estimating densty variance at the network level (default is TRUE).}
  \item{seed}{Optional specification of random seed (delfault is 1).}
}

\value{
Returns a fitted model of class b2ML, to be opened with the function summary().
}
\references{
Zijlstra, B. J., Van Duijn, M. A., & Snijders, T. A. (2006). The Multilevel p2 Model A random effects model for the analysis of multiple social networks. \emph{Methodology: European Journal of Research Methods for the Behavioral and Social Sciences}, 2(1), 42.
}
\author{
Bonne J.H. Zijlstra \email{b.j.h.zijlstra@uva.nl}
}

\examples{

# create two very small networks with covariates for illustrative purposes
Y1 <- matrix(c( 0,1,1,1,1,1,1,1,1,0,
                1,0,1,0,1,1,1,1,1,1,
                1,1,0,0,1,1,1,1,0,1,
                1,0,0,0,1,0,0,1,0,0,
                1,1,1,1,0,1,1,0,1,1,
                1,1,1,0,1,0,1,0,1,1,
                1,1,1,0,1,1,0,1,1,1,
                1,1,1,1,0,0,1,0,0,1,
                1,1,0,0,1,1,1,0,0,1,
                0,1,1,0,1,1,1,1,1,0), ncol=10) 
Y2 <- matrix(c( 0,0,1,0,1,1,0,1,0,0,
                0,0,0,0,0,0,0,1,1,0,
                1,0,0,1,0,1,0,1,0,0,
                0,0,1,0,0,0,1,1,0,0,
                1,0,0,0,0,0,1,1,0,0,
                1,0,1,0,0,0,1,1,0,0,
                0,0,0,1,1,1,0,1,0,0,
                1,1,1,1,1,1,1,0,0,1,
                0,1,0,0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,1,0,0), ncol=10) 
Y <- list(Y1, Y2)                
Aa1 <- c(1,0,1,0,1,1,0,1,0,1)
Aa2 <- c(1,0,0,1,0,0,1,1,0,1)
Aa <- list(Aa1, Aa2)
Aat <- do.call(plyr::rbind.fill.matrix, Aa)
Ab1 <- c(0,0,0,0,0,0,0,0,0,0)
Ab2 <- c(1,1,1,1,1,1,1,1,1,1)
Ab <- list(Ab1, Ab2)
Abt <- do.call(plyr::rbind.fill.matrix, Ab)
Da1 <- abs(matrix(rep(Aa1,10), byrow = FALSE, ncol= 10) -
            matrix(rep(Aa1,10), byrow = TRUE, ncol= 10))
Da2 <- abs(matrix(rep(Aa2,10), byrow = FALSE, ncol= 10) -
            matrix(rep(Aa2,10), byrow = TRUE, ncol= 10))
Da <- list(Da1, Da2)
Dat <- do.call(plyr::rbind.fill.matrix, Da)

# estimate b2 model for two networks
M1 <- b2ML(Y,actor= ~ Aat + Abt, density = ~ Dat, adapt = 10, burnin = 100, densVar = FALSE)
summary(M1)
# Notice: burn-in, and number of adaptive sequenses are 
# much smaller than recommended to keep computation time low.
# recommended code: 
\dontrun{
M1 <- b2ML(Y,actor= ~ Aat + Abt, density = ~ Dat, densVar = FALSE)
summary(M1)
}

# estimate b2 model for a single network
M2 <- b2ML(list(Y[[1]]),actor= ~ Aat[1:10,], density = ~ Dat[1:10,], adapt = 10, burnin = 100, 
densVar = FALSE)
summary(M2)
# Notice: burn-in, and number of adaptive sequenses are 
# much smaller than recommended to keep computation time low.
# recommended code: 
\dontrun{
M2 <- b2ML(list(Y[[1]]),actor= ~ Aat[1:10,], density = ~ Dat[1:10,], densVar = FALSE)
summary(M2)
}
}
\keyword{b2 model}

