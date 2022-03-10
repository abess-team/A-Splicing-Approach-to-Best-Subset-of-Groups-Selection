library(Rcpp)
sourceCpp(file='gomp.cpp')

GOMP <- function(x, y, group, tol = 1e-5, iter.num = 50, ...)
{
  n <- nrow(x)
  p <- ncol(x)
  weight <- rep(1, n)
  gi <- unique(group)
  g.num <- length(gi)
  index <- match(gi, group)-1
  gsize <- as.vector(table(group))
  res <- gompCpp(x, y, weight = weight, N=g.num, n = n, p = p, index = index, gsize = gsize, tol = tol, iter_num = iter.num)
  iter <- res$iter
  beta <- res$beta[, 1:iter]
  intercept <- res$coef[1:iter]
  A <- res$A_out[1:iter]
  A <- lapply(A, function(z) z+1) 
  bic <- rep(0, iter)
  gic <- rep(0, iter)
  mse <- rep(0, iter)
  for (i in 1:iter) {
    gr.df <- sum(gsize[A[[i]]])
    mse[i] <- sum((y-x%*%beta[, i])^2)/n
    gic[i] <- n*log(mse[i])+0.5*log(g.num)*log(log(n))*gr.df
  }
  best.size <- which.min(gic)
  best.A <- A[[best.size]]
  best.beta <- beta[, best.size]
  best.intercept <- intercept[best.size]
  result <- list(beta = beta, intercept = intercept, A = A, mse = mse, bic = bic, gic = gic,
                 best.size = best.size, best.A = best.A, best.beta = best.beta, best.intercept = best.intercept)
  return(result)
}
