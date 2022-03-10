library(abess)
library(grpreg)
library(splines)
library(Ball)
load("trim32.rda")

## Data pre-processing
t <- 2000
z <- apply(x, 2, function(x) bcor(x, y))
id <- which(z >= sort(z, decreasing = T)[t])
temp <- x[, id]
X <- ns(temp[, 1], df = 5)
for (i in 2:ncol(temp)) {
  X <- cbind(X, ns(temp[, i], df = 5))
}
group = rep(1:ncol(temp), each = 5)
grp1 <- function(x, y, group, penalty){
  fit <- grpreg(x, y, group = group, penalty = penalty, eps = 1e-5)
  lam <- grpreg::select(fit, crit = "BIC")$lambda
  return(list(model = fit, lam = lam))
}
grp2 <- function(x, y, group, penalty){
  fit <- grpreg(x, y, group = group, penalty = penalty, gmax = 8, eps = 1e-5)
  lam <- grpreg::select(fit, crit = "BIC",df.method = "active")$lambda
  return(list(model = fit, lam = lam))
}

##############################################################################
## Real-world dataset analysis
##############################################################################

result <- sapply(101:200, function(i) {
  set.seed(i)
  print(i)
  ind <- sample(1:120, 100)
  m <- 120-length(ind)
  fit1 <- grp1(X[ind, ], y[ind], group = group, penalty="grLasso")
  l1 <- predict(fit1$model, type = "ngroups", lambda = fit1$lam)
  coef <- predict(fit1$model, type = "coef", lambda = fit1$lam)
  res1 <- y[-ind] - X[-ind, ]%*%coef[-1]-coef[1]

  fit2 <- grp2(X[ind, ], y[ind], group = group, penalty="grMCP")
  l2 <- predict(fit2$model, type = "ngroups", lambda = fit2$lam)
  coef <- predict(fit2$model, type = "coef", lambda = fit2$lam)
  res2 <- y[-ind] - X[-ind, ]%*%coef[-1]-coef[1]

  fit3 <- cv.grpreg(X[ind, ], y[ind], group = group, penalty="grLasso", seed = i, nfolds = 5)
  l3 <- predict(fit3, type = "ngroups")
  coef <- predict(fit3, type = "coef")
  res3 <- y[-ind] - X[-ind, ]%*%coef[-1]-coef[1]

  fit4 <- cv.grpreg(X[ind, ], y[ind], group = group, penalty="grMCP", seed = i, nfolds = 5)
  l4 <- predict(fit4, type = "ngroups")
  coef <- predict(fit4, type = "coef")
  res4 <- y[-ind] - X[-ind, ]%*%coef[-1]-coef[1]
  
  fit5 <- GOMP(X[ind, ], y[ind], group = group)
  l5 <- fit5$best.size
  coef <- c(fit5$intercept[l5], fit5$beta[, l5])
  res5 <- y[-ind] - X[-ind, ]%*%coef[-1]-coef[1]
  
  fit6 <- abess(X[ind, ], y[ind], group.index = group, tune.type = "gic",  type.path = "gsection", gs.range = c(1, 10))
  l6 <- fit6$best.size
  coef <- coef(fit6, support.size = l6)
  res6 <- y[-ind] - X[-ind, ]%*%coef[-1]-coef[1]
  
  fit7 <- abess(X[ind, ], y[ind], group.index = group, tune.type = "gic", support.size = 1:10)
  l7 <- fit7$best.size
  coef <- coef(fit7, support.size = fit7$best.size)
  res7 <- y[-ind] - X[-ind, ]%*%coef[-1]-coef[1]
  
  return(c(5*sum(res1^2),
           5*sum(res2^2),
           5*sum(res3^2),
           5*sum(res4^2),
           5*sum(res5^2),
           5*sum(res6^2),
           5*sum(res7^2),
           l1, l2, l3, l4, l5, l6, l7))})