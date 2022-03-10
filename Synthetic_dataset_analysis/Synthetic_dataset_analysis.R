library(abess)
library(grpreg)
library(stringr)
library(mccr)
library(mvnfast)

## Evaluation criteria
f <- function(y, pre){
  tp <- sum(pre==1 & y==1)
  fp <- sum(pre==1 & y!=1)
  fn <- sum(pre!=1 & y==1)
  tn <- sum(pre!=1 & y!=1)
  mcc <- mccr(y, pre)
  result <- c(tp+fp, tp/(tp+fn), fp/(tn+fp), mcc)
  return(result)
}

## Generate coefficients
gen.coef <- function(k){
  coef <- rnorm(k+1, 0, 1)
  coef <- (coef - mean(coef))[-1]
  return(coef)
}
coeff <- function(Tn, k, ind, p){
  coef <- sapply(rep(k, Tn), gen.coef)
  beta <- rep(0, p)
  index <- as.vector(sapply(ind, function(x){
    return(((x-1)*k+1):(x*k))
  }))
  beta[index] <- coef
  return(beta)
}

## Settings for GLasso and GMCP
grp1 <- function(x, y, group, penalty){
  fit <- grpreg(x, y, group = group, penalty = penalty, eps = 1e-5)
  lam <- grpreg::select(fit, crit = "BIC")$lambda
  return(list(model = fit, lam = lam))
}
grp2 <- function(x, y, group, penalty){
  fit <- grpreg(x, y, group = group, penalty = penalty, gmax = 50, eps = 1e-5)
  lam <- grpreg::select(fit, crit = "BIC", df.method = "active")$lambda
  return(list(model = fit, lam = lam))
}

## Extract the selected groups from the output of abess
detect.group <- function(v, k){
  v <- as.numeric(unlist(stringr::str_extract_all(v, "[0-9]+")))
  temp <- v[which(v%%k == 1)]
  return((temp-1)/k+1)
}



##############################################################################
## Simulation 1:Influence of the correlation across groups
##############################################################################

n <- 500
p <- 4500
k <- 3
Tn <- 15

result <- sapply(1:100, function(i){
  print(i)
  set.seed(i)
  
  ind <- sort(sample(1:(p/k), Tn))
  beta <- coeff(Tn, k, ind, p)
  
  R <- matrix(rnorm(p*n, 0, 1), n, p)
  sigma <- matrix(0, (p/k), (p/k))
  for (i in 1:(p/k)) {
    for (j in 1:(p/k)) {
      sigma[i, j] <- 0.3^(abs(i-j))
      #sigma[i, j] <- ifelse(i==j, 1, 0.3)
    }
  }
  Z <- rmvn(n, rep(0, (p/k)), sigma)
  z <- sapply(1:p, function(x) {
    g <- floor((x-1)/k)+1
    return((Z[, g]+R[, x])/sqrt(2))
  })
  x <- matrix(unlist(z), n)
  pe <- x %*% beta
  y <- pe+rnorm(n, 0, 2)
  group <- rep(1:(p/k), each=k)
  true <- rep(0, length(unique(group)))
  true[ind] <- 1
  t <- as.numeric(system.time(fit <- grp1(x, y, group = group, penalty = "grLasso"))[3])
  lasso <- rep(0, length(unique(group)))
  gr <- predict(fit$model, type = "groups", lambda = fit$lam)
  coef <- predict(fit$model, type = "coef", lambda = fit$lam)
  lasso[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  lasso1 <- c(f(true, lasso), sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- grp2(x, y, group = group, penalty = "grMCP"))[3])
  MCP <- rep(0, length(unique(group)))
  gr <- predict(fit$model, type = "groups", lambda = fit$lam)
  coef <- predict(fit$model, type = "coef", lambda = fit$lam)
  MCP[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  MCP1 <- c(f(true, MCP), sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- cv.grpreg(x, y, group = group, penalty = "grLasso", seed = i, nfolds = 5))[3])
  lasso <- rep(0, length(unique(group)))
  gr <- predict(fit, type = "groups")
  coef <- predict(fit, type = "coef")
  lasso[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  lasso1.cv <- c(f(true, lasso), sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- cv.grpreg(x, y, group = group, penalty = "grMCP", seed = i, nfolds = 5))[3])
  MCP <- rep(0, length(unique(group)))
  gr <- predict(fit, type = "groups")
  coef <- predict(fit, type = "coef")
  MCP[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  MCP1.cv <- c(f(true, MCP), sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- GOMP(x, y,  group = group))[3])
  coef <- c(fit$best.intercept, fit$best.beta)
  best.group <- fit$best.A
  gomp <- rep(0, length(unique(group)))
  gomp[best.group] <- 1
  gomp1 <- c(f(true, gomp), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- abess(x, y,  group.index = group, tune.type = "gic"))[3])
  coef <- coef(fit, support.size = fit$best.size)
  best.group <- detect.group(extract(fit)[["support.vars"]], k)
  gsplicing <- rep(0, length(unique(group)))
  gsplicing[best.group] <- 1
  abess1 <- c(f(true, gsplicing), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- abess(x, y,  group.index = group, tune.type = "gic", tune.path = "gsection"))[3])
  coef <- coef(fit, support.size = fit$best.size)
  best.group <- detect.group(extract(fit)[["support.vars"]], k)
  gsplicing <- rep(0, length(unique(group)))
  gsplicing[best.group] <- 1
  abess.gs1 <- c(f(true, gsplicing), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
  
  R <- matrix(rnorm(p*n, 0, 1), n, p)
  sigma <- matrix(0, (p/k), (p/k))
  for (i in 1:(p/k)) {
    for (j in 1:(p/k)) {
      sigma[i, j] <- 0.6^(abs(i-j))
      #sigma[i, j] <- ifelse(i==j, 1, 0.6)
    }
  }
  Z <- rmvn(n, rep(0, (p/k)), sigma)
  z <- sapply(1:p, function(x) {
    g <- floor((x-1)/k)+1
    return((Z[, g]+R[, x])/sqrt(2))
  })
  x <- matrix(unlist(z), n)
  pe <- x %*% beta
  y <- pe+rnorm(n, 0, 2)
  group <- rep(1:(p/k), each=k)
  true <- rep(0, length(unique(group)))
  true[ind] <- 1
  t <- as.numeric(system.time(fit <- grp1(x, y, group = group, penalty = "grLasso"))[3])
  lasso <- rep(0, length(unique(group)))
  gr <- predict(fit$model, type = "groups", lambda = fit$lam)
  coef <- predict(fit$model, type = "coef", lambda = fit$lam)
  lasso[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  lasso2 <- c(f(true, lasso), sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- grp2(x, y, group = group, penalty = "grMCP"))[3])
  MCP <- rep(0, length(unique(group)))
  gr <- predict(fit$model, type = "groups", lambda = fit$lam)
  coef <- predict(fit$model, type = "coef", lambda = fit$lam)
  MCP[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  MCP2 <- c(f(true, MCP), sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- cv.grpreg(x, y, group = group, penalty = "grLasso", seed = i, nfolds = 5))[3])
  lasso <- rep(0, length(unique(group)))
  gr <- predict(fit, type = "groups")
  coef <- predict(fit, type = "coef")
  lasso[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  lasso2.cv <- c(f(true, lasso), sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- cv.grpreg(x, y, group = group, penalty = "grMCP", seed = i, nfolds = 5))[3])
  MCP <- rep(0, length(unique(group)))
  gr <- predict(fit, type = "groups")
  coef <- predict(fit, type = "coef")
  MCP[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  MCP2.cv <- c(f(true, MCP),  sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- GOMP(x, y,  group = group))[3])
  coef <- c(fit$best.intercept, fit$best.beta)
  best.group <- fit$best.A
  gomp <- rep(0, length(unique(group)))
  gomp[best.group] <- 1
  gomp2 <- c(f(true, gomp), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- abess(x, y,  group.index = group, tune.type = "gic"))[3])
  coef <- coef(fit, support.size = fit$best.size)
  best.group <- detect.group(extract(fit)[["support.vars"]], k)
  gsplicing <- rep(0, length(unique(group)))
  gsplicing[best.group] <- 1
  abess2 <- c(f(true, gsplicing), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- abess(x, y,  group.index = group, tune.type = "gic", tune.path = "gsection"))[3])
  coef <- coef(fit, support.size = fit$best.size)
  best.group <- detect.group(extract(fit)[["support.vars"]], k)
  gsplicing <- rep(0, length(unique(group)))
  gsplicing[best.group] <- 1
  abess.gs2 <- c(f(true, gsplicing), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
  
  R <- matrix(rnorm(p*n, 0, 1), n, p)
  sigma <- matrix(0, (p/k), (p/k))
  for (i in 1:(p/k)) {
    for (j in 1:(p/k)) {
      sigma[i, j] <- 0.9^(abs(i-j))
      #sigma[i, j] <- ifelse(i==j, 1, 0.9)
    }
  }
  Z <- rmvn(n, rep(0, (p/k)), sigma)
  z <- sapply(1:p, function(x) {
    g <- floor((x-1)/k)+1
    return((Z[, g]+R[, x])/sqrt(2))
  })
  x <- matrix(unlist(z), n)
  pe <- x %*% beta
  y <- pe+rnorm(n, 0, 2)
  group <- rep(1:(p/k), each=k)
  true <- rep(0, length(unique(group)))
  true[ind] <- 1
  t <- as.numeric(system.time(fit <- grp1(x, y, group = group, penalty = "grLasso"))[3])
  lasso <- rep(0, length(unique(group)))
  gr <- predict(fit$model, type = "groups", lambda = fit$lam)
  coef <- predict(fit$model, type = "coef", lambda = fit$lam)
  lasso[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  lasso3 <- c(f(true, lasso), sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- grp2(x, y, group = group, penalty = "grMCP"))[3])
  MCP <- rep(0, length(unique(group)))
  gr <- predict(fit$model, type = "groups", lambda = fit$lam)
  coef <- predict(fit$model, type = "coef", lambda = fit$lam)
  MCP[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  MCP3 <- c(f(true, MCP),  sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- cv.grpreg(x, y, group = group, penalty = "grLasso", seed = i, nfolds = 5))[3])
  lasso <- rep(0, length(unique(group)))
  gr <- predict(fit, type = "groups")
  coef <- predict(fit, type = "coef")
  lasso[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  lasso3.cv <- c(f(true, lasso),  sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- cv.grpreg(x, y, group = group, penalty = "grMCP", seed = i, nfolds = 5))[3])
  MCP <- rep(0, length(unique(group)))
  gr <- predict(fit, type = "groups")
  coef <- predict(fit, type = "coef")
  MCP[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  MCP3.cv <- c(f(true, MCP),  sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- GOMP(x, y,  group = group))[3])
  coef <- c(fit$best.intercept, fit$best.beta)
  best.group <- fit$best.A
  gomp <- rep(0, length(unique(group)))
  gomp[best.group] <- 1
  gomp3 <- c(f(true, gomp), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- abess(x, y,  group.index = group, tune.type = "gic"))[3])
  coef <- coef(fit, support.size = fit$best.size)
  best.group <- detect.group(extract(fit)[["support.vars"]], k)
  gsplicing <- rep(0, length(unique(group)))
  gsplicing[best.group] <- 1
  abess3 <- c(f(true, gsplicing), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- abess(x, y,  group.index = group, tune.type = "gic", tune.path = "gsection"))[3])
  coef <- coef(fit, support.size = fit$best.size)
  best.group <- detect.group(extract(fit)[["support.vars"]], k)
  gsplicing <- rep(0, length(unique(group)))
  gsplicing[best.group] <- 1
  abess.gs3 <- c(f(true, gsplicing), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)

    return(c(lasso1, MCP1, lasso1.cv, MCP1.cv, gomp1, abess1, abess.gs1,
           lasso2, MCP2, lasso2.cv, MCP2.cv, gomp2, abess2, abess.gs2,
           lasso3, MCP3, lasso3.cv, MCP3.cv, gomp3, abess3, abess.gs3))})



##############################################################################
## Simulation 2:Influence of the sample size
##############################################################################

n1 <- 500
n2 <- 600
n3 <- 700
p <- 5000
k <- 5
Tn <- 10

result <- sapply(1:100, function(i){
  print(i)
  set.seed(i)
  
  ind <- sort(sample(1:(p/k), Tn))
  beta <- coeff(Tn, k, ind, p)
  R <- matrix(rnorm(p*n1, 0, 1), n1, p)
  sigma <- matrix(0, (p/k), (p/k))
  for (i in 1:(p/k)) {
    for (j in 1:(p/k)) {
      sigma[i, j] <- ifelse(i==j, 1, 0.6)
    }
  }
  Z <- rmvn(n1, rep(0, (p/k)), sigma)
  z <- sapply(1:p, function(x) {
    g <- floor((x-1)/k)+1
    return((Z[, g]+R[, x])/sqrt(2))
  })
  x <- matrix(unlist(z), n1)
  pe <- x %*% beta
  y <- pe+rnorm(n1, 0, 3)
  group <- rep(1:(p/k), each=k)
  true <- rep(0, length(unique(group)))
  true[ind] <- 1
  t <- as.numeric(system.time(fit <- grp1(x, y, group = group, penalty = "grLasso"))[3])
  lasso <- rep(0, length(unique(group)))
  gr <- predict(fit$model, type = "groups", lambda = fit$lam)
  coef <- predict(fit$model, type = "coef", lambda = fit$lam)
  lasso[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  lasso1 <- c(f(true, lasso), sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- grp2(x, y, group = group, penalty = "grMCP"))[3])
  MCP <- rep(0, length(unique(group)))
  gr <- predict(fit$model, type = "groups", lambda = fit$lam)
  coef <- predict(fit$model, type = "coef", lambda = fit$lam)
  MCP[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  MCP1 <- c(f(true, MCP), sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- cv.grpreg(x, y, group = group, penalty = "grLasso", seed = i, nfolds = 5))[3])
  lasso <- rep(0, length(unique(group)))
  gr <- predict(fit, type = "groups")
  coef <- predict(fit, type = "coef")
  lasso[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  lasso1.cv <- c(f(true, lasso),  sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- cv.grpreg(x, y, group = group, penalty = "grMCP", seed = i, nfolds = 5))[3])
  MCP <- rep(0, length(unique(group)))
  gr <- predict(fit, type = "groups")
  coef <- predict(fit, type = "coef")
  MCP[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  MCP1.cv <- c(f(true, MCP),  sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- GOMP(x, y,  group = group))[3])
  coef <- c(fit$best.intercept, fit$best.beta)
  best.group <- fit$best.A
  gomp <- rep(0, length(unique(group)))
  gomp[best.group] <- 1
  gomp1 <- c(f(true, gomp), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- abess(x, y,  group.index = group, tune.type = "gic"))[3])
  coef <- coef(fit, support.size = fit$best.size)
  best.group <- detect.group(extract(fit)[["support.vars"]], k)
  gsplicing <- rep(0, length(unique(group)))
  gsplicing[best.group] <- 1
  abess1 <- c(f(true, gsplicing), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- abess(x, y,  group.index = group, tune.type = "gic", tune.path = "gsection"))[3])
  coef <- coef(fit, support.size = fit$best.size)
  best.group <- detect.group(extract(fit)[["support.vars"]], k)
  gsplicing <- rep(0, length(unique(group)))
  gsplicing[best.group] <- 1
  abess.gs1 <- c(f(true, gsplicing), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
  
  
  R <- matrix(rnorm(p*n2, 0, 1), n2, p)
  Z <- rmvn(n2, rep(0, (p/k)), sigma)
  z <- sapply(1:p, function(x) {
    g <- floor((x-1)/k)+1
    return((Z[, g]+R[, x])/sqrt(2))
  })
  x <- matrix(unlist(z), n2)
  pe <- x %*% beta
  y <- pe+rnorm(n2, 0, 3)
  group <- rep(1:(p/k), each=k)
  true <- rep(0, length(unique(group)))
  true[ind] <- 1
  t <- as.numeric(system.time(fit <- grp1(x, y, group = group, penalty = "grLasso"))[3])
  lasso <- rep(0, length(unique(group)))
  gr <- predict(fit$model, type = "groups", lambda = fit$lam)
  coef <- predict(fit$model, type = "coef", lambda = fit$lam)
  lasso[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  lasso2 <- c(f(true, lasso), sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- grp2(x, y, group = group, penalty = "grMCP"))[3])
  MCP <- rep(0, length(unique(group)))
  gr <- predict(fit$model, type = "groups", lambda = fit$lam)
  coef <- predict(fit$model, type = "coef", lambda = fit$lam)
  MCP[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  MCP2 <- c(f(true, MCP), sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- cv.grpreg(x, y, group = group, penalty = "grLasso", seed = i, nfolds = 5))[3])
  lasso <- rep(0, length(unique(group)))
  gr <- predict(fit, type = "groups")
  coef <- predict(fit, type = "coef")
  lasso[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  lasso2.cv <- c(f(true, lasso),  sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- cv.grpreg(x, y, group = group, penalty = "grMCP", seed = i, nfolds = 5))[3])
  MCP <- rep(0, length(unique(group)))
  gr <- predict(fit, type = "groups")
  coef <- predict(fit, type = "coef")
  MCP[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  MCP2.cv <- c(f(true, MCP),  sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- GOMP(x, y,  group = group))[3])
  coef <- c(fit$best.intercept, fit$best.beta)
  best.group <- fit$best.A
  gomp <- rep(0, length(unique(group)))
  gomp[best.group] <- 1
  gomp2 <- c(f(true, gomp), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
  
  
  t <- as.numeric(system.time(fit <- abess(x, y,  group.index = group, tune.type = "gic"))[3])
  coef <- coef(fit, support.size = fit$best.size)
  best.group <- detect.group(extract(fit)[["support.vars"]], k)
  gsplicing <- rep(0, length(unique(group)))
  gsplicing[best.group] <- 1
  abess2 <- c(f(true, gsplicing), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- abess(x, y,  group.index = group, tune.type = "gic", tune.path = "gsection"))[3])
  coef <- coef(fit, support.size = fit$best.size)
  best.group <- detect.group(extract(fit)[["support.vars"]], k)
  gsplicing <- rep(0, length(unique(group)))
  gsplicing[best.group] <- 1
  abess.gs2 <- c(f(true, gsplicing), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
  
  R <- matrix(rnorm(p*n3, 0, 1), n3, p)
  Z <- rmvn(n3, rep(0, (p/k)), sigma)
  z <- sapply(1:p, function(x) {
    g <- floor((x-1)/k)+1
    return((Z[, g]+R[, x])/sqrt(2))
  })
  x <- matrix(unlist(z), n3)
  pe <- x %*% beta
  y <- pe+rnorm(n3, 0, 3)
  group <- rep(1:(p/k), each=k)
  true <- rep(0, length(unique(group)))
  true[ind] <- 1
  t <- as.numeric(system.time(fit <- grp1(x, y, group = group, penalty = "grLasso"))[3])
  lasso <- rep(0, length(unique(group)))
  gr <- predict(fit$model, type = "groups", lambda = fit$lam)
  coef <- predict(fit$model, type = "coef", lambda = fit$lam)
  lasso[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  lasso3 <- c(f(true, lasso), sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- grp2(x, y, group = group, penalty = "grMCP"))[3])
  MCP <- rep(0, length(unique(group)))
  gr <- predict(fit$model, type = "groups", lambda = fit$lam)
  coef <- predict(fit$model, type = "coef", lambda = fit$lam)
  MCP[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  MCP3 <- c(f(true, MCP), sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- cv.grpreg(x, y, group = group, penalty = "grLasso", seed = i, nfolds = 5))[3])
  lasso <- rep(0, length(unique(group)))
  gr <- predict(fit, type = "groups")
  coef <- predict(fit, type = "coef")
  lasso[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  lasso3.cv <- c(f(true, lasso),  sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- cv.grpreg(x, y, group = group, penalty = "grMCP", seed = i, nfolds = 5))[3])
  MCP <- rep(0, length(unique(group)))
  gr <- predict(fit, type = "groups")
  coef <- predict(fit, type = "coef")
  MCP[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  MCP3.cv <- c(f(true, MCP),  sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- GOMP(x, y,  group = group))[3])
  coef <- c(fit$best.intercept, fit$best.beta)
  best.group <- fit$best.A
  gomp <- rep(0, length(unique(group)))
  gomp[best.group] <- 1
  gomp3 <- c(f(true, gomp), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- abess(x, y,  group.index = group, tune.type = "gic"))[3])
  coef <- coef(fit, support.size = fit$best.size)
  best.group <- detect.group(extract(fit)[["support.vars"]], k)
  gsplicing <- rep(0, length(unique(group)))
  gsplicing[best.group] <- 1
  abess3 <- c(f(true, gsplicing), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- abess(x, y,  group.index = group, tune.type = "gic", tune.path = "gsection"))[3])
  coef <- coef(fit, support.size = fit$best.size)
  best.group <- detect.group(extract(fit)[["support.vars"]], k)
  gsplicing <- rep(0, length(unique(group)))
  gsplicing[best.group] <- 1
  abess.gs3 <- c(f(true, gsplicing), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
  return(c(lasso1, MCP1, lasso1.cv, MCP1.cv, gomp1, abess1, abess.gs1,
           lasso2, MCP2, lasso2.cv, MCP2.cv, gomp2, abess2, abess.gs2,
           lasso3, MCP3, lasso3.cv, MCP3.cv, gomp3, abess3, abess.gs3))})



##############################################################################
## Simulation 3:Influence of the number of groups
##############################################################################

n <- 500
p1 <- 2000
p2 <- 4000
p3 <- 6000
k <- 4
Tn <- 10

result <- sapply(1:100, function(i){
  print(i)
  set.seed(i)
  
  ind <- sort(sample(1:(p1/k), Tn))
  coefff <- sapply(rep(k, Tn), gen.coef)
  beta <- rep(0, p1)
  index <- as.vector(sapply(ind, function(x){
    return(((x-1)*k+1):(x*k))
  }))
  beta[index] <- coefff
  R <- matrix(rnorm(p1*n, 0, 1), n, p1)
  sigma <- matrix(0, (p1/k), (p1/k))
  for (i in 1:(p1/k)) {
    for (j in 1:(p1/k)) {
      sigma[i, j] <- 0.9^(abs(i-j))
    }
  }
  Z <- rmvn(n, rep(0, (p1/k)), sigma)
  z <- sapply(1:p1, function(x) {
    g <- floor((x-1)/k)+1
    return((Z[, g]+R[, x])/sqrt(2))
  })
  x <- matrix(unlist(z), n)
  y <- x %*% beta+rnorm(n, 0, 3)
  group <- rep(1:(p1/k), each=k)
  true <- rep(0, length(unique(group)))
  true[ind] <- 1
  t <- as.numeric(system.time(fit <- grp1(x, y, group = group, penalty = "grLasso"))[3])
  lasso <- rep(0, length(unique(group)))
  gr <- predict(fit$model, type = "groups", lambda = fit$lam)
  coef <- predict(fit$model, type = "coef", lambda = fit$lam)
  lasso[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  lasso1 <- c(f(true, lasso), sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- grp2(x, y, group = group, penalty = "grMCP"))[3])
  MCP <- rep(0, length(unique(group)))
  gr <- predict(fit$model, type = "groups", lambda = fit$lam)
  coef <- predict(fit$model, type = "coef", lambda = fit$lam)
  MCP[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  MCP1 <- c(f(true, MCP), sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- cv.grpreg(x, y, group = group, penalty = "grLasso", seed = i, nfolds = 5))[3])
  lasso <- rep(0, length(unique(group)))
  gr <- predict(fit, type = "groups")
  coef <- predict(fit, type = "coef")
  lasso[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  lasso1.cv <- c(f(true, lasso), sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- cv.grpreg(x, y, group = group, penalty = "grMCP", seed = i, nfolds = 5))[3])
  MCP <- rep(0, length(unique(group)))
  gr <- predict(fit, type = "groups")
  coef <- predict(fit, type = "coef")
  MCP[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  MCP1.cv <- c(f(true, MCP), sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- GOMP(x, y,  group = group))[3])
  coef <- c(fit$best.intercept, fit$best.beta)
  best.group <- fit$best.A
  gomp <- rep(0, length(unique(group)))
  gomp[best.group] <- 1
  gomp1 <- c(f(true, gomp), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- abess(x, y,  group.index = group, tune.type = "gic"))[3])
  coef <- coef(fit, support.size = fit$best.size)
  best.group <- detect.group(extract(fit)[["support.vars"]], k)
  gsplicing <- rep(0, length(unique(group)))
  gsplicing[best.group] <- 1
  abess1 <- c(f(true, gsplicing), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- abess(x, y,  group.index = group, tune.type = "gic", tune.path = "gsection"))[3])
  coef <- coef(fit, support.size = fit$best.size)
  best.group <- detect.group(extract(fit)[["support.vars"]], k)
  gsplicing <- rep(0, length(unique(group)))
  gsplicing[best.group] <- 1
  abess.gs1 <- c(f(true, gsplicing), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
  
  ind <- sort(sample(1:(p2/k), Tn))
  beta <- rep(0, p2)
  index <- as.vector(sapply(ind, function(x){
    return(((x-1)*k+1):(x*k))
  }))
  beta[index] <- coefff
  R <- matrix(rnorm(p2*n, 0, 1), n, p2)
  sigma <- matrix(0, (p2/k), (p2/k))
  for (i in 1:(p2/k)) {
    for (j in 1:(p2/k)) {
      sigma[i, j] <- 0.9^(abs(i-j))
    }
  }
  Z <- rmvn(n, rep(0, (p2/k)), sigma)
  z <- sapply(1:p2, function(x) {
    g <- floor((x-1)/k)+1
    return((Z[, g]+R[, x])/sqrt(2))
  })
  x <- matrix(unlist(z), n)
  pe <- x %*% beta
  y <- pe+rnorm(n, 0, 3)
  group <- rep(1:(p2/k), each=k)
  true <- rep(0, length(unique(group)))
  true[ind] <- 1
  t <- as.numeric(system.time(fit <- grp1(x, y, group = group, penalty = "grLasso"))[3])
  lasso <- rep(0, length(unique(group)))
  gr <- predict(fit$model, type = "groups", lambda = fit$lam)
  coef <- predict(fit$model, type = "coef", lambda = fit$lam)
  lasso[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  lasso2 <- c(f(true, lasso), sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- grp2(x, y, group = group, penalty = "grMCP"))[3])
  MCP <- rep(0, length(unique(group)))
  gr <- predict(fit$model, type = "groups", lambda = fit$lam)
  coef <- predict(fit$model, type = "coef", lambda = fit$lam)
  MCP[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  MCP2 <- c(f(true, MCP), sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- cv.grpreg(x, y, group = group, penalty = "grLasso", seed = i, nfolds = 5))[3])
  lasso <- rep(0, length(unique(group)))
  gr <- predict(fit, type = "groups")
  coef <- predict(fit, type = "coef")
  lasso[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  lasso2.cv <- c(f(true, lasso),  sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- cv.grpreg(x, y, group = group, penalty = "grMCP", seed = i, nfolds = 5))[3])
  MCP <- rep(0, length(unique(group)))
  gr <- predict(fit, type = "groups")
  coef <- predict(fit, type = "coef")
  MCP[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  MCP2.cv <- c(f(true, MCP),  sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- GOMP(x, y,  group = group))[3])
  coef <- c(fit$best.intercept, fit$best.beta)
  best.group <- fit$best.A
  gomp <- rep(0, length(unique(group)))
  gomp[best.group] <- 1
  gomp2 <- c(f(true, gomp), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- abess(x, y,  group.index = group, tune.type = "gic"))[3])
  coef <- coef(fit, support.size = fit$best.size)
  best.group <- detect.group(extract(fit)[["support.vars"]], k)
  gsplicing <- rep(0, length(unique(group)))
  gsplicing[best.group] <- 1
  abess2 <- c(f(true, gsplicing), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- abess(x, y,  group.index = group, tune.type = "gic", tune.path = "gsection"))[3])
  coef <- coef(fit, support.size = fit$best.size)
  best.group <- detect.group(extract(fit)[["support.vars"]], k)
  gsplicing <- rep(0, length(unique(group)))
  gsplicing[best.group] <- 1
  abess.gs2 <- c(f(true, gsplicing), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
  
  ind <- sort(sample(1:(p3/k), Tn))
  beta <- rep(0, p3)
  index <- as.vector(sapply(ind, function(x){
    return(((x-1)*k+1):(x*k))
  }))
  beta[index] <- coefff
  ind <- sort(sample(1:(p3/k), Tn))
  beta <- coeff(Tn, k, ind, p3)
  R <- matrix(rnorm(p3*n, 0, 1), n, p3)
  sigma <- matrix(0, (p3/k), (p3/k))
  for (i in 1:(p3/k)) {
    for (j in 1:(p3/k)) {
      sigma[i, j] <- 0.9^(abs(i-j))
    }
  }
  Z <- rmvn(n, rep(0, (p3/k)), sigma)
  z <- sapply(1:p3, function(x) {
    g <- floor((x-1)/k)+1
    return((Z[, g]+R[, x])/sqrt(2))
  })
  x <- matrix(unlist(z), n)
  pe <- x %*% beta
  y <- pe+rnorm(n, 0, 3)
  group <- rep(1:(p3/k), each=k)
  true <- rep(0, length(unique(group)))
  true[ind] <- 1
  t <- as.numeric(system.time(fit <- grp1(x, y, group = group, penalty = "grLasso"))[3])
  lasso <- rep(0, length(unique(group)))
  gr <- predict(fit$model, type = "groups", lambda = fit$lam)
  coef <- predict(fit$model, type = "coef", lambda = fit$lam)
  lasso[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  lasso3 <- c(f(true, lasso), sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- grp2(x, y, group = group, penalty = "grMCP"))[3])
  MCP <- rep(0, length(unique(group)))
  gr <- predict(fit$model, type = "groups", lambda = fit$lam)
  coef <- predict(fit$model, type = "coef", lambda = fit$lam)
  MCP[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  MCP3 <- c(f(true, MCP), sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- cv.grpreg(x, y, group = group, penalty = "grLasso", seed = i, nfolds = 5))[3])
  lasso <- rep(0, length(unique(group)))
  gr <- predict(fit, type = "groups")
  coef <- predict(fit, type = "coef")
  lasso[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  lasso3.cv <- c(f(true, lasso),  sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- cv.grpreg(x, y, group = group, penalty = "grMCP", seed = i, nfolds = 5))[3])
  MCP <- rep(0, length(unique(group)))
  gr <- predict(fit, type = "groups")
  coef <- predict(fit, type = "coef")
  MCP[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  MCP3.cv <- c(f(true, MCP),  sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- GOMP(x, y,  group = group))[3])
  coef <- c(fit$best.intercept, fit$best.beta)
  best.group <- fit$best.A
  gomp <- rep(0, length(unique(group)))
  gomp[best.group] <- 1
  gomp3 <- c(f(true, gomp), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- abess(x, y,  group.index = group, tune.type = "gic"))[3])
  coef <- coef(fit, support.size = fit$best.size)
  best.group <- detect.group(extract(fit)[["support.vars"]], k)
  gsplicing <- rep(0, length(unique(group)))
  gsplicing[best.group] <- 1
  abess3 <- c(f(true, gsplicing), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- abess(x, y,  group.index = group, tune.type = "gic", tune.path = "gsection"))[3])
  coef <- coef(fit, support.size = fit$best.size)
  best.group <- detect.group(extract(fit)[["support.vars"]], k)
  gsplicing <- rep(0, length(unique(group)))
  gsplicing[best.group] <- 1
  abess.gs3 <- c(f(true, gsplicing), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
  
  return(c(lasso1, MCP1, lasso1.cv, MCP1.cv, gomp1, abess1, abess.gs1,
           lasso2, MCP2, lasso2.cv, MCP2.cv, gomp2, abess2, abess.gs2,
           lasso3, MCP3, lasso3.cv, MCP3.cv, gomp3, abess3, abess.gs3))})



##############################################################################
## Simulation 4:Influence of the group size
##############################################################################

n <- 1000
J <- 1000
k1 <- 5
k2 <- 10
k3 <- 15
Tn <- 5

result <- sapply(1:100, function(i){
  print(i)
  set.seed(i)
  
  ind <- sort(sample(1:J, Tn))
  coefff <- sapply(rep(k1, Tn), gen.coef)
  beta <- rep(0, J*k1)
  index <- as.vector(sapply(ind, function(x){
    return(((x-1)*k1+1):(x*k1))
  }))
  beta[index] <- coefff
  R <- matrix(rnorm(J*k1*n, 0, 1), n, (J*k1))
  sigma <- matrix(0, J, J)
  for (i in 1:J) {
    for (j in 1:J) {
      sigma[i, j] <- ifelse(i==j, 1, 0.9)
    }
  }
  Z <- rmvn(n, rep(0, J), sigma)
  z <- sapply(1:(J*k1), function(x) {
    g <- floor((x-1)/k1)+1
    return((Z[, g]+R[, x])/sqrt(2))
  })
  x <- matrix(unlist(z), n)
  y <- x %*% beta+rnorm(n, 0, 5)
  group <- rep(1:J, each=k1)
  true <- rep(0, length(unique(group)))
  true[ind] <- 1
  t <- as.numeric(system.time(fit <- grp1(x, y, group = group, penalty = "grLasso"))[3])
  lasso <- rep(0, length(unique(group)))
  gr <- predict(fit$model, type = "groups", lambda = fit$lam)
  coef <- predict(fit$model, type = "coef", lambda = fit$lam)
  lasso[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  lasso1 <- c(f(true, lasso), sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- grp2(x, y, group = group, penalty = "grMCP"))[3])
  MCP <- rep(0, length(unique(group)))
  gr <- predict(fit$model, type = "groups", lambda = fit$lam)
  coef <- predict(fit$model, type = "coef", lambda = fit$lam)
  MCP[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  MCP1 <- c(f(true, MCP), sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- cv.grpreg(x, y, group = group, penalty = "grLasso", seed = i, nfolds = 5))[3])
  lasso <- rep(0, length(unique(group)))
  gr <- predict(fit, type = "groups")
  coef <- predict(fit, type = "coef")
  lasso[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  lasso1.cv <- c(f(true, lasso), sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- cv.grpreg(x, y, group = group, penalty = "grMCP", seed = i, nfolds = 5))[3])
  MCP <- rep(0, length(unique(group)))
  gr <- predict(fit, type = "groups")
  coef <- predict(fit, type = "coef")
  MCP[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  MCP1.cv <- c(f(true, MCP), sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- GOMP(x, y,  group = group))[3])
  coef <- c(fit$best.intercept, fit$best.beta)
  best.group <- fit$best.A
  gomp <- rep(0, length(unique(group)))
  gomp[best.group] <- 1
  gomp1 <- c(f(true, gomp), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- abess(x, y,  group.index = group, tune.type = "gic"))[3])
  coef <- coef(fit, support.size = fit$best.size)
  best.group <- detect.group(extract(fit)[["support.vars"]], k1)
  gsplicing <- rep(0, length(unique(group)))
  gsplicing[best.group] <- 1
  abess1 <- c(f(true, gsplicing), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- abess(x, y,  group.index = group, tune.type = "gic", tune.path = "gsection"))[3])
  coef <- coef(fit, support.size = fit$best.size)
  best.group <- detect.group(extract(fit)[["support.vars"]], k1)
  gsplicing <- rep(0, length(unique(group)))
  gsplicing[best.group] <- 1
  abess.gs1 <- c(f(true, gsplicing), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
  
  coefff <- sapply(rep(k2, Tn), gen.coef)
  beta <- rep(0, J*k2)
  index <- as.vector(sapply(ind, function(x){
    return(((x-1)*k2+1):(x*k2))
  }))
  beta[index] <- coefff
  R <- matrix(rnorm(J*k2*n, 0, 1), n, (J*k2))
  z <- sapply(1:(J*k2), function(x) {
    g <- floor((x-1)/k2)+1
    return((Z[, g]+R[, x])/sqrt(2))
  })
  x <- matrix(unlist(z), n)
  y <- x %*% beta+rnorm(n, 0, 5)
  group <- rep(1:J, each=k2)
  true <- rep(0, length(unique(group)))
  true[ind] <- 1
  t <- as.numeric(system.time(fit <- grp1(x, y, group = group, penalty = "grLasso"))[3])
  lasso <- rep(0, length(unique(group)))
  gr <- predict(fit$model, type = "groups", lambda = fit$lam)
  coef <- predict(fit$model, type = "coef", lambda = fit$lam)
  lasso[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  lasso2 <- c(f(true, lasso), sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- grp2(x, y, group = group, penalty = "grMCP"))[3])
  MCP <- rep(0, length(unique(group)))
  gr <- predict(fit$model, type = "groups", lambda = fit$lam)
  coef <- predict(fit$model, type = "coef", lambda = fit$lam)
  MCP[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  MCP2 <- c(f(true, MCP), sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- cv.grpreg(x, y, group = group, penalty = "grLasso", seed = i, nfolds = 5))[3])
  lasso <- rep(0, length(unique(group)))
  gr <- predict(fit, type = "groups")
  coef <- predict(fit, type = "coef")
  lasso[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  lasso2.cv <- c(f(true, lasso), sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- cv.grpreg(x, y, group = group, penalty = "grMCP", seed = i, nfolds = 5))[3])
  MCP <- rep(0, length(unique(group)))
  gr <- predict(fit, type = "groups")
  coef <- predict(fit, type = "coef")
  MCP[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  MCP2.cv <- c(f(true, MCP),  sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- GOMP(x, y,  group = group))[3])
  coef <- c(fit$best.intercept, fit$best.beta)
  best.group <- fit$best.A
  gomp <- rep(0, length(unique(group)))
  gomp[best.group] <- 1
  gomp2 <- c(f(true, gomp), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- abess(x, y,  group.index = group, tune.type = "gic"))[3])
  coef <- coef(fit, support.size = fit$best.size)
  best.group <- detect.group(extract(fit)[["support.vars"]], k2)
  gsplicing <- rep(0, length(unique(group)))
  gsplicing[best.group] <- 1
  abess2 <- c(f(true, gsplicing), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- abess(x, y,  group.index = group, tune.type = "gic", tune.path = "gsection"))[3])
  coef <- coef(fit, support.size = fit$best.size)
  best.group <- detect.group(extract(fit)[["support.vars"]], k2)
  gsplicing <- rep(0, length(unique(group)))
  gsplicing[best.group] <- 1
  abess.gs2 <- c(f(true, gsplicing), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
  
  coefff <- sapply(rep(k3, Tn), gen.coef)
  beta <- rep(0, J*k3)
  index <- as.vector(sapply(ind, function(x){
    return(((x-1)*k3+1):(x*k3))
  }))
  beta[index] <- coefff
  R <- matrix(rnorm(J*k3*n, 0, 1), n, (J*k3))
  z <- sapply(1:(J*k3), function(x) {
    g <- floor((x-1)/k3)+1
    return((Z[, g]+R[, x])/sqrt(2))
  })
  x <- matrix(unlist(z), n)
  y <- x %*% beta+rnorm(n, 0, 5)
  group <- rep(1:J, each=k3)
  true <- rep(0, length(unique(group)))
  true[ind] <- 1
  t <- as.numeric(system.time(fit <- grp1(x, y, group = group, penalty = "grLasso"))[3])
  lasso <- rep(0, length(unique(group)))
  gr <- predict(fit$model, type = "groups", lambda = fit$lam)
  coef <- predict(fit$model, type = "coef", lambda = fit$lam)
  lasso[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  lasso3 <- c(f(true, lasso), sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- grp2(x, y, group = group, penalty = "grMCP"))[3])
  MCP <- rep(0, length(unique(group)))
  gr <- predict(fit$model, type = "groups", lambda = fit$lam)
  coef <- predict(fit$model, type = "coef", lambda = fit$lam)
  MCP[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  MCP3 <- c(f(true, MCP), sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- cv.grpreg(x, y, group = group, penalty = "grLasso", seed = i, nfolds = 5))[3])
  lasso <- rep(0, length(unique(group)))
  gr <- predict(fit, type = "groups")
  coef <- predict(fit, type = "coef")
  lasso[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  lasso3.cv <- c(f(true, lasso),  sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- cv.grpreg(x, y, group = group, penalty = "grMCP", seed = i, nfolds = 5))[3])
  MCP <- rep(0, length(unique(group)))
  gr <- predict(fit, type = "groups")
  coef <- predict(fit, type = "coef")
  MCP[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  MCP3.cv <- c(f(true, MCP),  sqrt(sum(sum((beta1-beta)^2)+intercept^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- GOMP(x, y,  group = group))[3])
  coef <- c(fit$best.intercept, fit$best.beta)
  best.group <- fit$best.A
  gomp <- rep(0, length(unique(group)))
  gomp[best.group] <- 1
  gomp3 <- c(f(true, gomp), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- abess(x, y,  group.index = group, tune.type = "gic"))[3])
  coef <- coef(fit, support.size = fit$best.size)
  best.group <- detect.group(extract(fit)[["support.vars"]], k3)
  gsplicing <- rep(0, length(unique(group)))
  gsplicing[best.group] <- 1
  abess3 <- c(f(true, gsplicing), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
  
  t <- as.numeric(system.time(fit <- abess(x, y,  group.index = group, tune.type = "gic", tune.path = "gsection"))[3])
  coef <- coef(fit, support.size = fit$best.size)
  best.group <- detect.group(extract(fit)[["support.vars"]], k3)
  gsplicing <- rep(0, length(unique(group)))
  gsplicing[best.group] <- 1
  abess.gs3 <- c(f(true, gsplicing), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
  
  return(c(lasso1, MCP1, lasso1.cv, MCP1.cv, gomp1, abess1, abess.gs1,
           lasso2, MCP2, lasso2.cv, MCP2.cv, gomp2, abess2, abess.gs2,
           lasso3, MCP3, lasso3.cv, MCP3.cv, gomp3, abess3, abess.gs3))})