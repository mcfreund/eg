## finding faster function for making CI ----
# https://www.alexejgossmann.com/benchmarking_r/

library(microbenchmark)
set.seed(2017)
n <- 1000
x <- rnorm(n)
X <- cbind(rep(1, length(x)), x)  ## add intercept
y <- rnorm(n)

check_for_equal_coefs <- function(values) {
  tol <- 1e-12
  max_error <- max(
    c(
      abs(values[[1]] - values[[2]]),
      abs(values[[2]] - values[[3]]),
      abs(values[[1]] - values[[3]]),
      abs(values[[1]] - values[[4]]),
      abs(values[[1]] - values[[5]])
    )
  )
  max_error < tol
}

mbm <- microbenchmark(
  
  "lm" = {
    
    y.hat <- predict(lm(y ~ x))
    
    },
  
  "cor" = {
    
    b1 <- cor(x, y) * sd(y) / sd(x)
    b0 <- mean(y) - b1 * mean(x)
    y.hat <- b1 * x + b0
    
    },
  
  "svd" = {
    
    y.hat <- svd(X)$u %*% t(svd(X)$u) %*% y
    
    },
  
  "pinv" = {
    
    y.hat <- X %*% solve(t(X) %*% X) %*% t(X) %*% y

    },
  
  "lin.sys" = {
    y.hat <- c(X %*% solve(t(X) %*% X, t(X) %*% y))
    
    },
  
  check = check_for_equal_coefs
  
  )

mbm  ## lin.sys wins! (but cor close behind)
library(ggplot2)
autoplot(mbm)

## test within sapply() loop ----

## first validate:

x <- rnorm(100)
y <- rnorm(100)
data <- data.frame(x = x, y = y)
grid <- data.frame(x = data$x)
set.seed(0)
p1 <- sapply(
  seq_len(10),
  function(.) {
    s <- sample.int(length(grid$x), replace = TRUE)
    b1 <- cor(x[s], y[s]) * sd(y[s]) / sd(x[s])
    b0 <- mean(y[s]) - b1 * mean(x[s])
    b1 * x + b0
  }
)
X <- cbind(rep(1, length(x)), x)
set.seed(0)
p2 <- sapply(
  seq_len(10),
  function(.) {
    samp <- sample.int(length(grid$x), replace = TRUE)
    Xsamp <- X[samp, ]
    X %*% solve(t(Xsamp) %*% Xsamp, t(Xsamp) %*% y[samp])
  }
)
all.equal(c(p2), c(p1))

## now test:

mbm <- microbenchmark(
  
  "cor" = {
    
    sapply(
      seq_len(1000),
      function(.) {
        s <- sample.int(length(grid$x), replace = TRUE)
        b1 <- cor(x[s], y[s]) * sd(y[s]) / sd(x[s])
        b0 <- mean(y[s]) - b1 * mean(x[s])
        b1 * x + b0
      }
    )
    
  },
  "lin.sys" = {
    
    sapply(
      seq_len(1000),
      function(.) {
        samp <- sample.int(length(grid$x), replace = TRUE)
        Xsamp <- X[samp, ]
        X %*% solve(t(Xsamp) %*% Xsamp, t(Xsamp) %*% y[samp])
      }
    )
    
    }
)

mbm
autoplot(mbm)