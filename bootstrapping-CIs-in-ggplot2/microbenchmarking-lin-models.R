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

## and test within full ggplot function:

BootCI <- ggplot2::ggproto(
  "BootCI", ggplot2::Stat, 
  required_aes = c("x", "y"),
  compute_group = function(data, scales, params, n = 1000, percent = 95, lmfun) {
    
    grid <- data.frame(x = data$x)
    
    X <- cbind(rep(1, length(x)), x)  ## design matrix (includes intercept)
    
    if (lmfun == "lin.sys") {
      
      predictions <- sapply(
        seq_len(n),
        FUN = function(.) {
          samp <- sample.int(length(grid$x), replace = TRUE)
          Xsamp <- X[samp, ]
          X %*% solve(t(Xsamp) %*% Xsamp, t(Xsamp) %*% y[samp])  ## get bs then dot with X
        }
      )
      
    } else if (lmfun == "cor") {
      
      predictions <- sapply(
        seq_len(n),
        function(.) {
          s <- sample.int(length(grid$x), replace = TRUE)
          b1 <- cor(x[s], y[s]) * sd(y[s]) / sd(x[s])
          b0 <- mean(y[s]) - b1 * mean(x[s])
          b1 * x + b0
        }
      )
      
    } else if (lmfun == "lm") {
      
      predictions <- sapply(
        seq_len(n),
        function(.)
          predict(
            lm(y ~ x, data[sample.int(nrow(data), replace = TRUE), ]),
            grid,
          )
      )
      
    }
    
    .alpha <- (100 - percent) / 200  ## 2 tailed
    grid$ymax <- apply(predictions, 1, quantile, 1 - .alpha)
    grid$ymin <- apply(predictions, 1, quantile, .alpha)
    
    grid
    
  }
  
)
stat_boot_ci <- function(mapping = NULL, data = NULL, geom = "ribbon",
                         position = "identity", na.rm = FALSE, show.legend = NA, 
                         inherit.aes = TRUE, n = 1000, percent = 95, ...) {
  ## see: https://cran.r-project.org/web/packages/ggplot2/vignettes/extending-ggplot2.html
  ggplot2::layer(
    stat = BootCI, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(n = n, percent = percent, na.rm = na.rm, ...)
  )
}



set.seed(0)
# n <- 1000  ## n datapoints
x <- 1:n  ## some measure x
y <- rnorm(x, x, sd = x)  ## some measure y; mean and variance depends on x
p <- ggplot(data.frame(x, y), aes(x, y), environment = environment()) +
  geom_point(fill = "black", color = "white", shape = 21, size = 3) +
  theme(
    axis.ticks = element_blank(),
    axis.text  = element_blank(),
    axis.line  = element_blank(),
    panel.background = element_blank()
  )

mbm <- microbenchmark(
  
  "lin.sys" = {
    
    p <- p +
      stat_boot_ci(
        alpha = 0.3, n = 1000, percent = 95,
        lmfun = "lin.sys"
      )
    print(p)
    
  },
  
  "cor" = {
    
    p <- p +
      stat_boot_ci(
        alpha = 0.3, n = 1000, percent = 95,
        lmfun = "cor"
      )
    print(p)
    
  },
  
  "lm" = {
    
    p <- p +
      stat_boot_ci(
        alpha = 0.3, n = 1000, percent = 95,
        lmfun = "lm"
      )
    print(p)
    
  }
  
)

mbm
autoplot(mbm)

walltime <- function(FUN) {
  starttime <- Sys.time()
  (out <- FUN)
  endtime <- Sys.time()
  endtime - starttime
}
walltime(
  p +
    stat_boot_ci(
      alpha = 0.3, n = 1E4, percent = 95,
      lmfun = "lm"
    )
)
walltime(f
  p +
    geom_smooth(method = "lm", se = TRUE, fill = "red", alpha = 1) +
    stat_boot_ci(
      alpha = 0.3, n = 1E4, percent = 95,
      lmfun = "cor"
    )
)

p +
  stat_boot_ci(
    alpha = 0.3, n = 1E4, percent = 95,
    lmfun = "lm"
  )