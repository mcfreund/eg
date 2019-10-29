## given set of regressors, fit all unique combinations and store results in list
## requires gtools for combinations() function
 
## create sample data ----
df <- as.data.frame(matrix(rnorm(1000), ncol = 10))
cols <- names(df)

## model comparison ----

lm.allcombs <- function(df, yname) {
  ## dependencies: gtools()
  ## warning: number of models fit increases exponentially with number of explanatory variables

  ## model comparison ----
  
  xnames <- names(df)[-grep(yname, names(df))]
  nvars  <- length(xnames)  ## num explanatory vars
  combs  <- lapply(1:nvars, function(x) gtools::combinations(n = nvars, r = x))  ## unique combos for each number of vars
  nmods  <- sum(sapply(combs, nrow))  ## total num models
  lmods  <- vector("list", nmods)  ## output list
  
  a <- 1
  for (ii in seq_along(combs)) {  ## along number of vars
    
    combs.ii <- combs[[ii]]
    
    for (jj in seq_len(nrow(combs.ii))) {  ## along unique combos
      
      xnames.jj <- xnames[combs.ii[jj, ]]
      
      lmform <- as.formula(paste0(yname, " ~ ", paste0(xnames.jj, collapse = " + ")))
      
      lmods[[a]] <- lm(lmform, df[, c(xnames.jj, yname)])
      
      a <- a + 1
      
    }
    
  }
  
  lmods
  
}


## use ----

l <- lm.allcombs(df, "V8")  ## to use

bics <- sapply(l, BIC)
aics <- sapply(l, AIC)

which.min(aics)
which.min(bics)

order(aics)
order(bics)

plot(bics, aics)  ## how great is spread between best and next-best models?

lmods[order(bics)[1:4]]  ## best 4
lmods[order(aics)[1:4]]  ## best 4

fit.best <- l[[which.min(bics)]]
fit.best.coef <- coef(summary(fit.best))
cols.best <- names(coef(fit.best))[-1]

