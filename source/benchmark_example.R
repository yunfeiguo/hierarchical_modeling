library('microbenchmark')
f <- function(n, x) for (i in 1:n) x = (1 + x)^(-1)
g <- function(n, x) for (i in 1:n) x = (1 + x)^(-1)

compare <- microbenchmark(f(1000,1),g(1000,1),times=1000)
compare
